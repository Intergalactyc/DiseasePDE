from lib.paths import DATA_PRODUCTS, CUMULATIVE_CSV, POPULATION_CSV, ADJACENCIES_FILE
import argparse
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import logging

LOGFILE = "dtc_debug.log"
with open(LOGFILE, "w") as f:
    pass
logging.basicConfig(filename=LOGFILE, level=logging.DEBUG)
logger = logging.getLogger(__name__)

parquet_path = os.path.join(DATA_PRODUCTS, "intermediate", "cumulatives_clean.parquet")

DEPARTMENTS = {
    84070015 : {49003, 49005, 49033}, # Bear River (UT)
    84070016 : {49023, 49027, 49039, 49041, 49031, 49055}, # Central Utah (UT)
    84070017 : {49007, 49015, 49019}, # Southeast Utah (UT)
    84070018 : {49001, 49017, 49021, 49025, 49053}, # Southwest Utah (UT)
    84070019 : {49009, 49013, 49047}, # Tricounty (UT)
    84070020 : {49057, 49029}, # Weber-Morgan (UT)
    84070002 : {25007, 25019}, # Dukes & Nantucket (MA)
    84070003 : {29095}, # Kansas City (MO)
# Missing Michigan Department of Corrections & Federal Correctional Institution (both MI), as well as all Unassigned
} # Map dept UID to set of FIPS codes for counties represented
ISLAND_COUNTIES = {
    25007, # Dukes, MA
    25019, # Nantucket, MA
    36047, # Kings, NY
    36059, # Nassau, NY
    36061, # New York, NY
    36081, # Queens, NY
    36085, # Richmond, NY
    36103, # Suffolk, NY
    53029, # Island, WA
    53055, # San Juan, WA
} # Set of FIPS codes for island counties

adjs = pd.read_csv(ADJACENCIES_FILE, delimiter = "|")
adjs = adjs[~adjs["Neighbor Name"].str.contains("Planning Region")] # Can't deal with CT planning regions
adjs = adjs[["County GEOID", "Neighbor GEOID"]]
def _get_neighbors(county: int, mainland: bool=True) -> set[int]:
    all_neighbors = {n for n in adjs.loc[adjs["County GEOID"] == county, "Neighbor GEOID"] if n != county}
    if mainland:
        result = {n for n in all_neighbors if n not in ISLAND_COUNTIES}
        if len(result) == 0:
            return {n for c in all_neighbors for n in _get_neighbors(c, False)}
        return result
    return all_neighbors

def get_cumulative_timeseries() -> pd.DataFrame:
    cumulative = pd.read_csv(CUMULATIVE_CSV)

    for i in cumulative.index: 
        admin2 = cumulative.loc[i, "Admin2"]
        if pd.isna(admin2) or "out of" in admin2.lower() or "unassigned" in admin2.lower():
            cumulative.drop(i, inplace = True)

    for index, row in cumulative.iterrows():
        if pd.isna(row["FIPS"]):
            cumulative.loc[index, "FIPS"] = row["UID"]
    timeseries = cumulative.dropna(subset = ["FIPS"]).drop(columns = ["UID","iso2","iso3","code3","Admin2","Province_State","Country_Region","Lat","Long_","Combined_Key"]).transpose()
    timeseries.rename(columns = timeseries.iloc[0].astype(int).astype(str), inplace = True)
    timeseries = timeseries[1:]
    timeseries.index = pd.to_datetime(timeseries.index, format = "%m/%d/%y")
    timeseries = pd.concat([pd.DataFrame(data = 0., columns = timeseries.columns, index = [timeseries.index[0] - pd.Timedelta("1d")]), timeseries])

    return timeseries

def clean_decreases(cumulatives: pd.DataFrame) -> pd.DataFrame:
    # Backward fill method to keep cumulative cases monotonically nondecreasing
    for col in tqdm(cumulatives.columns):
        s = cumulatives[col]
        for i in range(len(s)-1, 0, -1):
            if s.iloc[i] < s.iloc[i-1]:
                s.iloc[i-1] = s.iloc[i]

    return cumulatives

def smooth_batches(cumulatives: pd.DataFrame, min_flat_days: int = 2, jump_window: int = 7, only_after: int = 0, jump_factor: float = 1) -> pd.DataFrame:
    # Detect batch (e.g. weekly, semi-weekly) reporting and linearly interpolate reports within the batch period
    result = cumulatives.copy()
    
    for county in tqdm(result.columns):
        logger.debug(f"County {county}")
        series = result[county].astype(float)

        has_jumped = False
        i = only_after
        while i < len(series) - 1:
            if series.iloc[i+1] - series.iloc[i] == 0:
                start = i
                while i < len(series) - 1 and series.iloc[i+1] - series.iloc[i] == 0:
                    i += 1
                flat_days = i - start + 1

                if i >= len(series):
                    break

                if not has_jumped:
                    has_jumped = True
                    i += 1
                    logger.debug("Hasn't jumped yet")
                    continue

                if flat_days >= min_flat_days:
                    logger.debug(f" Batched report found of length {flat_days}")
                    # Batched report found, smooth it
                    start_val = series.iloc[start - 1] if start > 0 else series.iloc[start]
                    end_val = series.iloc[i]
                    interp_vals = np.linspace(start_val, end_val, flat_days + 1)[1:-1]
                    series.iloc[start:i] = interp_vals
            else:
                i += 1

        result[county] = series

    return result

def solve_for_IR(cumulatives: pd.DataFrame,
                 recovery_rate: float,
                 resolution_hours: int,
                 smoothing: int = 0)-> tuple[pd.DataFrame, pd.DataFrame]:  
    delta_cases = cumulatives.diff()

    # 0-initialize infected and recovered DataFrames for results
    infected = pd.DataFrame(data = 0., index = pd.date_range(start = cumulatives.index[1], end = cumulatives.index[-1], freq = f"{resolution_hours}h"), columns = cumulatives.columns)
    recovered = infected.copy()

    # Useful constants for ITR timestepping
    _alpha = 2 / (2 + recovery_rate * resolution_hours)
    _beta0 = (recovery_rate * resolution_hours)/2
    _beta = 1 - _beta0

    # ITR timestepping to solve I' = C' - (recovery_rate)*I
    for i in tqdm(range(1, len(infected.index))):
        current_index = int((i * resolution_hours / 24) // 1) + 1
        DC_DT = delta_cases.iloc[current_index] / 24 # Derivative of cumulative cases taken as rate of change between this and next timestep, no need to average
        infected.iloc[i] = _alpha * (infected.iloc[i-1] * _beta + resolution_hours * DC_DT)
        recovered.iloc[i] = recovered.iloc[i-1] + _beta0 * (infected.iloc[i-1] + infected.iloc[i])

    if smoothing > 0:
        infected = infected.rolling(window = f"{smoothing}h").mean()
        recovered = recovered.rolling(window = f"{smoothing}h").mean()

    return infected, recovered

def get_population_total_and_distribute(cumulatives: pd.DataFrame, fips: pd.Index, index: pd.Index) -> tuple[pd.DataFrame, pd.DataFrame]:
    cum_result = cumulatives.copy()
    populations = pd.read_csv(POPULATION_CSV, 
                              skiprows = [1]) # skip metadata row

    populations["FIPS"] = populations["GEO_ID"].apply(lambda x : str(int(x[-5:])))
    populations = populations[populations["FIPS"].isin(fips)]
    populations = populations[["FIPS","DP05_0001E"]].set_index("FIPS").rename(columns = {"DP05_0001E" : "population"})

    populations_redistributed = populations.copy().astype(float)
    for dept, counties in DEPARTMENTS.items():
        dept_data = cum_result[str(dept)]
        weights = {}
        total_weight = 0
        for c in counties:
            val = populations.loc[str(c), "population"]
            weights[c] = val
            total_weight += val
        weights = {c : w/total_weight for c, w in weights.items()}
        for c in counties:
            cum_result[str(c)] += dept_data * weights[c]

    for county in ISLAND_COUNTIES:
        pop_to_distribute = populations_redistributed.loc[str(county), "population"]
        cases_to_distribute = cum_result[str(county)]
        neighbors = _get_neighbors(county)
        weights = {}
        total_weight = 0
        for c in neighbors:
            val = populations.loc[str(c), "population"]
            weights[c] = val
            total_weight += val
        weights = {c : w/total_weight for c, w in weights.items()}
        for c, w in weights.items():
            populations_redistributed.loc[str(c), "population"] += w * pop_to_distribute
            cum_result[str(c)] += w * cases_to_distribute
    
    pop_result = pd.DataFrame(data = [populations_redistributed["population"] for _ in index], columns = populations_redistributed.index, index = index)

    return cum_result, pop_result

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reload", action = "store_true", help = "Force reloading + formatting of cumulative data")
    parser.add_argument("--recovery-period", type = float, default = 8, help = "Assumed recovery period (days)")
    parser.add_argument("--no-input-smoothing", action = "store_true", default = False, help = "Disable input smoothing (warning - this will result in oscillation artifacts)")
    parser.add_argument("--output-resolution", type = int, default = 24, help = "Time resolution of output (hours)")
    parser.add_argument("--no-output-smoothing", action = "store_true", default = False, help = "Disable output smoothing")
    return vars(parser.parse_args())

def main():
    args = parse()
    if args["reload"] or not os.path.exists(parquet_path):
        print("Recomputing cumulative timeseries...")
        cumulative_ts = get_cumulative_timeseries()
        print("Cleaning up case correction decreases using backward fill...")
        cumulative_ts = clean_decreases(cumulative_ts)
        if not args["no_input_smoothing"]:
            print("Smoothing batch reports...")
            cumulative_ts = smooth_batches(cumulative_ts)
        cumulative_ts.to_parquet(parquet_path)
        print("Saved to file.")
    else:
        print("Loading cumulative cases from file.")
        cumulative_ts = pd.read_parquet(parquet_path)
        print("Cases loaded.")

    assert cumulative_ts.diff().lt(0).sum().sum() == 0, "Uncorrected decrease detected in cumulative timeseries"  

    print("Loading total population data, distributing departments, and solving the island problem.")
    cumulative_ts, total_population = get_population_total_and_distribute(cumulatives=cumulative_ts, fips = cumulative_ts.columns, index = cumulative_ts.index)

    print(f"Solving for I/R compartments using recovery period of {args['recovery_period']:.1f} days with output resolution of {args['output_resolution']} hrs.")

    infected, recovered = solve_for_IR(cumulatives = cumulative_ts,
                                       recovery_rate = 1 / (args["recovery_period"] * 24),
                                       resolution_hours = args["output_resolution"],
                                       smoothing = 24 * 7 if not args["no_output_smoothing"] else 0) # smooth by 1 week (unless disabled)

    assert infected.lt(0).sum().sum() == 0, "Negative infected case count detected"
    assert recovered.lt(0).sum().sum() == 0, "Negative recovery count detected"

    print("Using total population data to compute S compartment.")
    notin = [f for f in infected.columns if f not in total_population.columns]
    infected.drop(columns = notin, inplace = True)
    recovered.drop(columns = notin, inplace = True)
    susceptible = total_population - infected - recovered # If there is an issue with mismatch, may need to 

    print("Saving results.")

    susceptible.to_csv(os.path.join(DATA_PRODUCTS, "compartments", "susceptible.csv"))
    infected.to_csv(os.path.join(DATA_PRODUCTS, "compartments", "infected.csv"))
    recovered.to_csv(os.path.join(DATA_PRODUCTS, "compartments", "recovered.csv"))

    print("Saved S/I/R compartments to CSV files.")

if __name__ == "__main__":
    main()
