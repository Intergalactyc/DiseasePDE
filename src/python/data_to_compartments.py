from lib.paths import DATA_PRODUCTS, CUMULATIVE_CSV, POPULATION_CSV
import argparse
import os
import pandas as pd
from tqdm import tqdm

parquet_path = os.path.join(DATA_PRODUCTS, "intermediate", "cumulatives_clean.parquet")

def get_cumulative_timeseries() -> pd.DataFrame:
    cumulative = pd.read_csv(CUMULATIVE_CSV)

    for i in cumulative.index: 
        admin2 = cumulative.loc[i, "Admin2"]
        if pd.isna(admin2) or "out of" in admin2.lower() or "unassigned" in admin2.lower():
            cumulative.drop(i, inplace = True)

    timeseries = cumulative.dropna(subset = ["FIPS"]).drop(columns = ["UID","iso2","iso3","code3","Admin2","Province_State","Country_Region","Lat","Long_","Combined_Key"]).transpose()
    timeseries.rename(columns = timeseries.iloc[0].astype(int).astype(str), inplace = True)
    timeseries = timeseries[1:]
    timeseries.index = pd.to_datetime(timeseries.index, format = "%m/%d/%y")
    timeseries = pd.concat([pd.DataFrame(data = 0., columns = timeseries.columns, index = [timeseries.index[0] - pd.Timedelta("1d")]), timeseries])

    return timeseries

def clean_decreases(cumulatives: pd.DataFrame) -> pd.DataFrame:
    # Backward fill method to keep cumulative cases monotonically decreasing
    for col in tqdm(cumulatives.columns):
        s = cumulatives[col]
        for i in range(len(s)-1, 0, -1):
            if s.iloc[i] < s.iloc[i-1]:
                s.iloc[i-1] = s.iloc[i]

    return cumulatives

def smooth_input(cumulatives: pd.DataFrame, smoothing: int = 7) -> pd.DataFrame:
    return cumulatives.rolling(window = f"{smoothing}h").mean() # TODO: this might not be the best way to do this

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

def get_population_total(fips: pd.Index, index: pd.Index) -> pd.DataFrame:
    populations = pd.read_csv(POPULATION_CSV, 
                              skiprows = [1]) # skip metadata row

    populations["FIPS"] = populations["GEO_ID"].apply(lambda x : str(int(x[-5:])))
    populations = populations[populations["FIPS"].isin(fips)]
    populations = populations[["FIPS","DP05_0001E"]].set_index("FIPS").rename(columns = {"DP05_0001E" : "population"})

    return pd.DataFrame(data = [populations["population"] for _ in index], columns = populations.index, index = index)

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
        print("Recomputing cumulative timeseries.")
        cumulative_ts = get_cumulative_timeseries()
        print("Cleaning up case correction decreases using backward fill...")
        cumulative_ts = clean_decreases(cumulative_ts)
        if not args["no_input_smoothing"]:
            cumulative_ts = smooth_input(cumulative_ts, smoothing = 24 * 7) # smooth input by 1 week
        cumulative_ts.to_parquet(parquet_path)
        print("Saved to file.")
    else:
        print("Loading cumulative cases from file.")
        cumulative_ts = pd.read_parquet(parquet_path)
        print("Cases loaded.")

    assert cumulative_ts.diff().lt(0).sum().sum() == 0, "Uncorrected decrease detected in cumulative timeseries"  

    print(f"Solving for I/R compartments using recovery period of {args['recovery_period']:.1f} days with output resolution of {args['output_resolution']} hrs.")

    infected, recovered = solve_for_IR(cumulatives = cumulative_ts,
                                       recovery_rate = 1 / (args["recovery_period"] * 24),
                                       resolution_hours = args["output_resolution"],
                                       smoothing = 24 * 7 if not args["no_output_smoothing"] else 0) # smooth by 1 week (unless disabled)

    assert infected.lt(0).sum().sum() == 0, "Negative infected case count detected"
    assert recovered.lt(0).sum().sum() == 0, "Negative recovery count detected"

    print("Loading total population data to determine S compartment.")

    total_population = get_population_total(fips = infected.columns, index = infected.index)

    notin = [f for f in infected.columns if f not in total_population.columns]

    infected.drop(columns = notin, inplace = True)
    recovered.drop(columns = notin, inplace = True)

    susceptible = total_population - infected - recovered

    print("Saving results.")

    susceptible.to_csv(os.path.join(DATA_PRODUCTS, "compartments", "susceptible.csv"))
    infected.to_csv(os.path.join(DATA_PRODUCTS, "compartments", "infected.csv"))
    recovered.to_csv(os.path.join(DATA_PRODUCTS, "compartments", "recovered.csv"))

    print("Saved S/I/R compartments to CSV files.")

if __name__ == "__main__":
    main()
