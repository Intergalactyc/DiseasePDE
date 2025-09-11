from lib.paths import DATA_PRODUCTS, STATES_SHP, LAKES_SHP, FIPS_CSV
import argparse
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
import lib.polygons as polygons
import os

fipsdf = pd.read_csv(FIPS_CSV, skipinitialspace=True)
consecutive_fips = list(fipsdf.loc[~fipsdf["stname"].isin(["Alaska", "Hawaii"]), "st"])

def load_all_states():
    states = gpd.read_file(STATES_SHP)
    states["COUNTRY"] = states.apply(lambda row : str(row["fips"])[:2], axis = 1)
    states = states[states["COUNTRY"] == "US"]
    states["FIPS"] = states.apply(lambda row : int(str(row["fips"])[2:]), axis = 1)
    return states

def load_lakes(scale):
    lakes = gpd.read_file(LAKES_SHP)
    lakes = lakes[lakes['scalerank'].astype(int) <= scale] # consider only the lakes at ScaleRank <= lake_scale
    return lakes.union_all()

def get_polygon(selected, lakes, tolerance, name):
    merged = selected.union_all()
    merged = merged.difference(lakes) if lakes else merged
    if isinstance(merged, Polygon):
        primary = merged
    else:
        primary = None
        for polygon in merged.geoms:
            if primary is None or polygon.area > primary.area:
                primary = polygon
    simplified = primary.simplify(tolerance = tolerance, preserve_topology = True) # Douglas-Peucker simplification
    polygons.save_polygon_as_poly(simplified, os.path.join(DATA_PRODUCTS, "polygons", f"{name.lower().replace(' ', '_')}.poly"), silent = True)

def get_us_polygon(states, lakes, tolerance):
    selected = states[states["FIPS"].isin(consecutive_fips)]
    get_polygon(selected, lakes, tolerance, "us")

def get_state_polygons(states, lakes, tolerance):
    failed = []
    for _, state in fipsdf.iterrows():
        name = state["stname"]
        fips = state["st"]
        try:
            selected = states[states["FIPS"] == fips]
            get_polygon(selected, lakes, tolerance, name)
        except Exception:
            failed.append(name)
    if len(failed) > 0:
        print("Failed to get polygons for the following states: " + ", ".join(failed))

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--us-tolerance", type = float, default = 0.01)
    parser.add_argument("--state-tolerance", type = float, default = 0.005)
    parser.add_argument("--us-lake-scale", type = int, default = 2)
    parser.add_argument("--state-lake-scale", type = int, default = 3) # higher value -> smaller scale
    parser.add_argument("--no-lakes", action = "store_true", default = False)

    return vars(parser.parse_args())

def main():
    args = parse()

    states = load_all_states()
    
    if args["no_lakes"]:
        lakes_us, lakes_states = None, None
    else:
        lakes_us = load_lakes(args["us_lake_scale"])
        lakes_states = lakes_us if args["us_lake_scale"] == args["state_lake_scale"] else load_lakes(args["state_lake_scale"])

    print("Getting US polygon...")
    get_us_polygon(states, lakes_us, args["us_tolerance"])

    print("Getting individual state polygons...")
    get_state_polygons(states, lakes_states, args["state_tolerance"])

    print("Complete!")

if __name__ == "__main__":
    main()
