import os
from lib.paths import DATA_PRODUCTS, COUNTIES_SHP, STATES_SHP, FIPS_CSV
import geopandas as gpd
import pandas as pd
from shapely import Polygon, prepare
import meshio
from tqdm import tqdm
import json
import numpy as np
from multiprocessing import Pool

COMPARTMENTS = {"I" : "infected", "R" : "recovered", "S" : "susceptible"}
OVERLAPS_PATH = os.path.join(DATA_PRODUCTS, "intermediate", "county_overlaps.json")
INPUT_MESH_PATH = os.path.join(DATA_PRODUCTS, "meshes", "us.exo")
OUTPUT_MESH_PATH = os.path.join(DATA_PRODUCTS, "datameshes")

LIMIT = [100,101]
NPROC = 12

def get_counties():
    # get FIPS codes of the 48 consecutive states (+DC)
    statefp = pd.read_csv(FIPS_CSV, skipinitialspace=True)
    consecutive_fips = list(statefp.loc[~statefp["stname"].isin(["Alaska", "Hawaii"]), "st"])
    # import counties and select the ones in the consecutive states
    counties = gpd.read_file(COUNTIES_SHP)
    counties["FIPS"] = counties["CODE_LOCAL"].astype(int)
    counties = counties[["FIPS", "REGION_COD", "AREA_SQKM", "geometry"]].rename(columns = {"AREA_SQKM" : "area"})
    counties["REGION_COD"] = counties["REGION_COD"].astype(int)
    counties = counties[counties["REGION_COD"].isin(consecutive_fips)]
    counties["geometry"].apply(lambda poly : prepare(poly)) # Prepare every polygon in geometries to make later checks faster. Maybe not working right? (No return values, supposed to be in-place)
    # counties["area_alt"] = counties["geometry"].apply(lambda poly : poly.area)
    # print(counties[["area","area_alt"]])
    return counties

def get_states():
    # get FIPS codes of the 48 consecutive states (+DC)
    statefp = pd.read_csv(FIPS_CSV, skipinitialspace=True)
    consecutive_fips = list(statefp.loc[~statefp["stname"].isin(["Alaska", "Hawaii"]), "st"])
    # import and select US consecutive states
    states = gpd.read_file(STATES_SHP)
    states["COUNTRY"] = states.apply(lambda row : str(row["fips"])[:2], axis = 1)
    states = states[states["COUNTRY"] == "US"]
    states["FIPS"] = states.apply(lambda row : int(str(row["fips"])[2:]), axis = 1)
    states = states[states["FIPS"].isin(consecutive_fips)]
    states = states[["geometry", "FIPS"]]
    return states

def get_data(names):
    result = dict()
    for shortname, filename in names.items():
        filepath = os.path.join(DATA_PRODUCTS, "compartments", f"{filename}.csv")
        result[shortname] = pd.read_csv(filepath)
    return result

def elements_from_mesh(mesh) -> list[Polygon]:
    return [Polygon([mesh.points[i] for i in triangle]) for triangle in tqdm(mesh.cells_dict["triangle"])]

def _interpolator_worker(args):
    interpolator, elem_counties, i = args
    return {
        k : np.array([interpolator._interpolate_element(cs, k, i)
            for cs in elem_counties])
        for k in interpolator.data.keys()
    }

class CountyInterpolator:
    def __init__(self, counties: gpd.GeoDataFrame, states: gpd.GeoDataFrame, data: dict, index: pd.Index|str = "infer"):
        self.counties = counties
        self.states = states
        self.state_county_map = {statefips : counties[counties["REGION_COD"].astype(int) == int(statefips)] for statefips in self.states["FIPS"]}
        self.county_areas = {fips : area for fips, area in zip(counties["FIPS"], counties["area"])}
        self.data = data
        if index == "infer": # The index is the time index
            ind = None
            for df in self.data.values():
                if ind is None:
                    ind = df.index
                else:
                    ind = ind.intersection(df.index)
            self.index = ind
        elif type(index) is pd.Index:
            self.index = index
        else:
            raise Exception("Specify valid index or leave to 'infer'")

    def _states_overlapping(self, element: Polygon) -> list:
        # Returns list of states(" FIPS codes) which given element intersects
        result = []
        for i in range(len(self.states)):
            s = self.states.iloc[i]
            if s["geometry"].intersects(element):
                result.append(s["FIPS"])
        # Polygon POLYGON ((-80.5996292794999 28.59437948250004, -80.59884599499989 28.588471788000042, -80.59721054540965 28.591694401622195, -80.5996292794999 28.59437948250004)) never matches
        return result

    def _counties_overlapping(self, element: Polygon) -> dict:
        # Returns dict of counties intersected as `FIPS : fraction overlap` pairs
        raw_result = dict()
        result = dict()
        states_fips = self._states_overlapping(element)
        # Technically the weighting is slightly inaccurate because we are in coords of lat-long
            # However, each element is so small that it should be basically correct
        sum_weights = 0.
        for sf in states_fips:
            counties = self.state_county_map[sf]
            for i in range(len(counties)):
                c = counties.iloc[i]
                if c["geometry"].intersects(element):
                    weight = c["geometry"].intersection(element).area
                    raw_result[int(c["FIPS"])] = weight
                    sum_weights += weight
        result = {cf : w / sum_weights for cf, w in raw_result.items()} # So that total weight is 1
        if len(result.keys()) == 0:
            print(f"Failed to find any counties matching {element}")
        return result
    
    def _interpolate_element(self, counties, compartment, index) -> float:
        # At given index, compute the interpolation of given data compartment onto given element `elem`
            # Requires passing known county overlap dictionary (previously obtained from self._counties_overlapping()) for elem
        result = 0.
        for county, weight in counties.items():
            result += weight * data[compartment][county].iloc[index] / self.county_areas[int(county)]
        return result
    
    def _mproc_interpolation_result(self, elem_counties, index, nproc: int):
        if type(nproc) is not int or nproc < 1:
            raise Exception(f"Invalid number of processors {nproc=}")
        elif nproc == 1:
            result = [ # List in which each element is a timestep
                { # Each timestep is a dictionary of `compartment : [data value by element]` pairs
                    k : np.array([self._interpolate_element(cs, k, i)
                        for cs in elem_counties]) # Data list follows same order as elems
                    for k in self.data.keys()
                }
                for i in tqdm(index)
            ]
        else:
            result = []
            taskqueue = [(self, elem_counties, i) for i in index]
            pbar = tqdm(total = len(taskqueue))
            pool = Pool(processes = nproc)
            for res in pool.imap(_interpolator_worker, taskqueue):
                result.append(res)
                pbar.update()
                pbar.refresh()
            pool.close()
            pool.join()

        return result
    
    def interpolate_onto_mesh(self, mesh, saveto: str = None, county_file: str = None, index_limit: tuple[int,int] = None, nproc: int = 1) -> list[dict]:        
        loaded = False
        if county_file is not None and os.path.exists(county_file): # Attempt to load county overlaps
            print("Loading county overlap data...")
            try:
                with open(county_file, "r") as file:
                    elem_counties = json.load(file)
                loaded = True
                print("Loaded successfully.")
            except Exception as e:
                print(f"Error encountered while loading: {e}")
        if not loaded: # If we couldn't load county overlaps, recompute them
            print("Getting elements from mesh...")
            elems = elements_from_mesh(mesh)
            print("Computing county overlaps...")
            elem_counties = [self._counties_overlapping(e) for e in tqdm(elems)]
            if county_file:
                print("Saving county overlap data...")
                with open(county_file, "w") as file:
                    json.dump(elem_counties, file)

        print("Interpolating...")
        limited_index = self.index if index_limit is None else self.index[index_limit[0]:index_limit[1]]
        result = self._mproc_interpolation_result(elem_counties, limited_index, nproc)

        if saveto is not None:
            print("Saving outputs...")
            meshio.write(
                os.path.join(saveto, "us_data.exo"),
                mesh,
                cell_data = result,
                file_format = "exodus"
            )
            for i, cell_data in zip(limited_index, result):
                filename = os.path.join(saveto, f"SIR_{i}.exo")
                meshio.write(
                    filename,
                    mesh,
                    timestep_arr = [i],
                    cell_data = cell_data,
                    file_format = "exodus"
                )
            print(f"Saved in {saveto}")
            return
        else:
            return result

if __name__ == "__main__":
    print("Loading data...")
    counties = get_counties()
    states = get_states()
    data = get_data(names = COMPARTMENTS)
    interpolator = CountyInterpolator(counties, states, data)

    mesh = meshio.read(INPUT_MESH_PATH)

    # Use index_limit to go through fewer timesteps
    interpolator.interpolate_onto_mesh(mesh, saveto = OUTPUT_MESH_PATH, county_file = OVERLAPS_PATH, index_limit = LIMIT, nproc = NPROC)
    print("Complete!")
