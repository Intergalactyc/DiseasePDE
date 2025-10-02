import pathlib
import os

PARENT = pathlib.Path(__file__).parent.parent.parent.parent.resolve()
DATA_SOURCES = os.path.join(PARENT, "data_sources")
DATA_PRODUCTS = os.path.join(PARENT, "data_products")

for product in ["compartments", "datameshes", "intermediate", "meshes", "polygons"]:
    dirpath = os.path.join(DATA_PRODUCTS, product)
    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)

COUNTIES_SHP = os.path.join(DATA_SOURCES, "natural_earth", "us_counties_10m_no-large-lakes", "ne_10m_admin_2_counties_lakes.shp")
STATES_SHP = os.path.join(DATA_SOURCES, "natural_earth", "states_provinces_10m_no-large-lakes", "ne_10m_admin_1_states_provinces_lakes.shp")
LAKES_SHP = os.path.join(DATA_SOURCES, "natural_earth", "ne_10m_lakes", "ne_10m_lakes.shp")
FIPS_CSV = os.path.join(DATA_SOURCES, "census", "state_fips.csv")
POPULATION_CSV = os.path.join(DATA_SOURCES, "census", "2020CensusDP05CountyLevel", "ACSDP5Y2020.DP05-Data.csv")
CUMULATIVE_CSV = os.path.join(DATA_SOURCES, "jhu", "time_series_covid19_confirmed_US.csv")
ADJACENCIES_FILE = os.path.join(DATA_SOURCES, "census", "county_adjacency2023.txt")
