import configparser
import pathlib
from warnings import warn
import os

PARENT = pathlib.Path(__file__).parent.parent.parent.parent.resolve()
DATA_PRODUCTS = os.path.join(PARENT, "data_products")
CONFIG_FILE = os.path.join(PARENT, "config.ini")

if not os.path.exists(CONFIG_FILE):
    raise OSError(f"Config file not found (expected {CONFIG_FILE}) - copy, update, and rename the config.ini.template file")

for product in ["compartments", "datameshes", "intermediate", "meshes", "polygons"]:
    dirpath = os.path.join(DATA_PRODUCTS, product)
    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)

parser = configparser.ConfigParser()
parser.read(CONFIG_FILE)

def get_item(section, option, default = None, cast = None):
    try:
        res = parser.get(section, option)
        if cast:
            return cast(res)
        return res
    except (KeyError, ValueError, configparser.NoSectionError, configparser.NoOptionError) as e:
        warn(f"Could not completely parse configuration file: {e}")
        return default

COUNTIES_SHP = get_item("paths", "counties_shapefile")
STATES_SHP = get_item("paths", "states_shapefile")
LAKES_SHP = get_item("paths", "lakes_shapefile")
FIPS_CSV = get_item("paths", "fips_csv")
POPULATION_CSV = get_item("paths", "county_population_csv")
CUMULATIVE_CSV = get_item("paths", "cumulative_cases_timeseries_csv")
NPROC = get_item("system", "nproc", cast = int)
