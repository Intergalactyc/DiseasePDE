import os
from lib.paths import DATA_PRODUCTS, COUNTIES_SHP, STATES_SHP, FIPS_CSV

COMPARTMENTS = {"I" : "infected", "R" : "recovered", "S" : "susceptible"}
OVERLAPS_PATH = os.path.join(DATA_PRODUCTS, "intermediate", "county_overlaps")

#def get_counties
