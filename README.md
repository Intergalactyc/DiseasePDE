## Data sources
*Names in parentheses are the corresponding options in the `[paths]` section of the config.ini file*

COVID-19 confirmed cases (`cumulative_cases_timeseries_csv`): [Johns Hopkins CSSE (GitHub)](https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv)

## Data transformation
Copy the `config.ini.template`, rename it to `config.ini`, and change the paths within it according to the locations of the data on your computer. Ensure that the ini file remains in the top level of this repository directory.

Install the required Python packages into your environment (a virtual environment is recommended) using `python -m pip install -r requirements.txt`.

Within `src/python` run these files in the following order:
1. `get_polygons.py` and `data_to_compartments.py` (order of these doesn't matter)
    - Data for all counties will be used, and polygons will be generated for every state as well as for the full contiguous US.
2. `mesh_polygon.py`
    - By default, the US polygon will be meshed. To instead mesh a state polygon, pass the argument `--state` followed by the name of the state.
3. `interpolate_over_mesh.py`
    - Right now this only works over the US polygon generated from `mesh_polygon.py` without a `--state` argument.

Final result is an `.exo` file located in `data_products/datameshes`.
