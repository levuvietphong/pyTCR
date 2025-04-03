(Examples)=

# Examples
`pyTCR` generates and provides access to a large set of synthetic rainfall events based on TC tracks. 

The repository includes downscaled TC datasets for the North Atlantic Ocean, derived from 26 CMIP6 models (`historical` and `ssp585` experiments) and `ERA5 reanalysis` data {cite:p}`Hersbach:2020` using the model developed by {cite:t}`Lin:2023`. These datasets, stored in the [Texas Advanced Computing Center (TACC)](https://tacc.utexas.edu/) [Corral storage](https://tacc.utexas.edu/systems/corral/), provide TC track and intensity information required as inputs for `pyTCR`.

To help users get started, we provide different example `Jupyter` notebooks.
These hands-on tutorials are designed for training purpose, guiding users through the key concepts and functions of `pyTCR`.

The notebooks include:

- [Downloading and preprocessing TCs data](./ex0_download_tracks_from_cmip6.ipynb)
- [Visualizing and analyzing TC tracks and densities](./ex1_tropical_cyclone_tracks.ipynb)
- [Generating TC rainfall timeseries](./ex2_rainfall_generation.ipynb)
- [Generating TC wind speed](./ex3_wind_speed_generation.ipynb)
- [Generating single rainfall event within polygons](./ex4_rainfall_polygons_generation.ipynb)
- [Generating multiple rainfall events within polygons](./ex5_multiple_rainfall_event_polygon.ipynb)
