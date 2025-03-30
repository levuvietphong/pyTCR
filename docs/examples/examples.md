(Examples)=

# Examples
`pyTCR` generates and provides access to a large set of synthetic rainfall events based on TC tracks. The repository includes downscaled TC datasets for the North Atlantic Ocean (NAO), derived from 26 CMIP6 models (historical and ssp585 experiments) and ERA5 reanalysis data {cite:p}`Hersbach:2020` using the model developed by {cite:t}`Lin:2023`. These datasets, stored in the Texas Advanced Computing Center (TACC) Corral storage, provide TC track and intensity information required as inputs for `pyTCR`.

To help users get started, we provide six example `Jupyter` notebooks.
These hands-on tutorials are designed for training purpose, guiding users through the key concepts and functions of `pyTCR`. The notebooks include:

- Downloading and preprocessing TCs data
- Visualizing and analyzing TC tracks and densities
- Generating TC rainfall timeseries
- Generating TC wind speed
- Generating single rainfall event within polygons
- Generating multiple rainfall events within polygons
