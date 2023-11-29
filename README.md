# A python-based Tropical Cyclone Rainfall model

An open-source, physics-based model for estimating tropical cylone rainfall (TCR; Zhu et al., 2013). TCR simulates convective TC rainfall by relating the precipitation rate to the total upward velocity within the TC vortex. The model takes inputs from TC tracks and generate spatial-temporal rainfall associated with tracks.

## Prerequisite

The following packages are required for running pyTCR:

* [pandas](https://pandas.pydata.org/)
* [geopandas](https://geopandas.org/en/stable/)
* [netcdf4](https://unidata.github.io/netcdf4-python/)
* [xarray](https://docs.xarray.dev/en/stable/)

## References
Zhu, L., Quiring, S. M., and Emanuel, K. A. (2013), Estimating tropical cyclone precipitation risk in Texas, Geophys. Res. Lett., 40, 6225–6230, doi:10.1002/2013GL058284.