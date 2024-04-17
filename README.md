# A python-based Tropical Cyclone Rainfall model

PyTCR is an open-source, physics-based model written in python for estimating tropical cylone rainfall (TCR; Zhu *et al.*, 2013; Lu *et al.*, 2018). The PyTCR simulates convective TC rainfall by relating the precipitation rate to the total upward velocity within the TC vortex. The PyTCR takes inputs from [a downscaled TC model](https://github.com/linjonathan/tropical_cyclone_risk) and generate spatial-temporal rainfall associated with hurricane tracks.

## Prerequisite

The following packages are required for running pyTCR:

* [pandas](https://pandas.pydata.org/)
* [geopandas](https://geopandas.org/en/stable/)
* [netcdf4](https://unidata.github.io/netcdf4-python/)
* [xarray](https://docs.xarray.dev/en/stable/)

## Getting Started
See notebooks for examples.

## Acknowledgements
This work is supported by the Southeast Texas Urban Integrated Field Lab (SETx-UIFL) project.

## References
1. Zhu, L., Quiring, S. M., and Emanuel, K. A. (2013), Estimating tropical cyclone precipitation risk in Texas, *Geophys. Res. Lett.*, 40, 6225–6230, doi:10.1002/2013GL058284.

2. Lu, P., N. Lin, K. Emanuel, D. Chavas, and J. Smith, 2018: Assessing Hurricane Rainfall Mechanisms Using a Physics-Based Model: Hurricanes Isabel (2003) and Irene (2011). *J. Atmos. Sci.*, 75, 2337–2358, doi:10.1175/JAS-D-17-0264.1.

