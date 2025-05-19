<p float="left">
<img src="images/logo.gif" alt="Logo" height="45">
<img src="images/pyTCR.png" alt="pyTCR logo" width="220px">
</p>

# A tropical cyclone rainfall model for python

![](https://img.shields.io/github/license/levuvietphong/pyTCR)
![](https://img.shields.io/github/issues/levuvietphong/pyTCR)
![](https://img.shields.io/github/forks/levuvietphong/pyTCR)
![](https://img.shields.io/github/last-commit/levuvietphong/pyTCR)
![](https://img.shields.io/github/v/release/levuvietphong/pyTCR)
![](https://img.shields.io/github/actions/workflow/status/levuvietphong/pyTCR/CI.yml)

![](images/Intro-hurricane.gif)


## Table of Contents

- [Overview](#book-overview)
- [Installation](#wrench-installation)
- [Getting Started](#arrow_forward-getting-started)
- [Gallery](#framed_picture-gallery)
- [Documentation](#books-documentation)
- [License](#page_facing_up-license)
- [Contributions](#rocket-contributions)
- [Acknowledgments](#people_hugging-acknowledgments)
- [Contacts](#mailbox-contacts)

## :book: Overview
pyTCR is a physics-based model developed in python to estimate rainfall induced by tropical cyclones (TCs). It is largely based on the TCR model described by [Zhu *et al.*, 2013](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013GL058284) and [Lu *et al.*, 2018](https://journals.ametsoc.org/view/journals/atsc/75/7/jas-d-17-0264.1.xml). PyTCR simulates convective TC rainfall by correlating the precipitation rate with the total upward velocity within the TC vortex. It integrates seamlessly with outputs from [a tropical cyclone downscaling model](https://github.com/linjonathan/tropical_cyclone_risk) (see [Lin *et al.,* 2023](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023MS003686)) to generate detailed spatial-temporal rainfall patterns that align with hurricane tracks.

## :wrench: Installation

1. To install the latest stable release from `PyPI`:
    ```sh
    pip install pyTCR
    ```
    OR `conda`:
    ```sh
    conda install -c conda-forge pyTCR
    ```

    > 💡 **Note:**
    > *Python 3.11 or higher is required for installing pyTCR.*
    
2. To get the latest (unstable) development version from GitHub:
    ```sh
    git clone https://github.com/levuvietphong/pyTCR.git
    cd pyTCR
    pip install .
    ```

    > 💡 **Note:**
    > *This installs the most recent updates, which may include new features and bug fixes that are not yet available in the stable release.*

## :arrow_forward: Getting Started
This repository provides a collection of [jupyter notebooks](https://github.com/levuvietphong/pyTCR/tree/main/notebooks) designed to facilitate the use of PyTCR for generating and visualizing rainfall and wind speed induced by TCs. To access the tutorials, clone the repository locally. These notebooks cover various use cases and data sources, providing comprehensive guidance for users. 

0. [Downloading Tropical Cyclone Data:](./notebooks/ex0_download_tracks_from_cmip6.ipynb) illustrates the process of downloading CMIP6 tropical cyclone tracks that have undergone downscaling using the tropical cyclone downscaling [model](https://github.com/linjonathan/tropical_cyclone_risk).

1. [Visualizing Tropical Cyclone Tracks:](./notebooks/ex1_tropical_cyclone_tracks.ipynb) showcases the visualization of tropical cyclone tracks and their density, utilizing both observational data and downscaled CMIP6 outputs.

2. [Rainfall Simulation:](./notebooks/ex2_rainfall_generation.ipynb) uses outputs from the tropical cyclone downscaling model to simulate rainfall driven by tropical cyclones. This includes rainfall at a specific location and over the entire domain of the tropical cyclone tracks.

3. [Wind Speed Generation:](./notebooks/ex3_wind_speed_generation.ipynb) demonstrates the generation of spatially and temporally varying wind speeds from downscaled CMIP6 outputs.

4. [Rainfall within a Polygon:](./notebooks/ex4_rainfall_polygons_generation.ipynb) outlines the process of generating total rainfall within polygons loaded from shapefiles.

5. [Multiple Rainfall Events within a Polygon:](./notebooks/ex5_multiple_rainfall_event_polygon.ipynb) generates multiple total rainfall events from multiple tropical cyclone tracks within a polygon.

> [!IMPORTANT]
> Downscaled tropical cyclone tracks for various CMIP6 models (including `historical` and `ssp585` experiments) are available for download [here](https://web.corral.tacc.utexas.edu/setxuifl/tropical_cyclones/downscaled_cmip6_tracks).

## :framed_picture: Gallery
The following figures illustrate the inputs used and results obtained from PyTCR:

- **Tropical cyclone tracks:** The figure below shows the tracks of tropical cyclones downscaled from outputs from E3SM v1.0 model, Reanalysis ERA5, and best track observations.

![](images/hurricane_tracks.gif)

- **Spatial Distribution of Total Rainfall:** The image below shows the spatial distribution of total rainfall along three storm tracks that made landfall in the Southeast Texas (SETx) region.

![](images/cumulative_rain_3storms.png)

- **Time Series of Precipitation:** The figure below depicts the time series of precipitation at a specific location (29.907°N, 94.900°W) within the SETx region.

![](images/rainfall_timeseries.png)


## :books: Documentation
The official documentation is hosted on [Read the Docs](https://pytcr.readthedocs.io/en/latest/).

## :page_facing_up: License
Distributed under the MIT License. See [LICENSE](LICENSE) for more information.

## :rocket: Contributions
All contributions, bug reports, bug fixes, documentation improvements, enhancements, and ideas are welcome. A detailed overview on how to contribute can be found in the [contributing guide](CONTRIBUTING.md).

## :people_hugging: Acknowledgments
This work was supported by the [Southeast Texas Urban Integrated Field Lab](https://setx-uifl.org/) <img src="https://setx-uifl.org/wp-content/uploads/2023/08/SETx-URBAN-IFL-Logo-Full-Color-Final-300x109.png" height="22" style="vertical-align: -8px" /> project, one of four [Urban Integrated Field Laboratories](https://ess.science.energy.gov/urban-ifls/) <img src="https://ess.science.energy.gov/urban-ifls/wp-content/uploads/sites/2/2023/04/UIFL-logo-final.jpg" height="22" style="vertical-align: -8px" /> (UIFLs). The software was mostly developed at <img src="https://map.ornl.org/art/logo.png" height="22" style="vertical-align: -10px"/>. We extend our gratitude to the SETx-UIFL team for their support and collaboration.


## :mailbox: Contacts
For questions and feedback, please contact:
- Phong Le (levuvietphong@gmail.com)

<hr>

[Go to Top](#table-of-contents)
