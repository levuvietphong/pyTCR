<p float="left">
<img src="notebooks/images/logo.png" alt="Logo" width="45">
<img src="notebooks/images/pyTCR.png" alt="pyTCR logo" width="220px">
</p>

# A tropical cyclone rainfall model for python

![](https://img.shields.io/github/license/levuvietphong/pyTCR)
![](https://img.shields.io/github/issues/levuvietphong/pyTCR)
![](https://img.shields.io/github/forks/levuvietphong/pyTCR)
![](https://img.shields.io/github/last-commit/levuvietphong/pyTCR)
![](https://img.shields.io/github/downloads/levuvietphong/pyTCR/total)
![](https://img.shields.io/github/v/release/levuvietphong/pyTCR)

![](notebooks/images/Intro-hurricane.gif)

## :books: Overview
**pyTCR** is a physics-based model designed for estimating rainfall driven by tropical cylones (TC) ([Zhu *et al.*, 2013](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013GL058284); [Lu *et al.*, 2018](https://journals.ametsoc.org/view/journals/atsc/75/7/jas-d-17-0264.1.xml)). The PyTCR simulates convective TC rainfall by relating the precipitation rate to the total upward velocity within the TC vortex. The `PyTCR` seamlessly integrates with outputs from [a tropical cylone downscaling model](https://github.com/linjonathan/tropical_cyclone_risk) (see [Lin *et al.,* 2023](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023MS003686)) to produce detailed spatial-temporal rainfall patterns aligned with hurricane tracks.

## :wrench: Installation

First, clone pyTCR to your local computer:

```shell
git clone https://github.com/levuvietphong/pyTCR.git
```

Next, navigate to the pyTCR directory and create a `conda` virtual environment using the provided `environment.yml` file:
```shell
conda env create -f environment.yml
```

Activate the environment:
```shell
conda activate pyTCR
```

Finally, install the PyTCR package:
```shell
pip install -e .
```

## :arrow_forward: Getting Started
This repository provides several [`jupyter notebooks`](https://github.com/levuvietphong/pyTCR/tree/main/notebooks) to help you get started with PyTCR for generating and visualizing tropical cyclone-driven rainfall. The notebooks demonstrate different use cases and data sources:
1. **Tropical Cyclone Tracks:** This notebook plots tropical cyclone tracks and density based on observations and outputs obtained from the tropical cyclone downscaling model -- View it [here](./notebooks/ex1_tropical_cyclone_tracks.ipynb).

2. **Rainfall Generation:** This notebook utilizes outputs from [a tropical cyclone downscaling model](https://github.com/linjonathan/tropical_cyclone_risk) to generate rainfall driven by tropical cyclones. The tropical cyclone downscaling model was modified to work with CMIP6 datasets -- Check it out [here](./notebooks/ex2_rainfall_generation.ipynb).


## :framed_picture: Gallery
The following figures illustrate the inputs used in pyTCR and results obtained from PyTCR:

- **Tropical cyclone tracks:** The figure below shows the tracks of tropical cyclones downscaled from outputs from E3SM v1.0 model, Reanalysis ERA5, and best track observations.

![](notebooks/images/hurricane_tracks.gif)

- **Spatial Distribution of Total Rainfall:** The image below shows the spatial distribution of total rainfall along three storm tracks that made landfall in the Southeast Texas (SETx) region.

![](notebooks/images/cumulative_rain_3storms.png)

- **Time Series of Precipitation:** The figure below depicts the time series of precipitation at a specific location (29°N, 95°W) within the SETx region.

![](notebooks/images/rainfall_timeseries.png)


## :page_facing_up: License
Distributed under the MIT License. See [LICENSE](LICENSE) for more information.


## :people_hugging: Acknowledgements
This work was supported by the [Southeast Texas Urban Integrated Field Lab](https://setx-uifl.org/) (<img src="https://setx-uifl.org/wp-content/uploads/2023/08/SETx-URBAN-IFL-Logo-Full-Color-Final-300x109.png" height="22" style="vertical-align: -5px" />) project, one of four [Urban Integrated Field Laboratories](https://ess.science.energy.gov/urban-ifls/) (UIFLs <img src="https://ess.science.energy.gov/urban-ifls/wp-content/uploads/sites/2/2023/04/UIFL-logo-final.jpg" height="22" style="vertical-align: -4px" />). The software was mostly developed at the <img src="https://map.ornl.org/art/logo.png" height="22" style="vertical-align: -6px"/> (ORNL). We extend our gratitude to the SETx-UIFL team for their support and collaboration.


## :mailbox: Contacts
Collaborators and contributions are very welcome! For questions and feedback, please contact:
- Phong Le (lepv@ornl.gov)