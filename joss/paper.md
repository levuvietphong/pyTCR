---
title: 'pyTCR: A tropical cyclone rainfall model for python'
tags:
  - Python
  - xarray
  - geopandas
  - cartopy
  - netcdf4
authors:
  - name: Phong V. V. Le
    orcid: 0000-0001-5558-1023
    equal-contrib: false
    corresponding: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Geeta Persad
    equal-contrib: false
    affiliation: 2
  - name: Gabriel Perez
    equal-contrib: false
    affiliation: 3
  - name: Ifeanyichukwu C. Nduka
    equal-contrib: false
    affiliation: 2
  - name: William Mobley
    equal-contrib: false
    affiliation: 4
  - name: Ethan T. Coon
    equal-contrib: false
    affiliation: 1

affiliations:
 - name: Environmental Sciences Division, Oak Ridge National Laboratory, Oak Ridge, TN, USA
   index: 1
 - name: Jackson School of Geosciences, The University of Texas at Austin, Austin, TX, USA 
   index: 2
 - name: Department of Civil and Environmental Engineering, Oklahoma State University, OK, USA
   index: 3
 - name: Texas Advanced Computing Center, The University of Texas at Austin, Austin, TX, USA 
   index: 4
date: 16 March 2025
bibliography: paper.bib
---

# Summary
`pyTCR` is a climatology software package developed in the Python programming language.
It integrates the capabilities of several legacy physical models and increases computational efficiency to allow rapid estimation of tropical cyclone (TC) rainfall consistent with the large-scale environment.
Specifically, `pyTCR` implements a horizontally distributed and vertically integrated model [@Zhu:2013] for simulating rainfall driven by TCs. Along storm tracks, rainfall is estimated by computing the cross-boundary-layer, upward water vapor transport caused by different mechanisms including frictional convergence, vortex stretching, large-scale baroclinic effect (i.e., wind shear), topographic forcing, and radiative cooling [@Lu:2018].
The package provides essential functionalities for modeling and interpreting spatio-temporal TC rainfall data. `pyTCR` requires a limited number of model input parameters, making it a convenient and useful tool for analyzing rainfall mechanisms driven by TCs.

To sample rare (most intense) rainfall events that are often of great societal interest, `pyTCR` adapts and leverages outputs from a statistical-dynamical TC downscaling model [@Lin:2023] capable of rapidly generating a large number of synthetic TCs given a certain climate. As a result, `pyTCR` significantly reduces computational effort and improves the efficiency in capturing extreme TC rainfall events at the tail of the distributions from limited datasets. Furthermore, the TC downscaling model is forced entirely by large-scale environmental conditions from reanalysis data or coupled General Circulation Models (GCMs), simplifying the projection of TC-induced rainfall and wind speed under future climate using `pyTCR`. Finally, `pyTCR` can be coupled with hydrological and wind models to assess risks associated with independent and compound events (e.g., storm surges and freshwater flooding).


# Statement of need
Tropical cyclones (TCs) -- that is, hurricanes and tropical storms -- are among the most destructive weather events, causing massive economic and human losses worldwide each year [@Krichene:2023]. In the United States, hurricanes can trigger a surge of deaths long after the storms through complex chains of lasting impacts [@Young:2024]. Much of the damage caused by TCs is done by water – particularly by torrential rainfall and subsequent flooding [@Shi:2024;@Zhang:2018]. Accurately capturing TC rainfall characteristics at high spatial (e.g., <10 km) and temporal (e.g, hourly) resolution is therefore of critical importance. Moreover, a growing body of evidence suggests that TC rainfall is becoming more intense under a warming climate due to the Clapeyron–Clausius scaling of water vapor in the atmosphere [@Held:2006], increasing the likelihood of extreme rainfall and flooding [@Zhu:2021]. Given the societal consequences of TCs, it is crucial to understand not only TC rainfall risk in the current climate, but also how the risk might evolve with warming. Advancing tools and models to accurately and efficiently quantify these risks is of great significance.

The ability of GCMs to simulate climate extremes has been substantially improved over the past few decades, primarily those participating in the Coupled Model Intercomparison Project phase 6 or CMIP6 [@Eyring:2016]. These climate models have become one of the main tools used to explore the effect of global warming on precipitation and climate variability [@Emanuel:2021;@Le:2021;@Le:2023]. While high-resolution GCMs have improved the representation of TCs [@Haarsma:2016;@Li:2018], they remain too computationally expensive for risk analysis, which requires robust sampling of extreme rainfall events. `pyTCR` responds to this need for an easy to use and efficient tool that facilitates TC-driven rainfall analysis across scales. It takes the advantage of a synthetic downscaling approach that can generate large ensembles of synthetic TCs at the basin (ocean) scale based on comprehensive climate conditions from observations and reanalysis data and a large number of GCM simulations [@Emanuel:2008;@Lin:2023].
`pyTCR` provides a fast and highly efficient tool for risk analysis related to TC rainfall.


# Mathematical Approach
`PyTCR` implements a TC rainfall model described in @Lu:2018 that simulates along-track convective rainfall by relating the precipitation rate to the total upward velocity within the TC vortex. Let $P_{TC}$ be the precipitation rate driven by TCs, calculated as:

\begin{equation}
P_{TC} = \epsilon_p \frac{\rho_{air}}{\rho_{liquid}} q_s \max(w,0)
\end{equation}

where $\epsilon_p$ is precipitation efficiency, $\rho_{air}$ and $\rho_{liquid}$ are the density of water vapor and liquid water, respectively, $q_s$ is the saturation specific humidity, and $w$ is the upward-positive vertical wind velocity that brings surface moisture into the upper atmosphere.
The key assumption here is that time-evolving TC rainfall is organized around the storm track and proportional to $w$. The core function of `pyTCR` includes estimating $w$ as a linear combination of five major components:

\begin{equation}
w = w_f + w_h + w_t + w_s + w_r
\end{equation}

where $w_f$ represents the velocity induced by surface frictional convergence, $w_h$ is the velocity driven by topographic forcing, $w_t$ denotes the velocity component arising from time dependence of the storm's vorticity, $w_s$ denotes the baroclinic/shear component velocity, and $w_r$ represents the velocity related to radiative cooling. We refer to @Lu:2018 for detailed formulations of these components.


# Examples

`pyTCR` generates and provides access to a large set of synthetic rainfall events based on TC tracks. To help users learn key concepts and functions of `pyTCR`, we provide six `Jupyter` notebooks designed for training purposes:

- Downloading and preprocessing TCs data
- Visualizing and analyzing TC tracks
- Generating TC rainfall timeseries
- Generating TC wind speed
- Generating single rainfall event within polygons
- Generating multiple rainfall events within polygons

We briefly present two notebooks here.
Figure 1 compares the tracks and mean power dissipation index (PDI) [@Emanuel:2005] of TCs downscaled from the outputs of the E3SM-1.0 model and ERA5 reanalysis data [@Hersbach:2020] using the TC downscaling model [@Lin:2023] with those obtained from IBTrACS observations in the North Atlantic Ocean (NAO) during the historical period (1964-2014). Overall, the results suggest that TCs downscaled from E3SM-1.0 is able to reproduce key aspects of TC behavior in the NAO over the historical period. However, the PDI maps indicate that E3SM-downscaled TCs tend to underestimate storm activity near the tropical eastern Atlantic region.

![(Top) Tracks of 200 example TCs in the North Atlantic. Color lines indicates wind speed and TC tracks that landfall in Texas, USA. (Bottom) Mean power dissipation index (PDI) per $2^{\circ} \times 2^{\circ}$ box per year. Plot was generated using the notebook `ex1_tropical_cyclone_tracks.ipynb` in the repository.\label{fig1}](Fig1.png)

Along each TC track from the downscaling model, `pyTCR` can generate time series and spatial patterns of rainfall events. Figure 2 illustrates the spatial distribution of total rainfall along a TC track in NAO. The TC originates in the central Atlantic (25$^\circ$N, 60$^\circ$W) and generally moves westward before making landfall in the US. Rainfall intensity increases significantly upon landfall in Texas compared to its intensity over the ocean. Time series of rainfall at any domain influenced by the TCs can be extracted in `pyTCR` for analyses.

![Illustration of the spatial distribution of total rainfall generated for a particular TC track that makes landfall on Texas, USA. Plot was generated using the notebook `ex2_rainfall_generation.ipynb` in the repository\label{fig2}](Fig2.png){width=90%}

# Acknowledgements
This work was supported by the U.S. Department of Energy, Office of Science, Biological and Environmental Research program and is a product of the Southeast Texas (SETx) Urban Integrated Field Laboratory (UIFL) project.

# References
