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
date: 16 February 2025
bibliography: paper.bib
---

# Summary
`pyTCR` is a climatology software package developed in the Python programming language.
It integrates the capabilities of several legacy physical models and increases computational efficiency to allow rapid estimation of tropical cyclone (TC) rainfall consistent with the large-scale environment.
Specifically, `pyTCR` implements a horizontally distributed and vertically integrated model [@Zhu:2013] for simulating rainfall driven by TCs. Along storm tracks, rainfall is estimated by computing the cross-boundary-layer, upward water vapor transport caused by different mechanisms including frictional convergence, vortex stretching, large-scale baroclinic effect (i.e., wind shear), topographic forcing, and radiative cooling [@Lu:2018].
The package provides essential functionalities for modeling and interpreting spatio-temporal TC rainfall data. `pyTCR` requires a limited number of model input parameters, making it a convenient and useful tool for analyzing rainfall mechanisms driven by TCs.

To sample rare (most intense) rainfall events that are often of great societal interest, `pyTCR` adapts and leverages outputs from a statistical-dynamical TC downscaling model [@Lin:2023] capable of rapidly generating a large number of synthetic TCs given a certain climate. As a result, `pyTCR` significantly reduces computational effort and improves the efficiency in capturing extreme TC rainfall events at the tail of the distributions from limited datasets. Furthermore, the TC downscaling model is forced entirely by large-scale environmental conditions from reanalysis data or coupled General Circulation Models (GCMs), simplifying the projection of TC-induced rainfall and wind speed under future climate using `pyTCR`. Finally, `pyTCR` can be coupled with hydrological and wind models to assess risks associated with independent and compound events (e.g., storm surges and freshwater flooding).


# Statement of need
Tropical cyclones (TCs) -- that is, hurricanes and tropical storms -- are among the most destructive weather events, causing massive economic and human losses worldwide each year [@Mendelsohn:2012;@Krichene:2023]. In the United States, hurricanes can trigger a surge of deaths long after the storms through complex chains of lasting impacts [@Young:2024]. Much of the damage caused by TCs is done by water – particularly by torrential rainfall (which is sensitive to TC structural characteristics, intensity, and movement) and subsequent flooding [@Shi:2024;@Zhang:2018]. Accurately capturing TC rainfall characteristics at high spatial (e.g., <10 km) and temporal (e.g, hourly) resolution is therefore of critical importance. Moreover, a growing body of evidence suggests that TC rainfall is becoming more intense under a warming climate due to the Clapeyron–Clausius scaling of water vapor in the atmosphere [@Held:2006], increasing the likelihood of extreme rainfall and flooding [@Emanuel:2021;@Zhu:2021]. Given the societal consequences of TCs, it is crucial to understand not only TC rainfall risk in the current climate, but also how the risk might evolve with warming. Advancing tools and models to accurately and efficiently quantify these risks is of great significance.

The ability of GCMs to simulate climate extremes has been substantially improved over the past few decades, primarily those participating in the Coupled Model Intercomparison Project phase 6 or CMIP6 [@Eyring:2016;@Kim:2020]. These climate models have become one of the main tools used for exploring the effect of global warming on precipitation and climate variability [@Emanuel:2021;@Le:2021;@Le:2023]. While high-resolution GCMs (e.g., those in the HighResMIP experiments) have improved the representation of TCs [@Haarsma:2016;@Li:2018;@Zhang:2024], they remain too computationally expensive for risk analysis, which requires robust sampling of extreme rainfall events. `pyTCR` responds to this need for an easy to use and efficient tool that facilitates TC-driven rainfall analysis across scales. It takes the advantage of a synthetic downscaling approach that can generate large ensembles of synthetic TCs at the basin (ocean) scale based on comprehensive climate conditions from observations and reanalysis data and a large number of GCM simulations [@Emanuel:2008;@Lin:2023].
`pyTCR` provides a fast and highly efficient tool for risk analysis related to TC rainfall.


# Mathematical Approach
`PyTCR` implements a TC rainfall model described in @Zhu:2013 and @Lu:2018 that simulates along-track convective rainfall by relating the precipitation rate to the total upward velocity within the TC vortex. We refer to the study by @Lu:2018 for detailed formulation of this model. Here, for convenience, we give a brief overview of the main rainfall mechanisms used in this model and implemented in `pyTCR`.

Let $P_{TC}$ be the precipitation rate driven by TCs, calculated as:

$$P_{TC} = \epsilon_p \frac{\rho_{air}}{\rho_{liquid}} q_s \max(w,0)$$

where $\epsilon_p$ is precipitation efficiency, $\rho_{air}$ and $\rho_{liquid}$ are the density of water vapor and liquid water, respectively (the ratio $\rho_{air}/\rho_{liquid}\approx 0.0012$), $q_s$ is the saturation specific humidity, and $w$ is the upward-positive vertical wind velocity that brings surface moisture into the upper atmosphere.
The key assumption used here is that time-evolving TC rainfall is organized around the storm track and proportional to $w$.
The core function of `pyTCR` includes estimating $w$ as a linear combination of five major components:

\begin{equation}
w = w_f + w_h + w_t + w_s + w_r
\end{equation}

First, $w_f$ represents the velocity induced by surface frictional convergence that depends on the boundary layer formulation. This process is critical for maintaining deep convection and sustaining the TC's core rainfall. It is the dominant factor of $w$, estimated as:

\begin{equation}
w_f = \frac{-1}{r} \frac{\partial}{\partial r} \left(r^2 \frac{\tau_{\theta s}}{\partial M/\partial r}\right)
\end{equation}

in which $r$ is the radius from the storm center, $\tau_{\theta s}$ is the azimuthal surface stress, $M=rV+\frac{1}{2} f r^{2}$ is the absolute angular momentum per unit mass, $V$ is the azimuthal wind speed at radius $r$, and $f$ is the Coriolis parameter. Second, $w_h$ is the surface vertical velocity driven by topographic forcing, representing the upward motion of air as it ascends over elevated terrain:

\begin{equation}
w_h = \mathbf{V} \cdot \nabla h
\end{equation}

where $h$ is the topographic height and $\mathbf{V}$ is the total horizontal wind velocity given as the vector sum of the gradient wind $V$ (azimuthal wind at the gradient level $\sim900$ $hPa$) and environmental background wind. Third, $w_t$ denotes the vertical wind velocity component arising from time dependence of the storm's vorticity (vortex stretching):

\begin{equation}
w_t = \int_b^H \frac{1}{r} \frac{\partial}{\partial r} \left(r \frac{\partial M / \partial t}{\partial M / \partial r} \right) dz \simeq H_b \frac{1}{r} \frac{\partial}{\partial r} \left(r \frac{\partial M / \partial t}{\partial M / \partial r} \right)
\end{equation}

in which $b$ is the height of the boundary layer, $H$ is the mid-level of the troposphere, and $H_b = H-b$ is a representative depth scale of the lower troposphere.
Fourth, $w_s$ denotes the baroclinic/shear component velocity approximated as:

\begin{equation}
w_s \simeq \frac{g}{c_p(T_s-T_t)(1-\epsilon_p)N^2} V \left(f+\frac{V}{r}+\frac{\partial V}{\partial r}\right)(\Delta \mathbf{V}_e \cdot \mathbf{j})
\end{equation}

where $c_p$ is the heat capacity of dry air, $g$ is the acceleration of gravity, $N$ is the buoyancy frequency for dry air, $T_s$ is the surface temperature, $T_t$ is the tropopause temperature, $\mathbf{j}$ is the unit vector pointing radially outward from the storm center, and $\Delta\mathbf{V}_e$ is the vector wind shear across the troposphere. Finally, $w_r$ represents wind velocity related to radiative cooling, a process in which heat is lost from the storm system due to infrared radiation emitted by clouds, the storm core, and surrounding environment.
For the sake of simplicity, it is set as a constant parameter (-0.005 $m/s$) in `pyTCR` .

The inputs to the TC rainfall model then are the gradient wind $V$, total horizontal wind velocity $\mathbf{V}$, vector wind shear $\Delta \mathbf{V}_e$, topographic height $h$, and saturation specific humidity $q_s$ along the storm track. The TC downscaling model [@Lin:2023] provides `pyTCR` with 3-hourly information on TCs such as track, intensity, and size to calculate these inputs. Specifically, the gradient wind $V$ is estimated in `pyTCR` using analytical wind profile models [@Chavas:2015] based on storm characteristics. $\mathbf{V}$ is approximated in `pyTCR` as the sum of the gradient wind $V$ and storm translation. $\Delta \mathbf{V}_e$ is estimated from the geostrophic wind at 200 and 850 $hPa$. The storm-centered specific humidity $q_s$ is calculated from the 600 $hPa$ atmospheric temperature and TC intensity at each time step following @Emanuel:2017. The outputs of `pyTCR` include high-resolution, time-evolving TC rainfall and wind fields, enabling detailed analysis of storm impacts.


# Examples

`pyTCR` generates and provides access to a large set of synthetic rainfall events based on TC tracks.
To demonstrate its use, the repository includes downscaled TC datasets for the North Atlantic Ocean (NAO), derived from 26 CMIP6 models (historical and ssp585 experiments) and ERA5 reanalysis data [@Hersbach:2020] using the model developed by @Lin:2023. These datasets, stored in the Texas Advanced Computing Center (TACC) Corral storage [@CorralTACCHPC], provide TC track and intensity information required as inputs for `pyTCR`. For downscaling TCs in other ocean basins, users are referred to @Lin:2023 for further details.
We note that `pyTCR` is not restricted to this datasets and can seamlessly integrate outputs from any other TC downscaling model.

To help users get started, we provide six example `Jupyter` notebooks.
These hands-on tutorials are designed for training purpose, guiding users through the key concepts and functions of `pyTCR`. The notebooks include:

- Downloading and preprocessing TCs data
- Visualizing and analyzing TC tracks and densities
- Generating TC rainfall timeseries
- Generating TC wind speed
- Generating single rainfall event within polygons
- Generating multiple rainfall events within polygons

We briefly present two notebooks in this paper.

Figure 1 compares the tracks and mean power dissipation index (PDI) of TCs downscaled from the low-resolution outputs of the Exascale Energy Earth Model ver 1.0 [@E3SMv1.0] and ERA5 reanalysis data [@Hersbach:2020] using the TC downscaling model [@Lin:2023] with those obtained from the International Best Track Archive for Climate Stewardship (IBTrACS) observations in the NAO during the historical period (1964-2014).
The PDI quantifies the total power dissipated annually by all TCs in the ocean basin and is defined as $PDI=\int_0^{\tau}V^3_{max}dt$ [@Emanuel:2005]. Here, $V_{max}$ is the maximum sustained wind speed at the conventional measurement altitude of \mbox{10 $m$}, and $\tau$ is the lifetime of the storm. The PDI captures not only TC frequency, but also duration and intensity.
Overall, the results suggest that TCs downscaled from E3SM is able to reproduce key aspects of TC behavior in the NAO over the historical period. However, the PDI maps indicate that E3SM-downscaled TCs tend to underestimate storm activity near the tropical eastern Atlantic region.

<!-- It highlights the ability of the TC downscaling model to reasonably reproduce the general behavior of TC observed over the past period, providing confidence for analyzing TC patterns in the future climate. -->

![(Top) Tracks of 200 example TCs in the North Atlantic. Color lines indicates wind speed and TC tracks that landfall in Texas.(Bottom) Mean power dissipation index (PDI) per $2^{\circ} \times 2^{\circ}$ box per year. Plot was generated using the notebook `ex1_tropical_cyclone_tracks.ipynb` in the repository.\label{fig1}](Fig1.png)

Along each TC track from the downscaling model, `pyTCR` can generate time series and spatial patterns of rainfall events. Figure 2 illustrates the spatial distribution of total rainfall along a TC track in NAO. The TC originates in the central Atlantic (25$^\circ$N, 60$^\circ$W) and generally moves westward before making landfall in the United States. Rainfall intensity increases significantly upon landfall in Texas compared to its intensity over the ocean. Time series of rainfall at any domain influenced by the TCs can be easily extracted in `pyTCR` for other analyses.

![Illustration of the spatial distribution of total rainfall generated for a particular TC track that makes landfall on Texas, USA. Plot was generated using the notebook `ex2_rainfall_generation.ipynb` in the repository\label{fig2}](Fig2.png){width=90%}

# Acknowledgements
This work was supported by the U.S. Department of Energy, Office of Science, Biological and Environmental Research program and is a product of the Southeast Texas (SETx) Urban Integrated Field Laboratory (UIFL) project.

# References
