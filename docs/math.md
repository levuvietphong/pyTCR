(math)=

# Background

`PyTCR` implements a tropical cyclone (TC) rainfall model described in {cite}`Zhu:2013` and {cite}`Lu:2018` that simulates along-track convective rainfall by relating the precipitation rate to the total upward velocity within the TC vortex. We refer to the study by @Lu:2018 for detailed formulation of this model. Here, for convenience, we give a brief overview of the main rainfall mechanisms used in this model and implemented in `pyTCR`.

Let $P_{TC}$ be the precipitation rate driven by TCs, calculated as:

$$P_{TC} = \epsilon_p \frac{\rho_{air}}{\rho_{liquid}} q_s \max(w,0)$$

where $\epsilon_p$ is precipitation efficiency, $\rho_{air}$ and $\rho_{liquid}$ are the density of water vapor and liquid water, respectively (the ratio $\rho_{air}/\rho_{liquid}\approx 0.0012$), $q_s$ is the saturation specific humidity, and $w$ is the upward-positive vertical wind velocity that brings surface moisture into the upper atmosphere.
The key assumption used here is that time-evolving TC rainfall is organized around the storm track and proportional to $w$.
The core function of `pyTCR` includes estimating $w$ as a linear combination of five major components:

$$w = w_f + w_h + w_t + w_s + w_r$$

First, $w_f$ represents the velocity induced by surface frictional convergence that depends on the boundary layer formulation. This process is critical for maintaining deep convection and sustaining the TC's core rainfall. It is the dominant factor of $w$, estimated as:

$$w_f = \frac{-1}{r} \frac{\partial}{\partial r} \left(r^2 \frac{\tau_{\theta s}}{\partial M/\partial r}\right)$$

in which $r$ is the radius from the storm center, $\tau_{\theta s}$ is the azimuthal surface stress, $M=rV+\frac{1}{2} f r^{2}$ is the absolute angular momentum per unit mass, $V$ is the azimuthal wind speed at radius $r$, and $f$ is the Coriolis parameter. 

Second, $w_h$ is the surface vertical velocity driven by topographic forcing, representing the upward motion of air as it ascends over elevated terrain:

$$w_h = \mathbf{V} \cdot \nabla h$$

where $h$ is the topographic height and $\mathbf{V}$ is the total horizontal wind velocity given as the vector sum of the gradient wind $V$ (azimuthal wind at the gradient level $\sim900$ $hPa$) and environmental background wind. 

Third, $w_t$ denotes the vertical wind velocity component arising from time dependence of the storm's vorticity (vortex stretching):

$$w_t = \int_b^H \frac{1}{r} \frac{\partial}{\partial r} \left(r \frac{\partial M / \partial t}{\partial M / \partial r} \right) dz \simeq H_b \frac{1}{r} \frac{\partial}{\partial r} \left(r \frac{\partial M / \partial t}{\partial M / \partial r} \right)$$

in which $b$ is the height of the boundary layer, $H$ is the mid-level of the troposphere, and $H_b = H-b$ is a representative depth scale of the lower troposphere.

Fourth, $w_s$ denotes the baroclinic/shear component velocity approximated as:

$$w_s \simeq \frac{g}{c_p(T_s-T_t)(1-\epsilon_p)N^2} V \left(f+\frac{V}{r}+\frac{\partial V}{\partial r}\right)(\Delta \mathbf{V}_e \cdot \mathbf{j})$$

where $c_p$ is the heat capacity of dry air, $g$ is the acceleration of gravity, $N$ is the buoyancy frequency for dry air, $T_s$ is the surface temperature, $T_t$ is the tropopause temperature, $\mathbf{j}$ is the unit vector pointing radially outward from the storm center, and $\Delta\mathbf{V}_e$ is the vector wind shear across the troposphere.

Finally, $w_r$ represents wind velocity related to radiative cooling, a process in which heat is lost from the storm system due to infrared radiation emitted by clouds, the storm core, and surrounding environment. For the sake of simplicity, it is set as a constant parameter (default: -0.005 $m/s$) in `pyTCR` .
