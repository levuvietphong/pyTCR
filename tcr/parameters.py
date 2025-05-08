"""
Physical constants and parameters used across functions in the TCR model.
"""
import math

# PHYSICAL CONSTANTS
OMEGA = math.pi / (12 * 3600)   # Earth angular velocity parameter
RHOA_OVER_RHOL = 0.00117        # ratio of air density over water density
RAD2DEG = math.pi / 180         # Conversion factor from radians to degrees (π/180)
DEG2KM = 60 * 1.852             # Conversion factor from degree to kilometers
                                # (1 degree = 60 nm, 1 nm = 1.852 km)
KNOTS2MPS = 1852.0 / 3600.0     # Conversion factor from knots to m/s (1 knots = 0.5144 m/s)

# RAINFALL PARAMETERS
eprecip = 0.9       # precipitation efficiency $\epsilon_p$ (dimensionless)
wrad = -0.005       # Background subsidence velocity under radiative cooling $w_r$ (m/s)
q900 = 0.01         # Default value of specific humidity at 900hPa if T600 is not available (gm/gm)
Htrop = 4000        # Depth of lower troposphere (m)
timeresw = 2        # Native time resolution of WRT output (hour)
deltar = 2          # Delta radius for calculating dM/dr (km)
wheight = 30        # Altitude above local terrain to estimate surface winds (m)

# TIME SERIES PARAMETERS
timeres = 0.5       # Time resolution, should divide evenly into 2 (hours)
timelength = 96     # Length of time series at fixed points (even number)

# STORM GEOMETRY PARAMETERS
wprofile = 3        # Wind profile 1=Holland, 2=Emanuel, 3=Emanuel & Rotunno 2011
magfac = 1.0        # Overall scale factor for storm size (dimensionless)

# PARAMETERS FOR STORM SNAPSHOTS AND SWATH MAPS
deltax = 5          # Longitudinal distance of map boundaries from storm center (degrees)
deltay = 4          # Latitudinal distance of map boundaries from storm center (degrees)
dellatlong = 0.05   # Horizontal resolution of field maps (degrees)
dellatlongq = 0.15  # Horizontal resolution of wind quiver map (degrees)
dellatlongs = 0.15  # (Lower) Horizontal resolution of swath maps (degrees)
radcity = 300       # Distance of storm from point of interest beyond which (km)
