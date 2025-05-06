"""
This module contains physical constans for use across the project.
All constans are in SI units unless otherwise specified.
"""
import math

#Physical constants
OMEGA = math.pi / (6 * 3600)      # Earth angular velocity parameter
RHOA_OVER_RHOL = 0.00117          # ratio of air density over water density
RAD_TO_DEG = math.pi / 180  # Conversion factor from radians to degrees (π/180)
KNOTS_TO_MPS = 1852.0 / 3600.0    # Conversion factor from knots to m/s (1 knots = 0.5144 m/s)

#Empirical constants
eprecip = 0.9     #precipitation efficiency $\epsilon_p$
wrad    = -0.005  #Background subsidence velocity under radiative cooling (m/s) $w_r$
magfac  = 1.0     #Overall scale factor for storm size



