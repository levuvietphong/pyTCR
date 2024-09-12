"""
Parameters for the TCR (Tropical Cyclone Risk) model

This module defines a dictionary of parameters used throughout the TCR model,
including settings for rain algorithms, storm geometry, mapping, and time series analysis.
"""

params = {
    # Rain algorithm parameters
    "wrad": 0.005,     # Background subsidence velocity under radiative cooling (m/s)
    "topores": 'high', # Topographic resolution ('high' = 0.1 deg, 'low' = 0.25 deg)
    "q900": 0.01,      # Default specific humidity at 900 hPa if T600 is unavailable (g/g)
    "eprecip": 0.9,    # Precipitation efficiency
    "Cdrag": 0.0015,   # Surface drag coefficient over water
    "Cdland": 0.003,   # Surface drag coefficient over land
    "Hcrit": 300,      # Altitude governing transition of Cd from sea to land (m)
    "Htrop": 4000,     # Depth of lower troposphere (m)
    "timeresw": 2,     # Native time resolution of WRT output (hours)
    "deltar": 2,       # Delta radius for calculating dM/dr (km)

    # Storm geometry parameters
    "randfac": 'n',    # Random factor (unused)
    "seyewall": 'n',   # Use secondary eyewall information?
    "wprofile": 3,     # Wind profile (1=Holland, 2=Emanuel, 3=Emanuel & Rotunno 2011)
    "magfac": 1.0,     # Overall scale factor for storm size

    # Storm snapshot and swath map parameters
    "dellatlong": 0.05,  # Horizontal resolution of field maps (degrees)
    "dellatlongq": 0.15, # Horizontal resolution of wind quiver map (degrees)
    "dellatlongs": 0.15, # (Lower) Horizontal resolution of swath maps (degrees)
    "radcity": 300,      # Distance from point of interest beyond which storm influence is ignored (km)
    "deltax": 5,         # Longitudinal distance of map boundaries from storm center (degrees)
    "deltay": 4,         # Latitudinal distance of map boundaries from storm center (degrees)
    "bxmin": 20,         # Minimum longitude of map (degrees)
    "bxmax": 380,        # Maximum longitude of map (degrees)
    "bymin": -60,        # Minimum latitude of map (degrees)
    "bymax": 60,         # Maximum latitude of map (degrees)

    # Wind reduction parameter
    "wheight": 30,     # Altitude above local terrain to estimate surface winds (m)

    # Mapping parameters
    "dellat": 2,       # Width of latitude buffers relative to storm track limits (degrees)
    "dellong": 2,      # Width of longitude buffers relative to storm track limits (degrees)
    "longmin": 0,      # Lower bound on longitude plotted (0 to 360 degrees)
    "longmax": 359,    # Upper bound on longitude plotted (0 to 360 degrees)
    "latmin": -50,     # Lower bound on latitude plotted (-90 to 90 degrees)
    "latmax": 50,      # Upper bound on latitude plotted (-90 to 90 degrees)

    # Event inclusion threshold
    "peakv": 40,       # Maximum wind speed that must be exceeded for event inclusion (knots)

    # Time series parameters
    "timeres": 0.5,    # Time resolution for time series at fixed points (hours)
    "timelength": 96,  # Length of time series at fixed points (hours, even number)
    "bthresh": 30,     # Wind speed threshold for best track duration calculation (knots)
}