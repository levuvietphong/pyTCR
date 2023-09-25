import os
import numpy as np

"""
Parameters for the TCR model
"""

# Parameters used by rain algorithms
wrad = 0.005       #    0.005  |  m/s  | Background subsidence velocity under radiative cooling
topores = 'high'   #   'high'  |   -   | Topgraphic resolution ('high'=0.1 degree, 'low'=0.25 degree)   
q900 = 0.01        #     0.01  | gm/gm | Default value of specific humidity at 900 hPa used if T600 is not available
eprecip = 0.9      #      0.9  |   -   | Precipitation efficiency
Cdrag = 0.0015     #    0.0015 |   -   | Surface drag coefficient over water
Cdland = 0.003     #    0.003  |   -   | Surface drag coefficient over land
Hcrit = 300        #      300  |   m   | Altitude governing transition of drag coefficient from sea to land
Htrop = 4000       #     4000  |   m   | Depth of lower troposphere
timeresw = 2       #        2  |  hrs  | Native time resolution of WRT output
deltar = 2         #        2  |   km  | Delta radius for calculating dM/dr

# Storm geometry
randfac = 'n'
seyewall = 'n'     #    'y'    |   -   | Use secondary eyewall information?
wprofile = 3       #     3     |   -   | Wind profile 1=Holland, 2=Emanuel (from Lilly volume), 3= Emanuel & Rotunno, 2011
magfac = 1.0       #    1.0    |   -   | Overall scale factor for storm size

# Parameters used in storm snapshots (windfield.m, windfield_quiver.m, and rainfield.m) and in swath maps   
plottrack = 'y'    #     'y'   |   -   | Plot storm track?
dellatlong = 0.05  #    0.05   |degrees| Horizontal resolution of field maps
dellatlongq = 0.15 #    0.15   |degrees| Horizontal resolution of wind quiver map
dellatlongs = 0.15 #    0.15   |degrees| (Lower) Horizontal resolution of swath maps (should not be less than about 0.15 for most machines)
radcity = 300      #    300    |   km  | Distance of storm from point of interest beyond which influence of storm is ignored   
bound = 'auto'     #   'auto'  |   -   | Automatic (auto) or manual (manu) specification of map bounds
deltax = 5         #      5    |degrees| Longitudinal distance of map boundaries from storm center if bound  =  'auto'
deltay = 4         #      4    |degrees| Latitudinal distance of map boundaries from storm center if bound  =  'auto'
bxmin = 20         #      -    |degrees| Minimum longitude of map if bounds = 'manu'
bxmax = 380        #      -    |degrees| Maximum longitude of map if bounds = 'manu'
bymin = -60        #      -    |degrees| Minimum latitude of map if bounds = 'manu' 
bymax = 60         #      -    |degrees| Maximum latitude of map if bounds='manu'

# Parameters for mapping
projection = 'Lambert'   #     |   -   | Map projection. See matlab mapping toolbox or m_map/map.html for documentation
gproject = 'Robinson'    #     |   -   | Map projection used for global plots
mapmode = 'auto'   #   'auto'  |   -   | Mode of determining lat-long bounds of map ('auto' or 'manu')
dellat = 2         #     2     |degrees| Width of latitude buffers relative to storm track limits (affects 'auto' mode only)
dellong = 2        #     2     |degrees| Width of longitude buffers relative to storm track limits (affects 'auto' mode only)
longmin = 0        #     -     |degrees| Lower bound on longitude plotted ('manu' mode only, 0 to 360)    
longmax = 359      #     -     |degrees| Upper bound on longitude plotted ('manu' mode only, 0 to 360)    
latmin = -50       #     -     |degrees| Lower bound on latitude plotted ('manu' mode only, -90 to 90)    
latmax = 50        #     -     |degrees| Upper bound on latitude plotted ('manu' mode only, -90 to 90)    

peakv = 40         #     40    | knots | Maximum wind speed that must be exceeded for event to be included
minv = 35          #     35    | knots | Minimum wind speed to be counted in track density plots 
nmax = 1e6         #    1e6    |   -   | Maximum number of points to include
mres = 2           #     2     |degrees| Resolution of track and power dissipation density maps
mresbest = 3       #     3     |degrees| Resolution of best track density maps

# Time series parameters
timeres = 0.5      #    0.5    | hours | Time resolution for time series at fixed points (should divide evenly into 2)
timelength = 96    #     96    | hours | Length of time series at fixed points (must be an even number)
bthresh = 30       #     30    | knots | Wind speed above which best tracks counted in calculating duration
