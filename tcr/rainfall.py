"""
Functions for rainfall in PyTCR
"""

import datetime
import math
import pandas as pd
import numpy as np
from tcr import wind as tcr_wind
from tcr import terrain_boundary as tcr_tb
from tcr import iodata as tcr_io
from tcr import params


def rainfieldx(nt, latstore, longstore, rmstore, vstore, rmsestore, vsestore, ut, vt,
               u850store, v850store, monthstore, daystore, hourstore, monthplot, dayplot, hourplot):
    """
    Calculates the distribution of surface rain rate (mm/hr) for a given storm at a specified time.

    Inputs:
    -------
        - nt : Index representing the storm track number.
        - latstore : latitude values for each time step and storm.
        - longstore : longitude values for each time step and storm.
        - rmstore : the radii of maximum wind for each time step and storm.
        - vstore : maximum circular wind speed (in knots) for each time step and storm.
        - rmsestore : secondary radii of maximum wind for each time step and storm.
        - vsestore : secondary maximum circular wind speed (in knots) for each time step and storm.
        - ut : west-east component of the storm translation velocity
        - vt : north-south component of the storm translation velocity
        - u850store : u-component wind speed at 850 hPa for each time step and storm.
        - v850store : v-component wind speed at 850 hPa for each time step and storm.
        - monthstore : month values (1-12) for each time step.
        - daystore : day values for each time step.
        - hourstore : hour values (in GMT) for each time step.
        - monthplot : The month (1-12) for which the rain rate is to be calculated.
        - dayplot : The day of the month (1-31) for which the rain rate is to be calculated.
        - hourplot : The hour of the day (0-23, GMT) for which the rain rate is to be calculated.

    Returns:
    --------
        - rainrate: surface rainfall rates (in mm/hr) for each storm storm track

    """

    magfac = params.magfac
    deltax = params.deltax
    deltay = params.deltay
    bxmin = params.bxmin                # boundary in x-y direction
    bxmax = params.bxmax
    bymin = params.bymin
    bymax = params.bymax
    dellatlong = params.dellatlong
    q900 = params.q900
    eprecip = params.eprecip            # Precipitation efficiency
    wrad = params.wrad                  # Background subsidence velocity

    m_to_mm = 1000                      # convert unit from meter to milimeter 
    rhoa_over_rhol = 0.00117            # ratio of air density over water density

    nrm, mrm = np.shape(rmstore)
    rfac = magfac * (1+np.zeros((nrm, mrm)))
    pifac = math.acos(-1)/180
    knotfac = 1852./3600

    # Get the length of each event by finding the first 0 element
    jmaxd = np.argmin(vstore, axis=1)

    dum = np.zeros((jmaxd[nt]))
    V = np.column_stack((dum+2000, monthstore[nt, 0:jmaxd[nt]].T,
                        daystore[nt, 0:jmaxd[nt]].T, hourstore[nt, 0:jmaxd[nt]].T, dum, dum))
    timev = pd.to_datetime(dict(
        year=V[:, 0], month=V[:, 1], day=V[:, 2], hour=V[:, 3], minute=V[:, 4], second=V[:, 5]))
    timeplot = datetime.datetime(2000, monthplot, dayplot, hourplot, 0, 0)

    diff0 = timeplot - timeplot
    timediff = timeplot - timev
    timediff[timediff < diff0] = diff0
    jplot = np.argmin(timediff)
    jplot = np.maximum(jplot, 1)
    jplot = np.minimum(jplot, 398)

    latstorm = latstore[nt, jplot-1:jplot+2]
    longstorm = longstore[nt, jplot-1:jplot+2]
    vstorm = vstore[nt, jplot-1:jplot+2]
    rmstorm = rfac[nt, jplot-1:jplot+2] * [rmstore[nt, jplot-1:jplot+2]]
    vsestorm = vsestore[nt, jplot-1:jplot+2]
    rmsestorm = rfac[nt, jplot-1:jplot+2] * [rmsestore[nt, jplot-1:jplot+2]]
    utstorm = ut[nt, jplot-1:jplot+2]
    vtstorm = vt[nt, jplot-1:jplot+2]
    ush = np.zeros((3))
    vsh = np.zeros((3))

    if 'u850store' in locals():
        vdrift = 1.5*3600/1852
        vdrift = vdrift*latstore[1, 1]/(abs(latstore[1, 1])+1e-8)
        u850storm = u850store[nt, jplot-1:jplot+2]
        v850storm = v850store[nt, jplot-1:jplot+2]
        ush = 5 * knotfac * (utstorm-u850storm)
        vsh = 5 * knotfac * (vtstorm-vdrift * np.cos(pifac*latstorm) - v850storm)

    bxmin = np.floor(longstorm[1]-deltax)
    bxmax = np.ceil(longstorm[1]+deltax)
    bymin = np.floor(latstorm[1]-deltay)
    bymax = np.ceil(latstorm[1]+deltay)

    h, hx, hy, x, y = tcr_tb.estimate_topographic_height(
        bxmin, bxmax, bymin, bymax, dellatlong
    )

    w = tcr_wind.pointwfield(latstorm, longstorm, vstorm, rmstorm, vsestorm,
                             rmsestorm, utstorm, vtstorm, ush, vsh, y, x, h, hx, hy)

    temp = eprecip * m_to_mm * 3600 * rhoa_over_rhol * q900 * np.maximum(w[0, 1, :, :]-wrad, 0)

    rainrate = temp
    rainrate = rainrate.transpose()

    return rainrate, x, y


def rainswathx(nt, latstore, longstore, rmstore, vstore, rmsestore, vsestore,
               ut, vt, u850store, v850store):
    """
    This script calculates the distribution of accumulated precipitation for a given
    individual storm.

    Inputs:
        - nt: Track number of the storm
        - latstore, longstore: Latitudes and longitudes along each track
        - vstore: Maximum circular wind along each storm track
        - rmstore: Radius (in km) of maximum circular wind along each track
        - vsestore: maximum circular wind of any secondary eyewalls that may be present
        - rmsestore: Radius (in km) of maximum circular wind of any secondary eyewalls
        - ut: West-east component of the storm translation velocity
        - u850store, v850store: Zonal & meridional components of the 850 hPa environmental
                wind speed (knots)
    Returns:
        - x, y: vectors containing the longitudes and latitudes of the grid
        - netrain: storm total rainfall (unit: mm) at each point on the grid
    """
    magfac = params.magfac              # overall scale factor for storm size
    deltax = params.deltax              # longitudinal distance of map boundaries from storm center
    deltay = params.deltay              # latitudinal distance of map boundaries from storm center
    bxmin = params.bxmin                # minimum longitude of map (degree)
    bxmax = params.bxmax                # maximum longitude of map (degree)
    bymin = params.bymin                # minimum latitude of map (degree)
    bymax = params.bymax                # maximum latitude of map (degree)
    dellatlongs = params.dellatlongs    # horizontal resolution of field maps
    q900 = params.q900                  # specific humidity at 900 hPa
    timeres = params.timeres            # time resolution for time series at fixed points
    wrad = params.wrad                  # background subsidence velocity under radiative cooling
    eprecip = params.eprecip            # precipitation efficiency

    m_to_mm = 1000
    rhoa_over_rhol = 0.00117            # $rho_{air} / rho_{liquid}$

    nrm, mrm = np.shape(rmstore)
    rfac = magfac * (1+np.zeros((nrm, mrm)))

    if bxmin < 0:
        bxmin += 360
    if bxmax < 0:
        bxmax += 360

    # Initialize variables
    latdata = latstore[nt, :]
    latdata = latdata[(latdata != 0) & ~np.isnan(latdata)]
    latsize = len(latdata)
    utd = ut[nt, :latsize].reshape((1, latsize))
    vtd = vt[nt, :latsize].reshape((1, latsize))
    ush = np.zeros_like(utd)
    vsh = np.zeros_like(vtd)
    vdrift = 1.5 * 3600 / 1852
    vdrift *= latstore[0, 0] / (np.abs(latstore[0, 0]) + 1e-8)

    if 'u850store' in locals():
        ush = 5 * 1852 / 3600 * (utd - u850store[nt, :latsize])
        vsh = (
            5 * 1852 / 3600
            * (
                vtd
                - vdrift * np.cos(np.pi / 180 * latstore[nt, :latsize])
                - v850store[nt, :latsize]
            )
        )

    lat = latstore[nt, :latsize].reshape((1, latsize))
    long = longstore[nt, :latsize].reshape((1, latsize))

    # Convert longitude to 0 to 360 degree east
    long[long < 0] += 360
    v = vstore[nt, :latsize].reshape((1, latsize))
    vse = vsestore[nt, :latsize].reshape((1, latsize))

    # Scale and randomize radii of maximum wind
    nrm, mrm = rmstore.shape
    jmaxd = latsize
    rfac = np.ones((nrm, mrm))

    temp = magfac * rfac[nt, :] * rmstore[nt, :]
    rm = temp[:latsize].reshape((1, latsize))
    nrm = np.shape(rm)[1]
    temp = magfac * rfac[nt, :] * rmsestore[nt, :]
    rmse = temp[:latsize].reshape((1, latsize))

    for i in range(jmaxd):
        if long[0, 0] > 200 and long[0, i] < 50:
            long[0, i] += 360

    bxmin = np.min(long[np.nonzero(long)]) - deltax
    bxmax = np.max(long[np.nonzero(long)]) + deltax
    bymin = np.min(lat[np.nonzero(lat)]) - deltay
    bymax = np.max(lat[np.nonzero(lat)]) + deltay

    # Get topographic height and its gradients
    h, hx, hy, x, y = tcr_tb. estimate_topographic_height(bxmin, bxmax, bymin, bymax, dellatlongs)

    # Calculate vertical velocity time series
    w = tcr_wind.pointwshortn(lat, long, v, rm, vse, rmse, utd, vtd,
                              ush, vsh, x, y, h, hx, hy, timeres)
    wq = np.maximum(w-wrad, 0) * q900
    netrain = (
        eprecip * m_to_mm * timeres * 3600 * rhoa_over_rhol * np.sum(wq, axis=(0, 1))
    )

    return x, y, netrain.T


def raingen(plat, plong, latstore, longstore, datestore, vstore, rmstore, vsestore, rmsestore,
            u850store, v850store, utrans, vtrans, T600=None):
    """
    Calculate accumulated rain and rainrates at specified locations for the active event set

    Inputs:
        - plat, plong: Latitude and longitude of points of interest.
        - latstore, longstore: Latitudes and longitude along each storm track
        - datestore: datetime info in integer format
        - vstore: Maximum circular wind along each storm track
        - rmstore: Radius (in km) of maximum circular wind along each track
        - vsestore: maximum circular wind of any secondary eyewalls that may be present
        - rmsestore: Radius (in km) of maximum circular wind of any secondary eyewalls
        - u850store: Zonal component of the 850 hPa environmental wind speed (in knots)
        - utrans: West-east component of the storm translation velocity
        - vtrans: North-south component of the storm translation velocity

    Returns:
        - rain: Total storm rainfall (unit: mm)
        - rainrate: Rain rate (unit: mm/hr)
        - dayk: Time in date format corresponding to rainrate
    """

    ut = np.nan_to_num(utrans)
    vt = np.nan_to_num(vtrans)

    # Convert point lat,lon to array
    plat = np.array([plat])
    plong = np.array([plong])

    magfac = params.magfac                  # overall scale factor for storm size
    q900_default = params.q900              # specific humidity at 900 hPa
    timeres = params.timeres                # time resolution for time series at fixed points
    wrad = params.wrad                      # background subsidence velocity under radiative cooling
    eprecip = params.eprecip                # precipitation efficiency

    sx = plong.size                         # get number of points
    # sy = plat.size

    # Load high-resolution bathymetry from matlab file
    mat = tcr_io.load_Matlab_data('data', 'bathymetry_high.mat')
    bathy = mat["bathy"]
    ntopo, _ = np.shape(bathy)
    topores = 360. / ntopo                  # topo resolution in degree
    toporesi = 1. / topores                 # inverse of topo resolution
    sfac = 1. / (topores * 60. * 1852)      # factor converting degree to m
    pifac = math.acos(-1) / 180             # pi number
    knotfac = 1852. / 3600                  # convert nautical mile to m/s
    m, n = ut.shape                         # m: num of storm; n: num of time steps

    if np.min(plong) < 0:
        plong += 360

    ush = np.zeros((n, m))
    vsh = np.zeros((n, m))
    vdrift = 1.5 * 3600 / 1852
    vdrift = vdrift * latstore[0, 0] / (abs(latstore[0, 0]) + 1e-8)

    if "u850store" in locals():
        ush = 5 * knotfac * (ut - u850store)
        vsh = 5 * knotfac * (vt - vdrift * np.cos(pifac * latstore) - v850store)

    lat = latstore.copy()
    long = longstore.copy()
    dates = datestore.copy()
    v = vstore.copy()
    vse = vsestore.copy()

    nrm, mrm = rmstore.shape
    rfac = magfac * np.ones((nrm, mrm))

    rm = rmstore * rfac
    rmse = rmsestore * rfac

    m_to_mm = 1000                          # meter to millimeter
    rhoa_over_rhol = 0.00117                # density of air over density of liquid
    bathy = np.maximum(bathy, -1)

    # Calculate topographic and its gradients
    h, hx, hy = tcr_tb.calculate_spatial_derivatives(
        bathy, plong, plat, sx, 1, sfac, pifac, ntopo, toporesi)

    # Calculates saturation specific humidity
    if T600 is not None:
        q900, _ = tcr_tb.calculate_qs900(T600, vstore)
    else:
        q900 = np.where(ut != 0, q900_default, 0)

    # Estimate vertical wind velocity
    wq, dayk = tcr_wind.pointwshortnqdx(
        lat, long, dates, q900, v, rm, vse, rmse, ut, vt, ush, vsh,
        plong, plat, h, hx, hy, timeres, wrad)

    # convert dayk to pandas format
    datetimes = pd.to_datetime(dayk.flatten(), unit='s')
    datetimes_numpy = np.array(
        [
            (
                pd.Timestamp(dt).strftime("%Y-%m-%d %H:%M:%S")
                if not pd.isnull(dt)
                else np.nan
            )
            for dt in datetimes
        ]
    )
    dayk = np.reshape(datetimes_numpy, wq.shape)

    rainrate = eprecip * m_to_mm * 3600 * rhoa_over_rhol * wq
    rainrate = np.nan_to_num(rainrate)
    rain = timeres * np.sum(rainrate, axis=1).reshape(-1, 1)
    return rain, rainrate, dayk
