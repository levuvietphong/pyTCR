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

    # Get the length of each event by finding the first 0 element, if non-zero, get all length
    jmaxd = np.argmin(vstore, axis=1)
    duration = jmaxd[nt]
    if duration == 0:
        duration = vstore.shape[1]

    dum = np.zeros((duration))
    V = np.column_stack((dum+2000, monthstore[nt, 0:duration].T,
                        daystore[nt, 0:duration].T, hourstore[nt, 0:duration].T, dum, dum))
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
    Calculate the distribution of accumulated precipitation for a given individual storm.

    Parameters:
    -----------
    nt : int
        Track number of the storm
    latstore, longstore : array_like
        Latitudes and longitudes along each track
    vstore : array_like
        Maximum circular wind along each storm track
    rmstore : array_like
        Radius (in km) of maximum circular wind along each track
    vsestore : array_like
        Maximum circular wind of any secondary eyewalls that may be present
    rmsestore : array_like
        Radius (in km) of maximum circular wind of any secondary eyewalls
    ut, vt : array_like
        West-east and north-south components of the storm translation velocity
    u850store, v850store : array_like
        Zonal & meridional components of the 850 hPa environmental wind speed (knots)

    Returns:
    --------
    x, y : array_like
        Vectors containing the longitudes and latitudes of the grid
    netrain : array_like
        Storm total rainfall (unit: mm) at each point on the grid
    """
    # Load parameters
    magfac = params.magfac
    deltax, deltay = params.deltax, params.deltay
    bxmin, bxmax = params.bxmin, params.bxmax
    bymin, bymax = params.bymin, params.bymax
    dellatlongs = params.dellatlongs
    q900 = params.q900
    timeres = params.timeres
    wrad = params.wrad
    eprecip = params.eprecip

    # Constants
    M_TO_MM = 1000
    RHOA_OVER_RHOL = 0.00117  # rho_air / rho_liquid

    # Initialize variables
    nrm, mrm = np.shape(rmstore)
    rfac = magfac * np.ones((nrm, mrm))

    bxmin = (bxmin + 360) if bxmin < 0 else bxmin
    bxmax = (bxmax + 360) if bxmax < 0 else bxmax

    latdata = latstore[nt, :]
    latdata = latdata[(latdata != 0) & ~np.isnan(latdata)]
    latsize = len(latdata)

    utd = ut[nt, :latsize].reshape((1, latsize))
    vtd = vt[nt, :latsize].reshape((1, latsize))
    ush = np.zeros_like(utd)
    vsh = np.zeros_like(vtd)

    vdrift = 1.5 * 3600 / 1852 * latstore[0, 0] / (np.abs(latstore[0, 0]) + 1e-8)

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

    # Scale radii of maximum wind
    rm = (magfac * rfac[nt, :] * rmstore[nt, :])[:latsize].reshape((1, latsize))
    rmse = (magfac * rfac[nt, :] * rmsestore[nt, :])[:latsize].reshape((1, latsize))

    # Adjust longitudes for date line crossing
    long[0, (long[0, 0] > 200) & (long[0, :] < 50)] += 360

    # Calculate map boundaries
    bxmin = np.min(long[np.nonzero(long)]) - deltax
    bxmax = np.max(long[np.nonzero(long)]) + deltax
    bymin = np.min(lat[np.nonzero(lat)]) - deltay
    bymax = np.max(lat[np.nonzero(lat)]) + deltay

    # Get topographic height and its gradients
    h, hx, hy, x, y = tcr_tb.estimate_topographic_height(bxmin, bxmax, bymin, bymax, dellatlongs)

    # Calculate vertical velocity time series
    w = tcr_wind.pointwshortn(lat, long, v, rm, vse, rmse, utd, vtd,
                              ush, vsh, x, y, h, hx, hy, timeres)
    wq = np.maximum(w - wrad, 0) * q900
    netrain = eprecip * M_TO_MM * timeres * 3600 * RHOA_OVER_RHOL * np.sum(wq, axis=(0, 1))

    return x, y, netrain.T


def raingen(plat, plong, latstore, longstore, datestore, vstore, rmstore, vsestore, rmsestore,
            u850store, v850store, utrans, vtrans, T600=None):
    """
    Calculate accumulated rain and rainrates at specified locations for the active event set.

    Parameters:
    -----------
    plat : array-like
        Latitude of points of interest.
    plong : array-like
        Longitude of points of interest.
    latstore : array-like
        Latitudes along each storm track.
    longstore : array-like
        Longitudes along each storm track.
    datestore : array-like
        Datetime info in integer format.
    vstore : array-like
        Maximum circular wind along each storm track.
    rmstore : array-like
        Radius (in km) of maximum circular wind along each track.
    vsestore : array-like
        Maximum circular wind of any secondary eyewalls that may be present.
    rmsestore : array-like
        Radius (in km) of maximum circular wind of any secondary eyewalls.
    u850store : array-like
        Zonal component of the 850 hPa environmental wind speed (in knots).
    v850store : array-like
        Meridional component of the 850 hPa environmental wind speed (in knots).
    utrans : array-like
        West-east component of the storm translation velocity.
    vtrans : array-like
        North-south component of the storm translation velocity.
    T600 : array-like, optional
        Temperature at 600 hPa level.

    Returns:
    --------
    rain : ndarray
        Total storm rainfall (unit: mm).
    rainrate : ndarray
        Rain rate (unit: mm/hr).
    date_record : ndarray
        Time in date format corresponding to rainrate.
    """
    ut = np.nan_to_num(utrans)
    vt = np.nan_to_num(vtrans)

    # Convert point lat,lon to array
    plat = np.array([plat])
    plong = np.array([plong])

    # Load parameters
    magfac = params.magfac                  # overall scale factor for storm size
    q900_default = params.q900              # specific humidity at 900 hPa
    timeres = params.timeres                # time resolution for time series at fixed points
    wrad = params.wrad                      # background subsidence velocity under radiative cooling
    eprecip = params.eprecip                # precipitation efficiency

    sx = plong.size                         # get number of points

    # Load high-resolution bathymetry from netcdf
    bathy = tcr_io.load_netcdf_2d_parameters('../data', 'surface_data.nc', 'bathymetry')

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
    vdrift = 1.5 * 3600 / 1852 * latstore[0, 0] / (np.abs(latstore[0, 0]) + 1e-8)

    if 'u850store' in locals():
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

    # Constants
    M_TO_MM = 1000
    RHOA_OVER_RHOL = 0.00117  # rho_air / rho_liquid

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
    wq, date_record = tcr_wind.pointwshortnqdx(
        lat, long, dates, q900, v, rm, vse, rmse, ut, vt, ush, vsh,
        plong, plat, h, hx, hy, timeres, wrad)

    # Convert date_record to pandas format
    datetimes = pd.to_datetime(date_record.flatten(), unit='s')
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
    date_record = np.reshape(datetimes_numpy, wq.shape)

    rainrate = eprecip * M_TO_MM * 3600 * RHOA_OVER_RHOL * wq
    rainrate = np.nan_to_num(rainrate)
    rain = timeres * np.sum(rainrate, axis=1).reshape(-1, 1)

    return rain, rainrate, date_record
