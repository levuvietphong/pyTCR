"""
Functions for rainfall in PyTCR
"""

import datetime
import pandas as pd
import numpy as np
from tcr import wind as tcr_wind
from tcr import terrain_boundary as tcr_tb
from tcr import iodata as tcr_io
from tcr import parameters as tcr_params
from tcr.datadir import BASE_DATA_DIR, DOWNSCALED_DATA_DIR


def calculate_rainfall_rate(
    nt, latitude, longitude, radius_storm, velocity, radius_storm_secondary,
    velocity_secondary, ut, vt, u850, v850, months, days, hours, monthplot,
    dayplot, hourplot, extent=None, shapefile=None, magfac=tcr_params.magfac,
    deltax=tcr_params.deltax, deltay=tcr_params.deltay, dellatlong=tcr_params.dellatlong,
    q900=tcr_params.q900, eprecip=tcr_params.eprecip, wrad=tcr_params.wrad
):
    """
    Computes the surface rain rate distribution (mm/hr) for a specified storm
    at a given time.

    Parameters:
    -----------
        nt : int
            Index representing the storm track number
        latitude : ndarray
            Latitude values for each time step and storm (degree)
        longitude : ndarray
            Longitude values for each time step and storm (degree)
        radius_storm : ndarray
            Radii of maximum wind for each time step and storm (km)
        velocity : ndarray
            Maximum circular wind speed for each time step and storm (knots)
        raidus_storm_secondary : ndarray
            Secondary radii of maximum wind for each time step and storm (km)
        velocity_secondary : ndarray
            Secondary maximum circular wind speed for each time step and storm (knots)
        ut : ndarray
            West-east component of the storm translation velocity (knots)
        vt : ndarray
            North-south component of the storm translation velocity (knots)
        u850 : ndarray
            U-component wind speed at 850 hPa for each time step and storm (knots)
        v850 : ndarray
            V-component wind speed at 850 hPa for each time step and storm (knots)
        months : ndarray
            Month values (1-12) for each time step.
        days : ndarray
            Day values for each time step
        hours : ndarray
            Hour values for each time step
        monthplot : int
            The month (1-12) for which the rain rate is to be calculated
        dayplot : int
            The day of the month (1-31) for which the rain rate is to be calculated
        hourplot : int
            The hour of the day (0-23, GMT) for which the rain rate is to be calculated
        extent : tuple, optional
            Bounding box and spacing in data coordinates (left, right, bottom, top)
            defines the spatial extent of the map, default is None
        shapefile : str or shapefile-like object, optional
            A shapefile to overlay on the map. Provides additional geographic context
            Default is None
        magfac : float, optional
            Overall scale factor for storm size, default is 1.0
        deltax : float, optional
            Longitudinal distance of map boundaries from storm center (degree)
            Default is 5
        deltay : float, optional
            Latitudinal distance of map boundaries from storm center (degree)
            Default is 4
        dellatlong : float, optional
            Horizontal resolution of field maps (degree), default is 0.15
        q900 : float, optional
            Specific humidity at 900 hPa if T600 is unavailable (g/g)
            Default is 0.01
        eprecip : float, optional
            Precipitation efficiency, default is 0.9
        wrad : float, optional
            Background subsidence velocity under radiative cooling (m/s)
            Default is -0.005

    Returns:
    --------
        rainrate : ndarray
            Surface rainfall rates for each storm track (mm/hr)

    """

    if extent is None and shapefile is None:
        bxmin, bxmax, bymin, bymax = 20, 380, -60, 60
    else:
        if shapefile is None:
            bxmin, bxmax, bymin, bymax = extent
        else:
            bxmin, bxmax, bymin, bymax = tcr_io.get_bbox_from_shapefile(
                shapefile)

    rhoa_over_rhol = tcr_params.RHOA_OVER_RHOL  # ratio of air density over water density
    pifac = tcr_params.RAD2DEG  # convert radians to degrees (pi/180)
    knotfac = tcr_params.KNOTS2MPS  # convert knots to m/s (1 knots = 0.5144 m/s)    
    nrm, mrm = np.shape(radius_storm)
    rfac = magfac * (1+np.zeros((nrm, mrm)))

    # Get the length of each event by finding the first 0 element,
    # if non-zero, get all length
    jmaxd = np.argmin(velocity, axis=1)
    duration = jmaxd[nt]
    if duration == 0:
        duration = velocity.shape[1]

    dum = np.zeros((duration))
    V = np.column_stack(
        (
            dum + 2000,
            months[nt, 0:duration].T,
            days[nt, 0:duration].T,
            hours[nt, 0:duration].T,
            dum,
            dum,
        ))
    timev = pd.to_datetime(
        dict(
            year=V[:, 0],
            month=V[:, 1],
            day=V[:, 2],
            hour=V[:, 3],
            minute=V[:, 4],
            second=V[:, 5],
        ))
    timeplot = datetime.datetime(2000, monthplot, dayplot, hourplot, 0, 0)

    diff0 = timeplot - timeplot
    timediff = timeplot - timev
    timediff[timediff < diff0] = diff0
    jplot = np.argmin(timediff)
    jplot = np.maximum(jplot, 1)
    jplot = np.minimum(jplot, 398)

    latstorm = latitude[nt, jplot - 1:jplot + 2]
    longstorm = longitude[nt, jplot - 1:jplot + 2]
    longstorm[longstorm < 0] += 360
    vstorm = velocity[nt, jplot - 1:jplot + 2]
    rmstorm = rfac[nt, jplot - 1:jplot + 2] * radius_storm[nt, jplot - 1:jplot + 2]
    vsestorm = velocity_secondary[nt, jplot - 1:jplot + 2]
    rmsestorm = rfac[nt, jplot - 1:jplot + 2] * radius_storm_secondary[
        nt, jplot - 1:jplot + 2]
    utstorm = ut[nt, jplot - 1:jplot + 2]
    vtstorm = vt[nt, jplot - 1:jplot + 2]
    ush = np.zeros(3)
    vsh = np.zeros(3)

    if u850 is not None:
        vdrift = 1.5 * 3600 / 1852
        vdrift = vdrift * latitude[0, 0] / (abs(latitude[0, 0]) + 1e-8)
        u850storm = u850[nt, jplot - 1:jplot + 2]
        v850storm = v850[nt, jplot - 1:jplot + 2]
        ush = 5 * knotfac * (utstorm - u850storm)
        vsh = 5 * knotfac * (vtstorm - vdrift * np.cos(pifac * latstorm) - v850storm)

    if extent is None and shapefile is None:
        bxmin = np.floor(longstorm[1] - deltax)
        bxmax = np.ceil(longstorm[1] + deltax)
        bymin = np.floor(latstorm[1] - deltay)
        bymax = np.ceil(latstorm[1] + deltay)

    h, hx, hy, x, y = tcr_tb.estimate_topographic_height(
        bxmin, bxmax, bymin, bymax, dellatlong
    )

    w = tcr_wind.calculate_upward_velocity_field(
        latstorm, longstorm, vstorm, rmstorm, vsestorm, rmsestorm, utstorm,
        vtstorm, ush, vsh, y, x, h, hx, hy
    )

    temp = (
        eprecip * 1000 * 3600 * rhoa_over_rhol * q900 *
        np.maximum(w[0, 1, :, :] + wrad, 0)
    )

    rainrate = temp
    rainrate = rainrate.transpose()

    return rainrate, x, y


def calculate_etr_swath(
    nt, latitude, longitude, radius_storm, velocity, radius_storm_secondary,
    velocity_secondary, ut, vt, u850, v850, extent=None, shapefile=None,
    magfac=1, deltax=tcr_params.deltax, deltay=tcr_params.deltay,
    dellatlongs=tcr_params.dellatlongs, q900=tcr_params.q900, timeres=tcr_params.timeres,
    wrad=tcr_params.wrad, eprecip=tcr_params.eprecip
):
    """
    Calculate the distribution of event total rainfall for a given individual
    storm.

    Parameters:
    -----------
    nt : int
        Track number of the storm
    latitude, longitude : array_like
        Latitudes and longitudes along each track (degree)
    velocity : array_like
        Maximum circular wind along each storm track (m/s)
    radius_storm : array_like
        Radius of maximum circular wind along each track (km)
    velocity_secondary : array_like
        Maximum circular wind of any secondary eyewalls that may be present (m/s)
    radius_storm_secondary : array_like
        Radius of maximum circular wind of any secondary eyewalls (km)
    ut, vt : array_like
        West-east and north-south components of the storm translation velocity (knots)
    u850, v850 : array_like
        Zonal & meridional components of the 850 hPa env wind speed (knots)
    extent : tuple, optional
        Bounding box and spacing in data coordinates (left, right, bottom, top, dx, dy)
        Defines the spatial extent of the map, default is None
    shapefile : str or shapefile-like object, optional
        A shapefile to overlay on the map. Provides additional geographic context
    magfac : float, optional
        Overall scale factor for storm size, default is 1.0
    deltax : float, optional
        Longitudinal distance of map boundaries from storm center (degree)
        Default is 5
    deltay : float, optional
        Latitudinal distance of map boundaries from storm center (degree)
        Default is 4
    dellatlongs : float, optional
        Horizontal resolution of swath maps (degree), default is 0.15
    q900 : float, optional
        Specific humidity at 900 hPa if T600 is unavailable (g/g)
        Default is 0.01
    timeres : float, optional
        Time resolution for time series at fixed points (hours), default is 2
    wrad : float, optional
        Background subsidence velocity under radiative cooling (m/s)
        Default is -0.005
    eprecip : float, optional
        Precipitation efficiency, default is 0.9

    Returns:
    --------
    x, y : array_like
        Vectors containing the longitudes and latitudes of the grid (degree)
    netrain : array_like
        Storm total rainfall at each point on the grid (mm)
    """

    if extent is None and shapefile is None:
        bxmin, bxmax, bymin, bymax = 20, 380, -60, 60
    else:
        if shapefile is None:
            bxmin, bxmax, bymin, bymax = extent
        else:
            bxmin, bxmax, bymin, bymax = tcr_io.get_bbox_from_shapefile(
                shapefile)

    # Constants
    rhoa_over_rhol = tcr_params.RHOA_OVER_RHOL  # rho_air / rho_liquid

    # Initialize variables
    nrm, mrm = np.shape(radius_storm)
    rfac = magfac * np.ones((nrm, mrm))

    bxmin = (bxmin + 360) if bxmin < 0 else bxmin
    bxmax = (bxmax + 360) if bxmax < 0 else bxmax

    latdata = latitude[nt, :].copy()
    latdata = latdata[(latdata != 0) & ~np.isnan(latdata)]
    latsize = len(latdata)
    # latsize = max(latsize, 193)

    utd = ut[nt, :latsize].copy().reshape((1, latsize))  # (knots)
    vtd = vt[nt, :latsize].copy().reshape((1, latsize))  # (knots)
    ush = np.zeros_like(utd)
    vsh = np.zeros_like(vtd)

    vdrift = 1.5 * 3600 / 1852 * latitude[0, 0] / (np.abs(latitude[0, 0]) + 1e-8)  # convert to (m/s)

    if u850 is not None:
        ush = 5 * 1852 / 3600 * (utd - u850[nt, :latsize])  # (m/s)
        vsh = (
            5 * 1852 / 3600
            * (
                vtd
                - vdrift * np.cos(np.pi / 180 * latitude[nt, :latsize])
                - v850[nt, :latsize]
            )
        )  # (m/s)

    lat = latitude[nt, :latsize].copy().reshape((1, latsize))
    long = longitude[nt, :latsize].copy().reshape((1, latsize))

    # Convert longitude to 0 to 360 degree east
    long[long < 0] += 360
    
    v = velocity[nt, :latsize].copy().reshape((1, latsize))  # (m/s)
    vse = velocity_secondary[nt, :latsize].copy().reshape((1, latsize))  # (m/s)

    # Scale radii of maximum wind (km)
    rm = (magfac * rfac[nt, :] * radius_storm[nt, :])[:latsize].reshape(
        (1, latsize)
    )  # (km)
    rmse = (magfac * rfac[nt, :] * radius_storm_secondary[nt, :])[
        :latsize].reshape((1, latsize))  # (km)

    # Adjust longitudes for date line crossing
    long[0, (long[0, 0] > 200) & (long[0, :] < 50) & (long[0, :] != 0)] += 360

    # Calculate map boundaries (degree) if not specified
    if extent is None and shapefile is None:
        bxmin = np.min(long[np.nonzero(long)]) - deltax
        bxmax = np.max(long[np.nonzero(long)]) + deltax
        bymin = np.min(lat[np.nonzero(lat)]) - deltay
        bymax = np.max(lat[np.nonzero(lat)]) + deltay

    # Get topographic height (m) and its gradients
    h, hx, hy, x, y = tcr_tb.estimate_topographic_height(
        bxmin, bxmax, bymin, bymax, dellatlongs
    )

    # Calculate vertical velocity time series (m/s)
    w = tcr_wind.calculate_upward_velocity_time_series(
        lat, long, v, rm, vse, rmse, utd, vtd, ush, vsh, x, y, h,
        hx, hy, timeres
    )
    wq = np.maximum(w + wrad, 0) * q900

    # Calculate event total rainfall (mm)
    netrain = (eprecip * 1000 * timeres * 3600 * rhoa_over_rhol *
               np.nansum(wq, axis=(0, 1)))

    return x, y, netrain.T


def generate_rainfall_point(
    plat, plong, latitude, longitude, datearray, velocity,
    radius_storm, velocity_secondary, radius_storm_secondary,
    u850, v850, utrans, vtrans, T600=None, magfac=tcr_params.magfac,
    q900_constant=tcr_params.q900, timeres=tcr_params.timeres,
    wrad=tcr_params.wrad, eprecip=tcr_params.eprecip
):
    """
    Calculate the accumulated rainfall and rain rates at a specified location
    for the given set of storm events.

    Parameters:
    -----------
    plat : array-like
        Latitudes of the points of interest (degree)
    plong : array-like
        Longitudes of the points of interest (degree)
    latitude : array-like
        Latitudes along each storm track (degree)
    longitude : array-like
        Longitudes along each storm track (degree)
    datearray : array-like
        Datetime information in integer format
    velocity : array-like
        Maximum circular wind speeds along each storm track (m/s)
    radius_storm : array-like
        Radii of maximum circular wind along each track (km)
    velocity_secondary : array-like
        Maximum circular wind speeds of any secondary eyewalls that may be
        present (m/s)
    radius_storm_secondary : array-like
        Radii of maximum circular wind of any secondary eyewalls (km)
    u850 : array-like
        Zonal component of the 850 hPa environmental wind speed (m/s)
    v850 : array-like
        Meridional component of the 850 hPa environmental wind speed (m/s)
    utrans : array-like
        West-east component of the storm translation velocity (knots)
    vtrans : array-like
        North-south component of the storm translation velocity (knots)
    T600 : array-like, optional
        Temperature at the 600 hPa level (K)
    magfac : float, optional
        Overall scale factor for storm size, default is 1.0
    q900_constant : float, optional
        Default specific humidity at 900 hPa if T600 is unavailable (g/g)
        Default is 0.01
    timeres : float, optional
        Time resolution for time series at fixed points (hours), default is 0.5
    wrad : float, optional
        Background subsidence velocity under radiative cooling (m/s)
        Default is -0.005
    eprecip : float, optional
        Precipitation efficiency, default is 0.9

    Returns:
    --------
    rain : ndarray
        Event total rainfall (mm)
    rainrate : ndarray
        Rate of rainfall (mm/hr)
    date_record : ndarray
        Time in date format corresponding to the rain rate.
    """

    ut = np.nan_to_num(utrans)
    vt = np.nan_to_num(vtrans)

    # Convert point lat,lon to array
    plat = np.array([plat])
    plong = np.array([plong])
    sx = plong.size  # get number of points

    # Load high-resolution bathymetry from netcdf
    bathy = tcr_io.load_netcdf_2d_parameters(
        BASE_DATA_DIR, 'surface_data.nc', 'bathymetry_high'
    )

    ntopo, _ = np.shape(bathy)  # get number of grid points in longitude direction
    topores = 360.0 / ntopo     # topo resolution in degree
    toporesi = 1.0 / topores    # inverse of topo resolution
    sfac = 1.0 / (tcr_params.DEG2KM * topores * 1000)  # topo resolution in degree to meter
    pifac = tcr_params.RAD2DEG   # convert radians to degrees (pi/180)
    knotfac = tcr_params.KNOTS2MPS  # convert nautical mile to m/s
    m, n = ut.shape  # m: num of storm, n: num of time steps

    if np.min(plong) < 0:
        plong += 360

    ush = np.zeros((n, m))
    vsh = np.zeros((n, m))
    vdrift = 1.5 * 3600 / 1852 * latitude[0, 0] / (np.abs(latitude[0, 0]) + 1e-8)

    # convert to unit m/s
    if u850 is not None:
        ush = 5 * knotfac * (ut - u850)
        vsh = 5 * knotfac * (vt - vdrift * np.cos(pifac * latitude) - v850)

    lat = latitude.copy()
    long = longitude.copy()
    dates = datearray.copy()
    v = velocity.copy()
    vse = velocity_secondary.copy()

    nrm, mrm = radius_storm.shape
    rfac = magfac * np.ones((nrm, mrm))

    rm = radius_storm * rfac
    rmse = radius_storm_secondary * rfac

    bathy = np.maximum(bathy, -1)

    # Calculate topographic and its gradients
    h, hx, hy = tcr_tb.calculate_spatial_derivatives(
        bathy, plong, plat, sx, 1, sfac, pifac, ntopo, toporesi)

    # Calculates saturation specific humidity
    if T600 is not None:
        q900, _ = tcr_tb.calculate_qs900(T600, velocity)
    else:
        q900 = np.where(ut != 0, q900_constant, 0)

    # Estimate vertical wind velocity
    wq, date_record = tcr_wind.calculate_upward_velocity_time_series(
        lat, long, v, rm, vse, rmse, ut, vt, ush, vsh, plong, plat,
        h, hx, hy, timeres, date_records=dates, dq=q900, wrad=wrad)

    # Convert date_record to pandas format
    datetimes = pd.to_datetime(np.where(date_record.flatten() < 0,
                                        np.nan, date_record.flatten()),
                               unit='s')
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

    # Constants
    rhoa_over_rhol = tcr_params.RHOA_OVER_RHOL  # rho_air / rho_liquid ratio

    # Calculate rain rate (mm/hr)
    rainrate = eprecip * 1000 * 3600 * rhoa_over_rhol * wq
    rainrate = np.nan_to_num(rainrate)

    # Total event rainfall (mm) 
    rain = timeres * np.sum(rainrate, axis=1).reshape(-1, 1)  # (mm)
    rain = rain.flatten()

    return rain, rainrate, date_record
