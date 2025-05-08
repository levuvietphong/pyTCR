"""
Functions for physics in PyTCR
"""

import numpy as np


def estimate_track_density(lat_trks, lon_trks, vmax_trks, num_trks, threshold, cellsize, interval):
    """
    Estimate the total number of tropical cyclone (TC) tracks crossing a latitude-longitude box
    at a specified resolution (cellsize x cellsize in degree).

    Parameters:
    -----------
    lat_trks : numpy.ndarray
        Array of latitudes for each TC track (degree)
    lon_trks : numpy.ndarray
        Array of longitudes for each TC track (degree)
    vmax_trks : numpy.ndarray
        Array of maximum wind speeds for each TC track (m/s)
    num_trks : int
        Number of TC tracks to randomly sample
    threshold : float
        Wind speed threshold to consider a TC track (m/s)
    cellsize : float
        Size of the grid cell (degree)
    interval : int
        Interval for sampling points along each track.

    Returns:
    --------
    latitude : numpy.ndarray
        Array of latitude values for the grid (degree)
    longitude : numpy.ndarray
        Array of longitude values for the grid (degree)
    density_all : numpy.ndarray
        2D array representing the density of TC tracks in each grid cell
    """

    num_tc = lat_trks.shape[0]
    ind_trks = np.random.choice(np.arange(num_tc), num_trks, replace=False)

    x0, x1 = cellsize / 2, 360 - cellsize / 2
    y0, y1 = cellsize / 2, 180 - cellsize / 2

    nrows, ncols = int(180 / cellsize), int(360 / cellsize)
    latitude = np.arange(y0, y1 + cellsize, cellsize)
    longitude = np.arange(x0, x1 + cellsize, cellsize)
    density_all = np.zeros((nrows, ncols))

    for tc_id in ind_trks:
        density = np.zeros((nrows, ncols))
        temp = lat_trks[tc_id, :]
        indmax = np.where(temp > -90)[0][-1] + 1
        lat = lat_trks[tc_id, 0:indmax:interval]
        lon = lon_trks[tc_id, 0:indmax:interval]
        lon[lon >= 360] -= 360

        vmax = vmax_trks[tc_id, 0:indmax:interval]
        valid_indices = np.where(vmax > threshold)
        lat, lon, vmax = lat[valid_indices], lon[valid_indices], vmax[valid_indices]

        jrow = np.floor(lat / cellsize).astype(int)
        icol = np.floor(lon / cellsize).astype(int)
        density[jrow, icol] = 1
        density_all += density

    return latitude, longitude, density_all


def estimate_pdi(lat_trks, lon_trks, vmax_trks, num_trks, cellsize, dt):
    """
    Estimate the Power Dissipation Index (PDI) of tropical cyclones.

    Parameters:
    -----------
    lat_trks : numpy.ndarray
        Array of latitude tracks for tropical cyclones (degree)
    lon_trks : numpy.ndarray
        Array of longitude tracks for tropical cyclones (degree)
    vmax_trks : numpy.ndarray
        Array of maximum wind speeds for tropical cyclones (m/s)
    num_trks : int
        Number of TC tracks to randomly sample.
    cellsize : float
        Size of the grid cell (degrees)
    dt : float
        Time interval in hours.

    Returns:
    --------
    latitude : numpy.ndarray
        Array of latitude values for the grid (degree)
    longitude : numpy.ndarray
        Array of longitude values for the grid (degree)
    pdi_all : numpy.ndarray
        2D array representing the PDI of TC tracks in each grid cell
        over a period of time (unit m^3/s^2)
    """
    num_tc = lat_trks.shape[0]
    ind_trks = np.random.choice(num_tc, num_trks, replace=False)

    x0, x1 = cellsize / 2, 360 - cellsize / 2
    y0, y1 = cellsize / 2, 180 - cellsize / 2

    nrows, ncols = int(180 / cellsize), int(360 / cellsize)
    latitude = np.arange(y0, y1 + cellsize, cellsize)
    longitude = np.arange(x0, x1 + cellsize, cellsize)
    pdi_all = np.zeros((nrows, ncols))

    for tc_id in ind_trks:
        temp = lat_trks[tc_id, :]
        indmax = np.where(temp > -90)[0][-1] + 1
        lat = lat_trks[tc_id, :indmax]
        lon = lon_trks[tc_id, :indmax]
        vmax = vmax_trks[tc_id, :indmax]
        lon[lon >= 360] -= 360

        valid_indices = ~np.isnan(vmax) & (vmax >= 0)
        lat, lon = lat[valid_indices], lon[valid_indices]
        vmax = vmax[valid_indices]

        jrows = np.floor(lat / cellsize).astype(int)
        icols = np.floor(lon / cellsize).astype(int)
        np.add.at(pdi_all, (jrows, icols), vmax**3 * dt * 3600)

    return latitude, longitude, pdi_all
