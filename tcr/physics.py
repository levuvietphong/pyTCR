"""
Functions for physics in PyTCR
"""

import numpy as np


def estimate_track_density(lat_trks, lon_trks, vmax_trks, num_trks, threshold, cellsize, interval):
    """
    Estimate the total number of TC tracks crossing a lat-long box at resolution
    cellsize x cellsize (degree).
    """
    num_tc = np.shape(lat_trks)[0]
    ind_trks = np.random.choice(np.arange(num_tc), num_trks, replace=False)
    num_tc = len(ind_trks)

    x0 = cellsize / 2
    x1 = 360 - cellsize / 2
    y0 = 0 + cellsize / 2
    y1 = 180 - cellsize / 2

    nrows = int(180. / cellsize)
    ncols = int(360. / cellsize)
    latitude = np.arange(y0, y1+cellsize, cellsize)
    longitude = np.arange(x0, x1+cellsize, cellsize)
    density_all = np.zeros((nrows, ncols))

    for tc_id in ind_trks:
        density = np.zeros((nrows, ncols))
        temp = lat_trks[tc_id, :]
        indmax = np.where(temp > -90)[0][-1]+1
        lat = lat_trks[tc_id, 0:indmax:interval]
        lon = lon_trks[tc_id, 0:indmax:interval]
        ind = np.where(lon >= 360)

        if len(ind) > 0:
            lon[ind] -= 360

        vmax = vmax_trks[tc_id, 0:indmax:interval]
        ind = np.where(vmax > threshold)
        lat = lat[ind]
        lon = lon[ind]
        vmax = vmax[ind]

        jrow = np.floor(lat / cellsize).astype(int)
        icol = np.floor(lon / cellsize).astype(int)
        density[jrow, icol] = 1
        density_all += density

    return latitude, longitude, density_all


def estimate_pdi(lat_trks, lon_trks, vmax_trks, num_trks, cellsize, dt):
    """
    Estimate the Power Dissipitation Index (PDI) of tropical cyclones
    """
    num_tc = np.shape(lat_trks)[0]
    ind_trks = np.random.choice(np.arange(num_tc), num_trks, replace=False)
    num_tc = len(ind_trks)
    x0 = cellsize/2
    x1 = 360-cellsize/2
    y0 = 0+cellsize/2
    y1 = 180-cellsize/2

    nrows = int(180/cellsize)
    ncols = int(360/cellsize)
    latitude = np.arange(y0, y1+cellsize, cellsize)
    longitude = np.arange(x0, x1+cellsize, cellsize)
    pdi_all = np.zeros((nrows, ncols))

    for tc_id in ind_trks:
        pdi = np.zeros((nrows, ncols))
        temp = lat_trks[tc_id, :]
        indmax = np.where(temp > -90)[0][-1]+1
        lat = lat_trks[tc_id, 0:indmax]
        lon = lon_trks[tc_id, 0:indmax]
        vmax = vmax_trks[tc_id, 0:indmax]
        ind = np.where(lon >= 360)
        if len(ind) > 0:
            lon[ind] -= 360

        for k in range(indmax):
            jrow = np.floor(lat[k] / cellsize).astype(int)
            icol = np.floor(lon[k] / cellsize).astype(int)
            if (np.isnan(vmax[k])) or (vmax[k] < 0):
                vmax[k] = 0
            pdi[jrow, icol] += vmax[k]**3 * dt * 3600         # unit L^3/L^2
        pdi_all += pdi

    return latitude, longitude, pdi_all
