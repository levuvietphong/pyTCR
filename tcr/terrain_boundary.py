"""
Functions for topography and boundary layers in PyTCR
"""

import math
import numpy as np
from tcr import iodata as tcr_io


def calculate_distance_POI_from_track(plat, plong, latstore, longstore, nn, m, sx, sy, ngrid, dfac):
    """
    This function calculates the distances of points of interest from the track

    Parameters:
        - plat, plong: point latitude & longitude
        - latstore, longstore: Latitude and longitude of points along each track
        - nn: number of storms
        - m: length of each storm
        - sx, sy: number of points for x & y direction?
        - ngrid: number of grid points???
        - dfac: factor converting degree to km
    Returns:
        radius: Euclidean distance from track
        dx, dy: distance components in x & y directions
    """
    pifac = math.acos(-1)/180                   # pi number
    dx = np.zeros((nn, m, sx, sy))              # 4d array of dx
    dy = np.zeros((nn, m, sx, sy))              # 4d array of dy
    for i in range(sx):
        for jj in range(sy):
            j = i if ngrid == 1 else jj

            if (np.ndim(longstore) == 1) & (np.ndim(latstore) == 1):
                dx[:, :, i, j] = dfac * np.cos(pifac*plat[j]) * (plong[i]-longstore)
                dy[:, :, i, j] = dfac * (plat[j]-latstore)
            elif (np.ndim(longstore) == 2) & (np.ndim(latstore) == 2):
                dx[:, :, i, j] = dfac * np.cos(pifac*plat[j]) * (plong[i]-longstore)
                dy[:, :, i, j] = dfac * (plat[j]-latstore)
            elif (np.ndim(longstore) == 4) & (np.ndim(latstore) == 4):
                dx[:, :, i, j] = dfac * np.cos(pifac*plat[j]) * (plong[i]-longstore[:, :, i, jj])
                dy[:, :, i, j] = dfac * (plat[j]-latstore[:, :, i, jj])

    radius = np.sqrt(dx * dx + dy * dy)
    return radius, dx, dy


def calculate_spatial_derivatives(bathy, x, y, sx, sy, sfac, pifac, ntopo, toporesi):
    """
        This function computes the spatial derivatives of a given topography

    Parameters:
    ----------
        - bathy : bathymetry
        - x : spatial coordinate in x
        - y : spatial coordinate in y
        - sx : size of x
        - sy : size of y
        - sfac : factor converting nautical to km
        - pifac : pi number
        - ntopo : number of row of bathymetry array
        - toporesi : inverse of topgraphic resolution (topores)

    Returns:
    -------
        - h : topographic heights
        - hx, hy : derivatives of h in x and y
    """

    h = np.zeros((sx, sy))
    hx = np.zeros((sx, sy))
    hy = np.zeros((sx, sy))
    bathy = np.maximum(bathy, -1)
    dhdx = sfac * (np.roll(bathy, -1, 0) - np.roll(bathy, 0, 0))
    dhdy = sfac * (np.roll(bathy, -1, 1) - np.roll(bathy, 0, 1))

    for i in range(sx):
        plong = x[i]
        if plong >= 360:
            plong -= 360

        for j in range(sy):
            plat = y[j]
            ib = np.floor(toporesi * plong).astype(int)
            ibp = ib + 1
            if ibp >= ntopo:
                ibp = ibp - ntopo

            ibs = np.floor(toporesi * plong - 0.5).astype(int)
            ibsp = ibs + 1
            plongs = plong

            if ibs < -1:
                ibs = ntopo - 1
                plongs = plong + 360

            if ibsp >= ntopo:
                ibsp = ibsp - ntopo

            jb = np.floor(toporesi * (plat + 90)).astype(int)
            jbs = np.floor(toporesi * (plat + 90) - 0.5).astype(int)
            b1 = bathy[ib, jb]
            b2 = bathy[ib, jb + 1]
            b3 = bathy[ibp, jb]
            b4 = bathy[ibp, jb + 1]
            dely = toporesi * (plat + 90) - jb
            delx = toporesi * plong - ib
            d1 = (1 - delx) * (1 - dely)
            d2 = dely * (1 - delx)
            d3 = delx * (1 - dely)
            d4 = delx * dely
            h[i, j] = (
                np.exp(
                    d1 * np.log(b1 + 11)
                    + d2 * np.log(b2 + 11)
                    + d3 * np.log(b3 + 11)
                    + d4 * np.log(b4 + 11)
                )
                - 11
            )

            b1 = dhdx[ibs, jbs]
            b2 = dhdx[ibs, jbs + 1]
            b3 = dhdx[ibsp, jbs]
            b4 = dhdx[ibsp, jbs + 1]
            dely = -0.5 + toporesi * (plat + 90) - jbs
            delx = -0.5 + toporesi * plongs - ibs
            d1 = (1 - delx) * (1 - dely)
            d2 = dely * (1 - delx)
            d3 = delx * (1 - dely)
            d4 = delx * dely
            hx[i, j] = (b1 * d1 + b2 * d2 + b3 * d3 + b4 * d4) / np.cos(pifac * plat)

            b1 = dhdy[ibs, jbs]
            b2 = dhdy[ibs, jbs + 1]
            b3 = dhdy[ibsp, jbs]
            b4 = dhdy[ibsp, jbs + 1]
            hy[i, j] = b1 * d1 + b2 * d2 + b3 * d3 + b4 * d4
    return h, hx, hy


def estimate_topographic_height(bxmin, bxmax, bymin, bymax, dellatlong):
    """
    Load topographic from file and estimate the spatial derivatives

    Inputs:
    -------
        - bxmin: bound min in x-axis
        - bxmax: bound max in x-axis
        - bymin: bound min in y-axis
        - bymax: bound max in y-axis
        - dellatlong: horizontal resolution of map (unit: degree)

    Returns:
    --------
        - h: topographic data
        - hx, hy: topographic gradient in x & y direction
        - x, y: vectors containing the longs and lats of the grid
    """

    # Load high-resolution bathymetry file in `data`` directory
    bathy = tcr_io.load_netcdf_2d_parameters('../data', 'surface_data.nc', 'bathymetry_high')
    ntopo, _ = np.shape(bathy)                      # number of grid point in x direction
    topores = 360 / ntopo                           # topo resolution in degree
    toporesi = 1 / topores                          # inverse of topo resolution
    sfac = 1.0 / (topores * 60.0 * 1852)            # factor converting degree to m
    pifac = math.acos(-1) / 180                     # pi number

    x = np.round(np.arange(bxmin, bxmax + 1e-8, dellatlong), 4)
    y = np.round(np.arange(bymin, bymax + 1e-8, dellatlong), 4)
    sx = np.max(np.shape(x))
    sy = np.max(np.shape(y))

    if sx == sy:
        sx += 1
        x = np.concatenate([x, [bxmax+dellatlong]])

    h, hx, hy = calculate_spatial_derivatives(
        bathy, x, y, sx, sy, sfac, pifac, ntopo, toporesi)

    return h, hx, hy, x, y


def estimate_drag_coefficients(plat, plong, sfac):
    """
    Load the drag coefficient from datasets

    Inputs:
    -------
        - plat: point latitude
        - plong: point longitude
        - sfac: factor converting nautical to km

    Returns:
    --------
        - cdrag: drag coefficient
        - cdx, cdy: derivatives of cd in x & y
    """

    pifac = math.acos(-1) / 180  # pi number

    # Load neutral drag coefficients
    cd = tcr_io.load_netcdf_2d_parameters('../data', 'surface_data.nc', 'cdrag')

    # This corrects the drag coefficient to be better applied to gradient wind
    cd = 0.9 * cd / (1 + 50 * cd)
    # see Esau et al. (2004)
    # Align over-water values with Fortran (including some wave drag effect)
    cd = np.maximum(cd, 1e-3)

    # Interpolate drag coefficient and its gradients to POI
    sy = np.max(np.shape(plat))
    sx = np.max(np.shape(plong))

    # Gradients of drag coefficients to POI
    dcddx = sfac * (np.roll(cd, -1, 0) - cd)
    dcddy = sfac * (np.roll(cd, -1, 1) - cd)

    cdrag = np.zeros((sx, sy))
    cdx = np.zeros((sx, sy))
    cdy = np.zeros((sx, sy))

    for i in range(sx):
        for j in range(sy):
            # This is applied specifically for 0.25x0.25 resolution
            ib = np.floor(4 * plong[i]).astype(int)
            if ib >= 1440:
                ib -= 1440

            ibp = ib + 1
            if ibp > 1440 - 1:
                ibp = 0

            ibs = np.floor(4 * plong[i] - 0.5).astype(int)
            plongs = plong
            if ibs < 0:
                ibs += 1440
                plongs = plong + 360

            if ibs >= 1440 - 1:
                ibs -= 1440
                plongs = plong - 360

            ibsp = ibs + 1
            jb = np.floor(4 * (plat[j] + 90)).astype(int)
            jbs = np.floor(4 * (plat[j] + 90) - 0.5).astype(int)
            b1 = cd[ib, jb]
            b2 = cd[ib, jb + 1]
            b3 = cd[ibp, jb]
            b4 = cd[ibp, jb + 1]
            b1x = dcddx[ibs, jbs]
            b2x = dcddx[ibs, jbs + 1]
            b3x = dcddx[ibsp, jbs]
            b4x = dcddx[ibsp, jbs + 1]
            b1y = dcddy[ibs, jbs]
            b2y = dcddy[ibs, jbs + 1]
            b3y = dcddy[ibsp, jbs]
            b4y = dcddy[ibsp, jbs + 1]
            dely = 4 * (plat[j] + 90) - jb
            delx = 4 * plong[i] - ib
            d1 = (1 - delx) * (1 - dely)
            d2 = dely * (1.0 - delx)
            d3 = delx * (1.0 - dely)
            d4 = delx * dely
            cdrag[i, j] = d1 * b1 + d2 * b2 + d3 * b3 + d4 * b4
            dely = -0.5 + 4 * (plat[j] + 90) - jbs
            delx = -0.5 + 4 * plongs[i] - ibs
            d1 = (1.0 - delx) * (1.0 - dely)
            d2 = dely * (1 - delx)
            d3 = delx * (1 - dely)
            d4 = delx * dely
            cdx[i, j] = (d1*b1x + d2*b2x + d3*b3x + d4*b4x) / np.cos(pifac * plat[j])
            cdy[i, j] = d1*b1y + d2*b2y + d3*b3y + d4*b4y
    return cdrag, cdx, cdy


def calculate_qs900(T600, vmax):
    """
    Calculates saturation specific humidity at 600 hPa and saturation
    specific humidity at pressure 'pref' given 600 hPa T (K) and assuming a
    moist adiabatic lapse rate.

    Parameters:
    ----------
    - T600: Temperature at 600 hPa (K)
    - vmax: Maximum wind speed (knots)

    Returns:
    --------
    - q900: Saturation specific humidity at pref (950 hPa)
    - q600: Saturation specific humidity at 600 hPa
    """
    pref = 950  # Pressure to find qs at (hPa)

    # Constants
    cp = 1005
    Rv = 491
    Rd = 287
    Lv = 2.5e6

    c1 = Lv / Rv
    c2 = Rd * np.log(pref / 600)
    c3 = 1.6 / 100

    # Convert vmax from knots to m/s
    vmax = vmax * 1852 / 3600

    # Ensure Tc is not below -50°C
    Tc = np.clip(T600 - 273.15, -50, None)

    # Calculate es (saturation vapor pressure)
    es = 6.112 * np.exp(17.67 * Tc / (243.5 + Tc))

    # Calculate q600 (saturation specific humidity at 600 hPa)
    q600 = 0.622 * es / (600 - es)

    # Initialize q900
    q900 = np.zeros_like(q600)

    # First guess for T
    T = T600 + 20

    # Mask to ignore zero fill values of T600
    valid_mask = T600 > 100

    # Only compute for valid entries
    for _ in range(5):  # Iterative refinement (5 iterations)
        Tc_valid = T[valid_mask] - 273.15
        es_valid = 6.112 * np.exp(17.67 * Tc_valid / (243.5 + Tc_valid))
        qs_valid = 0.622 * es_valid / (pref - es_valid)
        er = (cp * np.log(T[valid_mask] / T600[valid_mask]) +
              Lv * (qs_valid / T[valid_mask] - q600[valid_mask] / T600[valid_mask]) -
              c2 - c3 * vmax[valid_mask] ** 2)
        derdT = (cp * T[valid_mask] + Lv * qs_valid * (c1 / T[valid_mask] - 1)) / T[valid_mask] ** 2
        T[valid_mask] -= er / derdT

    # Update q900 only for valid entries
    q900[valid_mask] = (
        0.622
        * 6.112
        * np.exp(17.67 * (T[valid_mask] - 273.15) / (243.5 + (T[valid_mask] - 273.15)))
        / (
            pref
            - 6.112
            * np.exp(
                17.67 * (T[valid_mask] - 273.15) / (243.5 + (T[valid_mask] - 273.15))
            )
        )
    )

    return q900, q600
