"""
Functions for topography and boundary layers in PyTCR
"""

import math
import numpy as np
from tcr import iodata as tcr_io


def calculate_distance_to_track(
    point_lat,
    point_lon,
    track_lats,
    track_lons,
    num_storms,
    storm_length,
    num_x_points,
    num_y_points,
    num_grid_points,
    degree_to_km_factor,
):
    """
    Calculate the distances of points of interest (POI) from the track.

    Parameters:
    -----------
    point_lat : numpy.ndarray
        Array of latitudes for points of interest.
    point_lon : numpy.ndarray
        Array of longitudes for points of interest.
    track_lats : numpy.ndarray
        Array of latitudes along each track.
    track_lons : numpy.ndarray
        Array of longitudes along each track.
    num_storms : int
        Number of storms.
    storm_length : int
        Length of each storm.
    num_x_points : int
        Number of points in the x direction.
    num_y_points : int
        Number of points in the y direction.
    num_grid_points : int
        Number of grid points.
    degree_to_km_factor : float
        Factor for converting degrees to kilometers.

    Returns:
    --------
    radius : numpy.ndarray
        Euclidean distance from track.
    dx : numpy.ndarray
        Distance components in the x direction.
    dy : numpy.ndarray
        Distance components in the y direction.
    """
    pi_factor = math.acos(-1) / 180
    dx = np.zeros((num_storms, storm_length, num_x_points, num_y_points))
    dy = np.zeros((num_storms, storm_length, num_x_points, num_y_points))

    for i in range(num_x_points):
        for jj in range(num_y_points):
            j = i if num_grid_points == 1 else jj

            if np.ndim(track_lons) == 1 and np.ndim(track_lats) == 1:
                dx[:, :, i, j] = (
                    degree_to_km_factor
                    * np.cos(pi_factor * point_lat[j])
                    * (point_lon[i] - track_lons)
                )
                dy[:, :, i, j] = degree_to_km_factor * (point_lat[j] - track_lats)
            elif np.ndim(track_lons) == 2 and np.ndim(track_lats) == 2:
                dx[:, :, i, j] = (
                    degree_to_km_factor
                    * np.cos(pi_factor * point_lat[j])
                    * (point_lon[i] - track_lons)
                )
                dy[:, :, i, j] = degree_to_km_factor * (point_lat[j] - track_lats)
            elif np.ndim(track_lons) == 4 and np.ndim(track_lats) == 4:
                dx[:, :, i, j] = (
                    degree_to_km_factor
                    * np.cos(pi_factor * point_lat[j])
                    * (point_lon[i] - track_lons[:, :, i, jj])
                )
                dy[:, :, i, j] = degree_to_km_factor * (
                    point_lat[j] - track_lats[:, :, i, jj]
                )

    radius = np.sqrt(dx**2 + dy**2)
    return radius, dx, dy


def calculate_spatial_derivatives(
    bathymetry,
    x_coords,
    y_coords,
    x_size,
    y_size,
    scale_factor,
    pi_factor,
    num_topo_rows,
    topo_resolution_inv,
):
    """
    Compute the spatial derivatives of a given topography.

    Parameters:
    ----------
    bathymetry : numpy.ndarray
        Array representing the bathymetry.
    x_coords : numpy.ndarray
        Array of spatial coordinates in the x direction.
    y_coords : numpy.ndarray
        Array of spatial coordinates in the y direction.
    x_size : int
        Size of the x dimension.
    y_size : int
        Size of the y dimension.
    scale_factor : float
        Factor for converting nautical miles to kilometers.
    pi_factor : float
        Value of pi divided by 180.
    num_topo_rows : int
        Number of rows in the bathymetry array.
    topo_resolution_inv : float
        Inverse of the topographic resolution.

    Returns:
    -------
    h : numpy.ndarray
        Array of topographic heights.
    hx : numpy.ndarray
        Array of derivatives of h in the x direction.
    hy : numpy.ndarray
        Array of derivatives of h in the y direction.
    """

    h = np.zeros((x_size, y_size))
    hx = np.zeros((x_size, y_size))
    hy = np.zeros((x_size, y_size))
    bathymetry = np.maximum(bathymetry, -1)
    dhdx = scale_factor * (np.roll(bathymetry, -1, 0) - np.roll(bathymetry, 0, 0))
    dhdy = scale_factor * (np.roll(bathymetry, -1, 1) - np.roll(bathymetry, 0, 1))

    for i in range(x_size):
        longitude = x_coords[i]
        if longitude >= 360:
            longitude -= 360

        for j in range(y_size):
            latitude = y_coords[j]
            ib = np.floor(topo_resolution_inv * longitude).astype(int)
            ibp = ib + 1
            if ibp >= num_topo_rows:
                ibp -= num_topo_rows

            ibs = np.floor(topo_resolution_inv * longitude - 0.5).astype(int)
            ibsp = ibs + 1
            if ibs < 0:
                ibs = num_topo_rows - 1
                longitude += 360

            if ibsp >= num_topo_rows:
                ibsp -= num_topo_rows

            jb = np.floor(topo_resolution_inv * (latitude + 90)).astype(int)
            jbs = np.floor(topo_resolution_inv * (latitude + 90) - 0.5).astype(int)
            b1, b2, b3, b4 = (
                bathymetry[ib, jb],
                bathymetry[ib, jb + 1],
                bathymetry[ibp, jb],
                bathymetry[ibp, jb + 1],
            )
            dely, delx = (
                topo_resolution_inv * (latitude + 90) - jb,
                topo_resolution_inv * longitude - ib,
            )
            d1, d2, d3, d4 = (
                (1 - delx) * (1 - dely),
                dely * (1 - delx),
                delx * (1 - dely),
                delx * dely,
            )
            h[i, j] = (
                np.exp(
                    d1 * np.log(b1 + 11)
                    + d2 * np.log(b2 + 11)
                    + d3 * np.log(b3 + 11)
                    + d4 * np.log(b4 + 11)
                )
                - 11
            )

            b1, b2, b3, b4 = (
                dhdx[ibs, jbs],
                dhdx[ibs, jbs + 1],
                dhdx[ibsp, jbs],
                dhdx[ibsp, jbs + 1],
            )
            dely, delx = (
                -0.5 + topo_resolution_inv * (latitude + 90) - jbs,
                -0.5 + topo_resolution_inv * longitude - ibs,
            )
            d1, d2, d3, d4 = (
                (1 - delx) * (1 - dely),
                dely * (1 - delx),
                delx * (1 - dely),
                delx * dely,
            )
            hx[i, j] = (b1 * d1 + b2 * d2 + b3 * d3 + b4 * d4) / np.cos(
                pi_factor * latitude
            )

            b1, b2, b3, b4 = (
                dhdy[ibs, jbs],
                dhdy[ibs, jbs + 1],
                dhdy[ibsp, jbs],
                dhdy[ibsp, jbs + 1],
            )
            hy[i, j] = b1 * d1 + b2 * d2 + b3 * d3 + b4 * d4

    return h, hx, hy


def estimate_topographic_height(bxmin, bxmax, bymin, bymax, dellatlong):
    """
    Estimate the topographic height and its spatial derivatives.

    Parameters:
    -----------
    bxmin : float
        Minimum bound in the x-axis (longitude).
    bxmax : float
        Maximum bound in the x-axis (longitude).
    bymin : float
        Minimum bound in the y-axis (latitude).
    bymax : float
        Maximum bound in the y-axis (latitude).
    dellatlong : float
        Horizontal resolution of the map in degrees.

    Returns:
    --------
    h : numpy.ndarray
        Topographic height data.
    hx : numpy.ndarray
        Topographic gradient in the x direction.
    hy : numpy.ndarray
        Topographic gradient in the y direction.
    x : numpy.ndarray
        Array of longitude values for the grid.
    y : numpy.ndarray
        Array of latitude values for the grid.
    """

    # Load high-resolution bathymetry data
    bathymetry = tcr_io.load_netcdf_2d_parameters(
        "../data", "surface_data.nc", "bathymetry_high"
    )
    ntopo, _ = np.shape(bathymetry)  # number of grid point in x direction
    topo_resolution = 360 / ntopo
    topo_resolution_inv = 1 / topo_resolution
    scale_factor = 1.0 / (
        topo_resolution * 60.0 * 1852
    )  # factor converting degree to m
    pi_factor = math.pi / 180  # pi number

    x = np.round(np.arange(bxmin, bxmax + 1e-8, dellatlong), 4)
    y = np.round(np.arange(bymin, bymax + 1e-8, dellatlong), 4)
    x_size = len(x)
    y_size = len(y)

    if x_size == y_size:
        x_size += 1
        x = np.append(x, bxmax + dellatlong)

    h, hx, hy = calculate_spatial_derivatives(
        bathymetry,
        x,
        y,
        x_size,
        y_size,
        scale_factor,
        pi_factor,
        ntopo,
        topo_resolution_inv,
    )

    return h, hx, hy, x, y


def estimate_drag_coefficients(plat, plong, sfac):
    """
    Estimate the drag coefficients and their gradients at points of interest (POI).

    Parameters:
    -----------
    plat : numpy.ndarray
        Array of latitudes for points of interest.
    plong : numpy.ndarray
        Array of longitudes for points of interest.
    sfac : float
        Factor for converting nautical miles to kilometers.

    Returns:
    --------
    cdrag : numpy.ndarray
        Array of drag coefficients at POI.
    cdx : numpy.ndarray
        Array of drag coefficient gradients in the x direction.
    cdy : numpy.ndarray
        Array of drag coefficient gradients in the y direction.
    """

    pi_factor = math.pi / 180  # pi number

    # Load neutral drag coefficients
    cd = tcr_io.load_netcdf_2d_parameters("../data", "surface_data.nc", "cdrag")

    # Correct the drag coefficient for gradient wind application
    cd = 0.9 * cd / (1 + 50 * cd)
    cd = np.maximum(cd, 1e-3)  # Ensure minimum drag coefficient value

    # Calculate gradients of drag coefficients
    dcddx = sfac * (np.roll(cd, -1, axis=0) - cd)
    dcddy = sfac * (np.roll(cd, -1, axis=1) - cd)

    # Initialize arrays for drag coefficients and their gradients at POI
    sx, sy = plong.size, plat.size
    cdrag = np.zeros((sx, sy))
    cdx = np.zeros((sx, sy))
    cdy = np.zeros((sx, sy))

    for i in range(sx):
        for j in range(sy):
            # Calculate indices for interpolation
            ib = int(np.floor(4 * plong[i])) % 1440
            ibp = (ib + 1) % 1440
            ibs = int(np.floor(4 * plong[i] - 0.5)) % 1440
            ibsp = (ibs + 1) % 1440
            jb = int(np.floor(4 * (plat[j] + 90)))
            jbs = int(np.floor(4 * (plat[j] + 90) - 0.5))

            # Interpolate drag coefficient and its gradients
            delx = 4 * plong[i] - ib
            dely = 4 * (plat[j] + 90) - jb
            weights = [
                (1 - delx) * (1 - dely),
                dely * (1 - delx),
                delx * (1 - dely),
                delx * dely,
            ]

            cdrag[i, j] = (
                weights[0] * cd[ib, jb]
                + weights[1] * cd[ib, jb + 1]
                + weights[2] * cd[ibp, jb]
                + weights[3] * cd[ibp, jb + 1]
            )

            delx = -0.5 + 4 * plong[i] - ibs
            dely = -0.5 + 4 * (plat[j] + 90) - jbs
            weights = [
                (1 - delx) * (1 - dely),
                dely * (1 - delx),
                delx * (1 - dely),
                delx * dely,
            ]

            cdx[i, j] = (
                weights[0] * dcddx[ibs, jbs]
                + weights[1] * dcddx[ibs, jbs + 1]
                + weights[2] * dcddx[ibsp, jbs]
                + weights[3] * dcddx[ibsp, jbs + 1]
            ) / np.cos(pi_factor * plat[j])

            cdy[i, j] = (
                weights[0] * dcddy[ibs, jbs]
                + weights[1] * dcddy[ibs, jbs + 1]
                + weights[2] * dcddy[ibsp, jbs]
                + weights[3] * dcddy[ibsp, jbs + 1]
            )

    return cdrag, cdx, cdy


def calculate_qs900(T600, vmax):
    """
    Calculate saturation specific humidity at 600 hPa and 950 hPa given the
    temperature at 600 hPa and assuming a moist adiabatic lapse rate.

    Parameters:
    ----------
    T600 : numpy.ndarray
        Temperature at 600 hPa (K).
    vmax : numpy.ndarray
        Maximum wind speed (knots).

    Returns:
    --------
    q900 : numpy.ndarray
        Saturation specific humidity at 950 hPa.
    q600 : numpy.ndarray
        Saturation specific humidity at 600 hPa.
    """
    pref = 950  # Pressure to find qs at (hPa)

    # Constants
    cp = 1005  # Specific heat capacity of dry air at constant pressure (J/(kg·K))
    Rv = 491  # Gas constant for water vapor (J/(kg·K))
    Rd = 287  # Gas constant for dry air (J/(kg·K))
    Lv = 2.5e6  # Latent heat of vaporization (J/kg)

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
        er = (
            cp * np.log(T[valid_mask] / T600[valid_mask])
            + Lv * (qs_valid / T[valid_mask] - q600[valid_mask] / T600[valid_mask])
            - c2
            - c3 * vmax[valid_mask] ** 2
        )
        derdT = (cp * T[valid_mask] + Lv * qs_valid * (c1 / T[valid_mask] - 1)) / T[
            valid_mask
        ] ** 2
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
