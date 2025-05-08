"""
Functions for topography and boundary layers in PyTCR
"""

import math
import numpy as np
from tcr import iodata as tcr_io
from tcr import parameters as tcr_params
from tcr.datadir import BASE_DATA_DIR


def calculate_distance_to_track(
    point_lat, point_lon, track_lats, track_lons, num_storms, storm_length,
    num_x_points, num_y_points, num_grid_points, degree_to_km
):
    """
    Calculate the distances of points of interest (POI) from the track.

    Parameters:
    -----------
    point_lat : numpy.ndarray
        Array of latitudes for points of interest (degree)
    point_lon : numpy.ndarray
        Array of longitudes for points of interest (degree)
    track_lats : numpy.ndarray
        Array of latitudes along each track (degree)
    track_lons : numpy.ndarray
        Array of longitudes along each track (degree)
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
    degree_to_km : float
        Factor for converting degrees to kilometers.

    Returns:
    --------
    radius : numpy.ndarray
        Euclidean distance from track (km).
    dx : numpy.ndarray
        Distance components in the x direction (km).
    dy : numpy.ndarray
        Distance components in the y direction (km).
    """
    pi_factor = tcr_params.RAD2DEG  # Use numpy's pi constant
    shape = (num_storms, storm_length, num_x_points, num_y_points)
    dx = np.zeros(shape)
    dy = np.zeros(shape)

    for i in range(num_x_points):
        for jj in range(num_y_points):
            j = i if num_grid_points == 1 else jj

            if np.ndim(track_lons) == 1 and np.ndim(track_lats) == 1:
                dx[:, :, i, j] = np.cos(pi_factor * point_lat[j]) * (point_lon[i] - track_lons)
                dy[:, :, i, j] = point_lat[j] - track_lats
            elif np.ndim(track_lons) == 2 and np.ndim(track_lats) == 2:
                dx[:, :, i, j] = np.cos(pi_factor * point_lat[j]) * (point_lon[i] - track_lons)
                dy[:, :, i, j] = point_lat[j] - track_lats
            elif np.ndim(track_lons) == 4 and np.ndim(track_lats) == 4:
                dx[:, :, i, j] = (np.cos(pi_factor * point_lat[j]) *
                                  (point_lon[i] - track_lons[:, :, i, jj]))
                dy[:, :, i, j] = point_lat[j] - track_lats[:, :, i, jj]

    dx *= degree_to_km
    dy *= degree_to_km
    radius = np.sqrt(dx**2 + dy**2)
    return radius, dx, dy


def calculate_spatial_derivatives(
    bathymetry, x_coords, y_coords, x_size, y_size, scale_factor,
    pi_factor, num_topo_rows, topo_resolution_inv
):
    """
    Compute the spatial derivatives of a given topography.

    Parameters:
    -----------
    bathymetry : numpy.ndarray
        Array representing the bathymetry (m)
    x_coords : numpy.ndarray
        Array of spatial coordinates in the x direction (degree)
    y_coords : numpy.ndarray
        Array of spatial coordinates in the y direction (degree)
    x_size : int
        Size of the x dimension.
    y_size : int
        Size of the y dimension.
    scale_factor : float
        Factor for converting nautical miles to meters.
    pi_factor : float
        Value of pi divided by 180.
    num_topo_rows : int
        Number of rows in the bathymetry array.
    topo_resolution_inv : float
        Inverse of the topographic resolution.

    Returns:
    --------
    h : numpy.ndarray
        Array of topographic heights (m).
    hx : numpy.ndarray
        Array of derivatives of h in the x direction.
    hy : numpy.ndarray
        Array of derivatives of h in the y direction.
    """

    def calculate_indices(lon, lat):
        """ Helper function for calculating indices """
        ib = int(topo_resolution_inv * lon)
        ibp = (ib + 1) % num_topo_rows
        ibs = int(topo_resolution_inv * lon - 0.5)
        ibsp = (ibs + 1) % num_topo_rows
        jb = int(topo_resolution_inv * (lat + 90))
        jbs = int(topo_resolution_inv * (lat + 90) - 0.5)
        return ib, ibp, ibs, ibsp, jb, jbs

    def weighted_sum(data, indices, delx, dely, logfac=True):
        """ Helper function for weighted sum """
        ib, ibp, _, _, jb, _ = indices
        b = [data[ib, jb], data[ib, jb + 1], data[ibp, jb], data[ibp, jb + 1]]
        d = [(1 - delx) * (1 - dely), dely * (1 - delx), delx * (1 - dely), delx * dely]
        if logfac:
            return sum(np.log(b[i]+11) * d[i] for i in range(4))
        else:
            return sum(b[i] * d[i] for i in range(4))

    # since bathymetry is in meters, need to convert space resolution from degree to meters
    # for topographic derivatives
    h, hx, hy = [np.zeros((x_size, y_size)) for _ in range(3)]
    bathymetry = np.maximum(np.float32(bathymetry), -1)
    dhdx = scale_factor * (np.roll(bathymetry, -1, 0) - np.roll(bathymetry, 0, 0))
    dhdy = scale_factor * (np.roll(bathymetry, -1, 1) - np.roll(bathymetry, 0, 1))

    for i, lon in enumerate(x_coords):
        lon = lon if lon < 360 else lon - 360

        for j, lat in enumerate(y_coords):
            indices = calculate_indices(lon, lat)
            ib, _, ibs, ibsp, jb, jbs = indices
            delx, dely = topo_resolution_inv * lon - ib, topo_resolution_inv * (lat + 90) - jb
            h[i, j] = np.exp(weighted_sum(bathymetry, indices, delx, dely, logfac=True)) - 11

            delx = -0.5 + topo_resolution_inv * lon - ibs
            dely = -0.5 + topo_resolution_inv * (lat + 90) - jbs
            hx[i, j] = weighted_sum(dhdx, (ibs, ibsp, None, None, jbs, None),
                                    delx, dely, logfac=False) / np.cos(pi_factor * lat)
            hy[i, j] = weighted_sum(dhdy, (ibs, ibsp, None, None, jbs, None),
                                    delx, dely, logfac=False)
    return h, hx, hy


def estimate_topographic_height(bxmin, bxmax, bymin, bymax, dellatlong):
    """
    Estimate the topographic height and its spatial derivatives.

    Parameters:
    -----------
    bxmin : float
        Minimum bound in the x-axis (longitude in degree)
    bxmax : float
        Maximum bound in the x-axis (longitude in degree)
    bymin : float
        Minimum bound in the y-axis (latitude in degree)
    bymax : float
        Maximum bound in the y-axis (latitude in degree)
    dellatlong : float
        Horizontal resolution of the map in degrees.

    Returns:
    --------
    h : numpy.ndarray
        Topographic height data (m)
    hx : numpy.ndarray
        Topographic gradient in the x direction
    hy : numpy.ndarray
        Topographic gradient in the y direction
    x : numpy.ndarray
        Array of longitude values for the grid (degree)
    y : numpy.ndarray
        Array of latitude values for the grid (degree)
    """

    # Load high-resolution bathymetry data
    bathymetry = tcr_io.load_netcdf_2d_parameters(
        BASE_DATA_DIR, "surface_data.nc", "bathymetry_high"
    )
    ntopo, _ = np.shape(bathymetry)  # number of grid point in x direction
    topo_resolution = 360 / ntopo
    topo_resolution_inv = 1 / topo_resolution
    pi_factor = tcr_params.RAD2DEG
    # scale factor converting degree to meter
    scale_factor = 1.0 / (tcr_params.DEG2KM * topo_resolution * 1000)

    x = np.round(np.arange(bxmin, bxmax + 1e-8, dellatlong), 4)
    y = np.round(np.arange(bymin, bymax + 1e-8, dellatlong), 4)
    x_size = len(x)
    y_size = len(y)

    if x_size == y_size:
        x_size += 1
        x = np.append(x, bxmax + dellatlong)

    h, hx, hy = calculate_spatial_derivatives(
        bathymetry, x, y, x_size, y_size, scale_factor,
        pi_factor, ntopo, topo_resolution_inv,
    )

    return h, hx, hy, x, y


def estimate_drag_coefficients(plat, plong, sfac):
    """
    Estimate the drag coefficients and their gradients at points of interest (POI).

    Parameters:
    -----------
    plat : numpy.ndarray
        Array of latitudes for points of interest (degree)
    plong : numpy.ndarray
        Array of longitudes for points of interest (degree)
    sfac : float
        Factor for converting spatial resolution in degree to meters

    Returns:
    --------
    cdrag : numpy.ndarray
        Array of drag coefficients at POI (-)
    cdx : numpy.ndarray
        Array of drag coefficient gradients in the x direction (1/m)
    cdy : numpy.ndarray
        Array of drag coefficient gradients in the y direction (1/m)
    """

    def _interpolate(data, ix, iy, dx, dy):
        """Helper function for bilinear interpolation"""
        weights = [(1 - dx) * (1 - dy), dy * (1 - dx), dx * (1 - dy), dx * dy]
        return (weights[0] * data[ix, iy] +
                weights[1] * data[ix, iy + 1] +
                weights[2] * data[(ix + 1) % 1440, iy] +
                weights[3] * data[(ix + 1) % 1440, iy + 1])

    pi_factor = tcr_params.RAD2DEG

    # Load neutral drag coefficients
    cd = tcr_io.load_netcdf_2d_parameters(BASE_DATA_DIR, "surface_data.nc", "cdrag")

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
            ibs = int(np.floor(4 * plong[i] - 0.5)) % 1440
            jb = int(np.floor(4 * (plat[j] + 90)))
            jbs = int(np.floor(4 * (plat[j] + 90) - 0.5))

            # Interpolate drag coefficient and its gradients
            delx = 4 * plong[i] - ib
            dely = 4 * (plat[j] + 90) - jb
            cdrag[i, j] = _interpolate(cd, ib, jb, delx, dely)

            delx = -0.5 + 4 * plong[i] - ibs
            dely = -0.5 + 4 * (plat[j] + 90) - jbs
            cdx[i, j] = _interpolate(dcddx, ibs, jbs, delx, dely) / np.cos(pi_factor * plat[j])
            cdy[i, j] = _interpolate(dcddy, ibs, jbs, delx, dely)
    return cdrag, cdx, cdy


def calculate_qs900(T600, vmax):
    """
    Calculate saturation specific humidity at 600 hPa and 950 hPa given the
    temperature at 600 hPa and assuming a moist adiabatic lapse rate.

    Parameters:
    -----------
    T600 : numpy.ndarray
        Temperature at 600 hPa (K).
    vmax : numpy.ndarray
        Maximum wind speed (knots).

    Returns:
    --------
    q900 : numpy.ndarray
        Saturation specific humidity at 950 hPa (g/g)
    q600 : numpy.ndarray
        Saturation specific humidity at 600 hPa (g/g)
    """
    # Constants
    CONSTANTS = {
        'pref': 950,  # Pressure to find qs at (hPa)
        'cp': 1005,  # Specific heat capacity of dry air at constant pressure (J/(kg·K))
        'Rv': 491,  # Gas constant for water vapor (J/(kg·K))
        'Rd': 287,  # Gas constant for dry air (J/(kg·K))
        'Lv': 2.5e6,  # Latent heat of vaporization (J/kg)
    }

    def calc_es_qs(T, p):
        Tc = np.clip(T - 273.15, -50, None)
        es = 6.112 * np.exp(17.67 * Tc / (243.5 + Tc))
        q600 = 0.622 * es / (p - es)
        return es, q600

    c1 = CONSTANTS['Lv'] / CONSTANTS['Rv']
    c2 = CONSTANTS['Rd'] * np.log(CONSTANTS['pref'] / 600)
    c3 = 1.6 / 100

    # Convert vmax from knots to m/s
    vmax = vmax * 1852 / 3600
    _, q600 = calc_es_qs(T600, 600)

    # First guess for T
    T = T600 + 20

    # Mask to ignore zero fill values of T600
    valid_mask = T600 > 100

    # Only compute for valid entries
    for _ in range(5):  # Iterative refinement (5 iterations)
        _, qs = calc_es_qs(T[valid_mask], CONSTANTS['pref'])
        er = (CONSTANTS['cp'] * np.log(T[valid_mask] / T600[valid_mask])
              + CONSTANTS['Lv'] * (qs / T[valid_mask] - q600[valid_mask] / T600[valid_mask])
              - c2 - c3 * vmax[valid_mask] ** 2)
        derdT = (CONSTANTS['cp'] * T[valid_mask] +
                 CONSTANTS['Lv'] * qs * (c1 / T[valid_mask] - 1)) / T[valid_mask] ** 2
        T[valid_mask] -= er / derdT

    q900 = np.zeros_like(q600)
    _, q900[valid_mask] = calc_es_qs(T[valid_mask], CONSTANTS['pref'])

    return q900, q600
