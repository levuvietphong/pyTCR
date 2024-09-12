"""
Input/Output functions for PyTCR
"""

import os
import glob
import numpy as np
import scipy.io as sio
import xarray as xr


def convert_to_mps(values, conversion_factor=0.514444):
    """Convert values to meters per second using the given conversion factor."""
    return values * conversion_factor


def filter_by_year(data, years, start_year, end_year):
    """Filter data by year range and return filtered data with indices."""
    indices = np.where((years >= start_year) & (years <= end_year))[0]
    return data[indices], indices


def load_Matlab_data(directory, filename):
    """
    Load data in Matlab format.

    Parameters:
    -----------
    directory : str
        Name of the directory where data is stored.
    filename : str
        Name of the Matlab file.

    Returns:
    --------
    dict
        Dictionary containing the loaded Matlab data.

    Raises:
    -------
    FileNotFoundError
        If the specified file does not exist.
    """
    try:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(current_dir)
        data_folder = os.path.join(parent_dir, directory)
        file_path = os.path.join(data_folder, filename)

        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The file {file_path} does not exist.")

        mat = sio.loadmat(file_path)
        return mat
    except Exception as e:
        raise RuntimeError(f"An error occurred while loading the Matlab file: {e}") from e


def load_netcdf_track_data(directory, filename):
    """
    Load data in NetCDF format from Tropical cyclone downscaling.

    Parameters:
    -----------
    directory : str
        Name of the directory where data is stored.
    filename : str
        Name of the NetCDF file.

    Returns:
    --------
    ds : xarray.Dataset
        The loaded NetCDF dataset.
    lat_trks : numpy.ndarray
        Array of latitude tracks.
    lon_trks : numpy.ndarray
        Array of longitude tracks.
    n_trk : numpy.ndarray
        Array of track numbers.
    v_trks : numpy.ndarray
        Array of track velocities.
    vmax_trks : numpy.ndarray
        Array of maximum track velocities.
    u850_trks : numpy.ndarray
        Array of 850 hPa zonal wind components.
    v850_trks : numpy.ndarray
        Array of 850 hPa meridional wind components.
    tc_month : numpy.ndarray
        Array of tropical cyclone months.
    tc_years : numpy.ndarray
        Array of tropical cyclone years.
    tc_time : numpy.ndarray
        Array of tropical cyclone times.
    """

    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_folder = os.path.join(current_dir, directory)
    file_path = os.path.join(data_folder, filename)

    ds = xr.open_dataset(file_path)
    lat_trks = np.nan_to_num(ds['lat_trks'].values)
    lon_trks = np.nan_to_num(ds['lon_trks'].values)
    n_trk = np.nan_to_num(ds['n_trk'].values)
    v_trks = np.nan_to_num(ds['v_trks'].values)
    vmax_trks = np.nan_to_num(ds['vmax_trks'].values)
    u850_trks = np.nan_to_num(ds['u850_trks'].values)
    v850_trks = np.nan_to_num(ds['v850_trks'].values)
    tc_month = np.nan_to_num(ds['tc_month'].values)
    tc_years = np.nan_to_num(ds['tc_years'].values)
    tc_time = np.nan_to_num(ds['time'].values)

    return (
        ds, lat_trks, lon_trks, n_trk, v_trks, vmax_trks, 
        u850_trks, v850_trks, tc_month, tc_years, tc_time
    )


def load_netcdf_2d_parameters(directory, filename, varname):
    """
    Load a specific variable from a NetCDF file.

    Parameters:
    -----------
    directory : str
        The directory where the NetCDF file is stored.
    filename : str
        The name of the NetCDF file.
    varname : str
        The name of the variable to be loaded from the NetCDF file.

    Returns:
    --------
    numpy.ndarray
        The data of the specified variable from the NetCDF file.
    """

    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_folder = os.path.join(current_dir, directory)
    file_path = os.path.join(data_folder, filename)

    ds = xr.open_dataset(file_path)
    data = ds[varname].values

    return data


def load_best_tracks_obs(fname, year_start, year_end):
    """
    Load the best tracks of observed tropical cyclones.

    Parameters:
    -----------
    fname : str
        The name of the observed tropical cyclone file.
    year_start : int
        The first year to load data.
    year_end : int
        The last year to load data.

    Returns:
    --------
    lat_tc : numpy.ndarray
        Latitudes of tropical cyclone tracks.
    lon_tc : numpy.ndarray
        Longitudes of tropical cyclone tracks.
    time_tc : numpy.ndarray
        Times of tropical cyclones.
    ind_tc : numpy.ndarray
        Indices of tropical cyclones.
    name_tc : numpy.ndarray
        Names of tropical cyclones.
    basin_tc : numpy.ndarray
        Names of ocean basins.
    wind_tc : numpy.ndarray
        Circular wind speeds of tropical cyclones (in m/s).
    speed_tc : numpy.ndarray
        Translational velocities of tropical cyclones (in m/s).
    """

    ds = xr.open_dataset(fname)
    name_tc = ds['name'].values
    lat_tc = ds['lat'].values
    lon_tc = ds['lon'].values
    time_tc = ds['time'].values
    yeartc = ds['season'].values
    wind_tc = convert_to_mps(ds['usa_wind'].values)
    speed_tc = convert_to_mps(ds['storm_speed'].values)
    basin_tc = ds['basin'].values

    # Filter TCs from year_start to year_end
    name_tc, ind = filter_by_year(name_tc, yeartc, year_start, year_end)
    lat_tc = lat_tc[ind]
    lon_tc = lon_tc[ind]
    wind_tc = wind_tc[ind]
    speed_tc = speed_tc[ind]
    time_tc = time_tc[ind]
    lon_tc[lon_tc < 0] += 360  # convert to 0-360 coordinate
    num_storms = name_tc.shape[0]  # get number of storms
    ind_tc = np.arange(num_storms)  # id from 0 to num_storms

    return lat_tc, lon_tc, time_tc, ind_tc, name_tc, basin_tc, wind_tc, speed_tc


def load_tracks_GCMs(pathdir, model='E3SM-1-0', basin='NA', expmnt='historical'):
    """
    Load the downscaled tracks of tropical cyclones from CMIP6 models.

    Parameters:
    -----------
    pathdir : str
        Path to the data directory.
    model : str
        Name of the CMIP6 model.
    basin : str
        Name of the ocean basin.
    expmnt : str
        Name of the model experiment.

    Returns:
    --------
    lat_trks : numpy.ndarray
        Latitudes of tropical cyclone tracks.
    lon_trks : numpy.ndarray
        Longitudes of tropical cyclone tracks.
    year_trks : numpy.ndarray
        Years of tropical cyclones.
    id_trks : numpy.ndarray
        Indices of tropical cyclones.
    vmax_trks : numpy.ndarray
        Maximum wind speeds of tropical cyclones.
    """
    try:
        ncfile = glob.glob(f"{pathdir}/{expmnt}/tracks_{basin}_{model}_*.nc")[0]
        ds = xr.open_dataset(ncfile)
        lon_trks = ds['lon_trks'].values
        lat_trks = ds['lat_trks'].values
        vmax_trks = ds['vmax_trks'].values
        year_trks = ds['year'].values
        id_trks = ds['n_trk'].values
    except IndexError as exc:
        raise FileNotFoundError(f"No files found for the given path: {pathdir}/{expmnt}/tracks_{basin}_{model}_*.nc") from exc
    except KeyError as e:
        raise KeyError(f"Missing expected data in the dataset: {e}") from e
    except Exception as e:
        raise RuntimeError(f"An error occurred while loading the dataset: {e}") from e

    return lat_trks, lon_trks, year_trks, id_trks, vmax_trks
