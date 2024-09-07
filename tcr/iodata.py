"""
Functions for I/O in PyTCR
"""

import os
import glob
import numpy as np
import scipy.io as sio
import xarray as xr
from netCDF4 import Dataset


def load_Matlab_data(directory, filename):
    """
    Load data in Matlab format

    Inputs:
    -------
        - directory: name of directory where data is stored
        - filename: name of the Matlab file

    Returns:
    --------
        - mat: tuple data
    """
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(current_dir)
    data_folder = os.path.join(parent_dir, directory)
    file_path = os.path.join(data_folder, filename)
    mat = sio.loadmat(file_path)
    return mat


def load_NetCDF_data(directory, filename):
    """
    Load data in NetCDF format from Tropical cyclone downscaling

    Inputs:
    -------
        - directory: name of directory where data is stored
        - filename: name of the net file

    Returns:
    --------
        - ds: tuple data
        - list variables
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

    return ds, lat_trks, lon_trks, n_trk, v_trks, vmax_trks, u850_trks, v850_trks, \
        tc_month, tc_years, tc_time


def load_best_tracks_obs(fname, year_start, year_end):
    """
    Load the best tracks of observed tropical cyclones

    Inputs:
    -------
        - fname : name of the observed tropical cyclone file
        - year_start : first year to load data
        - year_end : last year to load data

    Returns:
    --------
        - lat_tc : latitudes of TC tracks
        - lon_tc : longitudes of latitude of TC tracks
        - time_tc : time of TCs
        - ind_tc : index of TCs
        - name_tc : TC names
        - basins : names of ocean basins
        - windtc : circular wind speed of TCs
        - speed_tc : translational velocity of TCs
    """
    nc_file = Dataset(fname, 'r')
    variables = nc_file.variables
    nametc = variables['name'][:]
    lat_tc = variables['lat'][:]
    lon_tc = variables['lon'][:]
    time_tc = variables['time'][:]
    yeartc = variables['season'][:]
    usa_wind = variables['usa_wind'][:]
    wind_tc = usa_wind * 0.514444              # convert to unit m/s
    speed_tc = variables['storm_speed'][:]
    speed_tc = speed_tc * 0.514444        # convert to unit m/s
    basin = variables['basin'][:]

    # Filter TCs from year_start to year_end
    ind = np.where((yeartc >= year_start) & (yeartc <= year_end))[0]
    nametc = nametc[ind, :]
    lat_tc = lat_tc[ind, :]
    lon_tc = lon_tc[ind, :]
    wind_tc = wind_tc[ind, :]
    speed_tc = speed_tc[ind, :]
    time_tc = time_tc[ind]
    ind = lon_tc < 0
    lon_tc[ind] += 360                           # convert to 0-360 coordinate

    num_storms = nametc.shape[0]
    name_tc = []
    basins = []
    ind_tc = np.arange(0, num_storms, 1)

    # Decode storm and basin names from utf-8
    for k in range(num_storms):
        sd = ''.join(list(map(lambda x: x.decode('utf-8'), nametc.data[k, :]))).strip()
        name_tc.append(sd)
        bd = ''.join(list(map(lambda x: x.decode('utf-8'), basin.data[k, 1, :]))).strip()
        basins.append(bd)

    return lat_tc, lon_tc, time_tc, ind_tc, name_tc, basins, wind_tc, speed_tc


def load_tracks_GCMs(pathdir, modelid, basin, expmnt):
    """
    Load the downscaled tracks of tropical cyclones from CMIP6 models.

    Inputs:
    -------
        - pathdir : path to data directory
        - modelid : name of CMIP6 model
        - basin : name of ocean basin
        - expmnt : name of model experiments

    Returns:
    --------
        - lat_trks : latitudes of TC tracks
        - lon_trks : longitudes of latitude of TC tracks
        - year_trks : time of TCs
        - id_trks : index of TCs
        - vmax_trks : TC names
    """
    ncfile = glob.glob(f"{pathdir}/{expmnt}/tracks_{basin}_{modelid}_*.nc")[0]
    nc_fid = Dataset(ncfile, 'r')
    lon_trks = nc_fid.variables['lon_trks'][:]
    lat_trks = nc_fid.variables['lat_trks'][:]
    vmax_trks = nc_fid.variables['vmax_trks'][:]
    year_trks = nc_fid.variables['year'][:]
    id_trks = nc_fid.variables['n_trk'][:]

    return lat_trks, lon_trks, year_trks, id_trks, vmax_trks
