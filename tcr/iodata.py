"""
Input/Output functions for PyTCR
"""

import os
import glob
import numpy as np
import scipy.io as sio
import xarray as xr
import requests
import shapefile
from bs4 import BeautifulSoup
from tcr import datadir as tcr_data
BASE_DATA_DIR = tcr_data.BASE_DATA_DIR
DOWNSCALED_DATA_DIR = tcr_data.get_downscaled_data_dir()


def convert_to_mps(values, conversion_factor=0.514444):
    """
    Convert the provided values to meters per second (m/s) using the specified conversion factor.

    Parameters:
    -----------
    values : array-like
        The input values to be converted.
    conversion_factor : float, optional
        The factor used for conversion to meters per second (default is 0.514444).

    Returns:
    --------
    array-like
        The converted values in meters per second (m/s).
    """
    return values * conversion_factor


def filter_by_year(data, years, start_year, end_year):
    """
    Filter the input data based on a specified year range and return the filtered data along with their indices.

    Parameters:
    -----------
    data : array-like
        The input data to be filtered.
    years : array-like
        Array of years corresponding to the input data.
    start_year : int
        The starting year of the desired range (inclusive).
    end_year : int
        The ending year of the desired range (inclusive).

    Returns:
    --------
    tuple
        A tuple containing the filtered data and the indices of the filtered elements.
    """
    indices = np.where((years >= start_year) & (years <= end_year))[0]
    return data[indices], indices


def load_Matlab_data(directory, filename):
    """
    Load data from a Matlab (.mat) file.

    Parameters:
    -----------
    directory : str
        The directory containing the Matlab file.
    filename : str
        The name of the Matlab file to be loaded.

    Returns:
    --------
    dict
        A dictionary containing the data loaded from the Matlab file.

    Raises:
    -------
    FileNotFoundError
        If the specified file does not exist in the given directory.
    RuntimeError
        If an error occurs while loading the Matlab file.
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


def load_tracks_GCMs(
    data_directory=None, model="E3SM-1-0", basin="NA", expmnt="historical", dropna=True,
):
    """
    Load and extract data from a NetCDF file containing tropical cyclone downscaling information.

    Parameters:
    -----------
    data_directory : str, optional
        Path to the directory containing the NetCDF files.
        If not provided, defaults to the downscaled data directory.
    model : str, optional
        CMIP6 model name (e.g., 'E3SM-1-0'). Default is 'E3SM-1-0'.
    basin : str, optional
        Name of the ocean basin (e.g., 'NA' for North Atlantic, 'WP' for Western Pacific).
        Default is 'NA'.
    expmnt : str, optional
        Experiment name (e.g., 'historical', 'ssp585'). Default is 'historical'.
    dropna : bool, optional
        Whether to convert NaNs to zeros. Default is True.        

    Returns:
    --------
    ds : xarray.Dataset
        The loaded NetCDF dataset containing tropical cyclone data.
    lat_trks : numpy.ndarray
        Latitude coordinates of tropical cyclone tracks (in degrees).
    lon_trks : numpy.ndarray
        Longitude coordinates of tropical cyclone tracks (in degrees).
    n_trk : numpy.ndarray
        Track numbers identifying individual tropical cyclones.
    v_trks : numpy.ndarray
        Azimuthal wind speeds along the tracks (in m/s).
    vmax_trks : numpy.ndarray
        Maximum wind speeds or intensities along the tracks (in m/s).
    u850_trks : numpy.ndarray
        Zonal wind components at 850 hPa along the tracks (in m/s).
    v850_trks : numpy.ndarray
        Meridional wind components at 850 hPa along the tracks (in m/s).
    tc_month : numpy.ndarray
        Months corresponding to the tropical cyclone occurrences.
    tc_years : numpy.ndarray
        Years corresponding to the tropical cyclone occurrences.
    tc_time : numpy.ndarray
        Timestamps of the tropical cyclone occurrences.

    Raises:
    -------
    FileNotFoundError
        If no matching NetCDF file is found in the specified directory.
    KeyError
        If the expected variables are missing in the NetCDF dataset.
    RuntimeError
        If an error occurs while loading or processing the dataset.
    """
    if data_directory is None:
        data_directory = os.path.join(DOWNSCALED_DATA_DIR, 'downscaled')

    ncfilename = os.path.join(data_directory, expmnt, f"tracks_{basin}_{model}_*.nc")
    try:
        ncfile = glob.glob(ncfilename)[0]
        with xr.open_dataset(ncfile) as file_handle:
            ds = file_handle.load()

        def get_vals(var):
            return np.nan_to_num(ds[var].values) if dropna else ds[var].values

        var_names = ['lat_trks', 'lon_trks', 'n_trk', 'v_trks', 'vmax_trks', 'year',
                     'u850_trks', 'v850_trks', 'tc_month', 'tc_years', 'time']

        lat_trks, lon_trks, n_trk, v_trks, vmax_trks, year_trks, \
            u850_trks, v850_trks, tc_month, tc_years, tc_time = map(get_vals, var_names)

    except IndexError as e:
        raise FileNotFoundError(f"No files found for the given path: {ncfilename}") from e
    except KeyError as e:
        raise KeyError(f"Missing expected data in the dataset: {e}") from e
    except Exception as e:
        raise RuntimeError(f"An error occurred while loading the dataset: {e}") from e

    return (
        ds, lat_trks, lon_trks, year_trks, n_trk, vmax_trks, v_trks,
        u850_trks, v850_trks, tc_month, tc_years, tc_time
    )


def load_netcdf_2d_parameters(directory, filename, varname):
    """
    Load a specific variable from a NetCDF file.

    Parameters:
    -----------
    directory : str
        Path to the directory containing the NetCDF file.
    filename : str
        Name of the NetCDF file to be loaded.
    varname : str
        Name of the variable to extract from the NetCDF file.

    Returns:
    --------
    numpy.ndarray
        Array containing the data of the specified variable.

    Raises:
    -------
    FileNotFoundError
        If the specified NetCDF file does not exist in the given directory.
    KeyError
        If the specified variable is not found in the NetCDF file.
    RuntimeError
        If an error occurs while loading the NetCDF file.
    """

    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_folder = os.path.join(current_dir, directory)
    file_path = os.path.join(data_folder, filename)

    ds = xr.open_dataset(file_path)
    data = ds[varname].values

    return data


def load_best_tracks_obs(fname, year_start, year_end):
    """
    Load the best tracks of observed tropical cyclones (IBTrACS).

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
        Latitudes of tropical cyclone tracks (degree).
    lon_tc : numpy.ndarray
        Longitudes of tropical cyclone tracks (degree).
    time_tc : numpy.ndarray
        Times of tropical cyclones (YYYY-MM-DD HH:MM:SS)
    ind_tc : numpy.ndarray
        Indices of tropical cyclones.
    name_tc : numpy.ndarray
        Names of tropical cyclones.
    basin_tc : numpy.ndarray
        Names of ocean basins.
    wind_tc : numpy.ndarray
        Circular wind speeds of tropical cyclones (m/s).
    speed_tc : numpy.ndarray
        Translational velocities of tropical cyclones (m/s).
    """

    ds = xr.open_dataset(fname)
    name_tc = ds['name'].values
    lat_tc = ds['lat'].values
    lon_tc = ds['lon'].values
    time_tc = ds['time'].values
    yeartc = ds['season'].values
    # unit in
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


def fetch_directory_tree(url=None, depth=0, max_depth=2):
    """
    Fetch and print the directory tree from a specified URL up to a given depth.

    Parameters:
    -----------
    url : str, optional
        URL to the directory containing the downscaled tracks (default is TACC url).
    depth : int, optional
        Current depth of the directory tree (default is 0).
    max_depth : int, optional
        Maximum depth to traverse the directory tree (default is 2).

    Returns:
    --------
    None
    """
    if url is None:
        url = "https://web.corral.tacc.utexas.edu/setxuifl/tropical_cyclones/downscaled_cmip6_tracks"

    # Make a GET request to fetch the HTML content
    response = requests.get(url, timeout=10)
    if response.status_code != 200:
        print(f"Failed to retrieve contents of {url}")
        return

    # Parse the HTML content using BeautifulSoup
    soup = BeautifulSoup(response.text, 'html.parser')

    # Find all links in the directory page
    links = soup.find_all('a')

    # Iterate through the links
    for link in links:
        href = link.get('href')

        # Ignore parent directory links and base URL itself
        if href == "../" or href == "/":
            continue

        # Define colors for different depths
        colors = ["\033[1m\033[94m", "\033[93m", "\033[92m", "\033[91m", "\033[95m", "\033[96m"]
        color = colors[depth % len(colors)]

        # Indent for subdirectory depth
        print(f"{color}{'    ' * depth + '|-- ' + link.text}\033[0m")

        # If it's a directory (ends with '/') and within max_depth, fetch its contents recursively
        if href.endswith('/') and depth < max_depth - 1:
            fetch_directory_tree(os.path.join(url, href), depth + 1, max_depth)


def download_tracks_data_cmip6(
    url=None, experiments='historical', models='E3SM-1-0', target_directory=None
):
    """
    Download the downscaled tropical cyclone tracks from CMIP6 models.

    Parameters:
    -----------
    url : str, optional
        URL to the directory containing the downscaled tracks (default is TACC url).
    experiment : str, optional
        Name of the model experiment (default is 'historical').
    model : str, optional
        Name of the CMIP6 model (default is 'E3SM-1-0').
    target_directory : str, optional
        Path to the target directory to save the downloaded files.
        Default is '{BASE_DATA_DIR}/downscaled/'

    Returns:
    --------
    None
    """
    if isinstance(experiments, str):
        experiments = [experiments]

    if isinstance(models, str):
        models = [models]

    if url is None:
        url = "https://web.corral.tacc.utexas.edu/setxuifl/tropical_cyclones/downscaled_cmip6_tracks"

    if target_directory is None:
        target_directory = os.path.join(DOWNSCALED_DATA_DIR, 'downscaled')

    for experiment in experiments:
        for model in models:
            # Create a local directory to save the downloaded files
            os.makedirs(os.path.join(target_directory, experiment), exist_ok=True)

            # Get the list of netcdf files in the folder
            response = requests.get(f'{url}/{experiment}/{model}/', timeout=10)
            if response.status_code != 200:
                print(
                    f"Error: Unable to access the directory for model '{model}' in experiment '{experiment}'. Please check the model name and try again."
                )
                continue
            soup = BeautifulSoup(response.text, 'html.parser')

            # Extract the links to the netcdf files
            links = soup.find_all('a')

            for link in links:
                href = link.get('href')
                if href.endswith('.nc'):
                    file_path = os.path.join(target_directory, experiment, href)
                    if os.path.exists(file_path):
                        print(f"File {href} already exists. Skipping download.")
                        continue
                    print(f"Downloading {href}...")
                    # Download the file
                    file_url = f'{url}/{experiment}/{model}/{href}'
                    response = requests.get(file_url, stream=True, timeout=10)
                    total_size = int(response.headers.get('content-length', 0))
                    block_size = 1024  # 1 Kibibyte
                    downloaded_size = 0
                    # Save the file to the target_directory/experiment folder
                    with open(file_path, 'wb') as f:
                        for data in response.iter_content(block_size):
                            downloaded_size += len(data)
                            f.write(data)
                            mb_downloaded = downloaded_size / (1024 * 1024)
                            mb_total = total_size / (1024 * 1024)
                            percentage = (downloaded_size / total_size) * 100
                            print(f"{mb_downloaded:.2f} MB of {mb_total:.2f} MB ({percentage:.2f}%)", end='\r')
                    print()  # Move to the next line after download is complete


def get_bbox_from_shapefile(shapefile_path):
    """
    Extract the bounding box coordinates from a shapefile.

    Parameters:
    -----------
    shapefile_path : str
        The file path to the shapefile.

    Returns:
    --------
    tuple
        A tuple representing the bounding box in the format (xmin, xmax, ymin, ymax).

    Raises:
    -------
    FileNotFoundError
        If the specified shapefile does not exist.
    RuntimeError
        If an error occurs while reading the shapefile.
    """
    if not os.path.exists(shapefile_path):
        raise FileNotFoundError(f"Shapefile not found: {shapefile_path}")

    try:
        with shapefile.Reader(shapefile_path) as sf:
            bbox = sf.bbox  # [xmin, ymin, xmax, ymax]
            return bbox[0], bbox[2], bbox[1], bbox[3]
    except Exception as e:
        raise RuntimeError(f"Failed to read shapefile: {e}") from e
