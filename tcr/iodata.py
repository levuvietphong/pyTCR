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
from tcr.datadir import DATA_DIR

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


def load_netcdf_track_data(data_directory=os.path.join(DATA_DIR, 'downscaled'),
                           model='E3SM-1-0', basin='NA', expmnt='historical'):
    """
    Load data in NetCDF format from Tropical cyclone downscaling.

    Parameters:
    -----------
    data_directory : str
        Path to the data directory.
    model : str
        Name of the CMIP6 model.
    basin : str
        Name of the ocean basin (e.g., NA, WP, EP, SH).
    expmnt : str
        Name of the model experiment (e.g., historical, ssp585).

    Returns:
    --------
    ds : xarray.Dataset
        The loaded NetCDF dataset.
    lat_trks : numpy.ndarray
        Array of track latitude (degree).
    lon_trks : numpy.ndarray
        Array of track longitude (degree).
    n_trk : numpy.ndarray
        Array of track numbers.
    v_trks : numpy.ndarray
        Array of track azimuthal wind (m/s).
    vmax_trks : numpy.ndarray
        Array of track maximum wind or intensity (m/s).
    u850_trks : numpy.ndarray
        Array of 850 hPa zonal wind components (m/s).
    v850_trks : numpy.ndarray
        Array of 850 hPa meridional wind components (m/s).
    tc_month : numpy.ndarray
        Array of tropical cyclone months.
    tc_years : numpy.ndarray
        Array of tropical cyclone years.
    tc_time : numpy.ndarray
        Array of tropical cyclone times.
    """
    ncfilename = os.path.join(data_directory, expmnt, f"tracks_{basin}_{model}_*.nc")
    try:
        ncfile = glob.glob(ncfilename)[0]
        ds = xr.open_dataset(ncfile)
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
    except IndexError as exc:
        raise FileNotFoundError(f"No files found for the given path: {ncfilename}") from exc
    except KeyError as e:
        raise KeyError(f"Missing expected data in the dataset: {e}") from e
    except Exception as e:
        raise RuntimeError(f"An error occurred while loading the dataset: {e}") from e

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


def load_tracks_GCMs(data_directory=os.path.join(DATA_DIR, 'downscaled'),
                     model='E3SM-1-0', basin='NA', expmnt='historical'):
    """
    Load the downscaled tracks of tropical cyclones from CMIP6 models.

    Parameters:
    -----------
    data_directory : str
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
        Latitudes of tropical cyclone tracks (degree).
    lon_trks : numpy.ndarray
        Longitudes of tropical cyclone tracks (degree).
    year_trks : numpy.ndarray
        Years of tropical cyclones.
    id_trks : numpy.ndarray
        Indices of tropical cyclones.
    vmax_trks : numpy.ndarray
        Maximum wind speeds of tropical cyclones (m/s).
    """
    ncfilename = os.path.join(data_directory, expmnt, f"tracks_{basin}_{model}_*.nc")
    try:
        ncfile = glob.glob(ncfilename)[0]
        ds = xr.open_dataset(ncfile)
        lon_trks = ds['lon_trks'].values
        lat_trks = ds['lat_trks'].values
        vmax_trks = ds['vmax_trks'].values
        year_trks = ds['year'].values
        id_trks = ds['n_trk'].values
    except IndexError as exc:
        raise FileNotFoundError(f"No files found for the given path: {ncfilename}") from exc
    except KeyError as e:
        raise KeyError(f"Missing expected data in the dataset: {e}") from e
    except Exception as e:
        raise RuntimeError(f"An error occurred while loading the dataset: {e}") from e

    return lat_trks, lon_trks, year_trks, id_trks, vmax_trks


def fetch_directory_tree(
    url="https://web.corral.tacc.utexas.edu/setxuifl/tropical_cyclones/downscaled_cmip6_tracks/",
    depth=0,
    max_depth=2,
):
    """
    Fetch and print the directory tree from a specified URL up to a given depth.

    Parameters:
    -----------
    depth : int, optional
        Current depth of the directory tree (default is 0).
    max_depth : int, optional
        Maximum depth to traverse the directory tree (default is 2).

    Returns:
    --------
    None
    """
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
    url="https://web.corral.tacc.utexas.edu/setxuifl/tropical_cyclones/downscaled_cmip6_tracks/",
    experiments=None,
    models=None,
    target_directory=os.path.join(DATA_DIR, "downscaled"),
):
    """
    Download the downscaled tropical cyclone tracks from CMIP6 models.

    Parameters:
    -----------
    experiment : str, optional
        Name of the model experiment (default is 'historical').
    model : str, optional
        Name of the CMIP6 model (default is 'E3SM-1-0').
    target_directory : str, optional
        Path to the target directory to save the downloaded files (default is '{DATA_DIR}/downscaled/').

    Returns:
    --------
    None
    """
    if isinstance(experiments, str):
        experiments = [experiments]
    if isinstance(models, str):
        models = [models]

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
    Get the bounding box from a shapefile.

    Parameters:
    -----------
    shapefile_path : str
        Path to the shapefile.

    Returns:
    --------
    bbox : tuple
        Bounding box in the format (xmin, xmax, ymin, ymax).
    """
    sf = shapefile.Reader(shapefile_path)
    bbox = sf.bbox  # returns [xmin, ymin, xmax, ymax]
    return bbox[0], bbox[2], bbox[1], bbox[3]
