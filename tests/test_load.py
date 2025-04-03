"""
Test suite for loading functions in PyTCR.
"""

from unittest.mock import patch
import xarray as xr
import pytest
import numpy as np
from tcr.iodata import load_netcdf_track_data  # Replace with actual module name


@pytest.fixture
def mock_netcdf_data():
    """Fixture to create a mock xarray.Dataset similar to a real TC NetCDF file."""
    data = {
        "lat_trks": (["tracks"], np.array([10.0, 20.0, 30.0])),
        "lon_trks": (["tracks"], np.array([-50.0, -60.0, -70.0])),
        "n_trk": (["tracks"], np.array([1, 2, 3])),
        "v_trks": (["tracks"], np.array([5.5, 6.5, 7.5])),
        "vmax_trks": (["tracks"], np.array([50, 60, 70])),
        "u850_trks": (["tracks"], np.array([1.0, 2.0, 3.0])),
        "v850_trks": (["tracks"], np.array([0.5, 1.5, 2.5])),
        "tc_month": (["tracks"], np.array([6, 7, 8])),
        "tc_years": (["tracks"], np.array([2000, 2001, 2002])),
        "time": (["tracks"], np.array([100, 200, 300])),
    }
    return xr.Dataset(data)


@patch("glob.glob", return_value=["/path/tracks_basin_model_time.nc"])
@patch("xarray.open_dataset")
def test_load_netcdf_track_data(mock_open_dataset, mock_glob, mock_netcdf_data):
    """Test successful loading of a NetCDF file with expected variables."""
    mock_open_dataset.return_value = mock_netcdf_data

    (ds, lat_trks, lon_trks, n_trk, v_trks, vmax_trks, u850_trks, v850_trks,
     tc_month, tc_years, tc_time) = load_netcdf_track_data("/path")

    # Assertions
    assert isinstance(ds, xr.Dataset)
    assert isinstance(lat_trks, np.ndarray) and lat_trks.shape == (3,)
    assert isinstance(lon_trks, np.ndarray) and lon_trks.shape == (3,)
    assert isinstance(n_trk, np.ndarray) and n_trk.shape == (3,)
    assert isinstance(v_trks, np.ndarray) and v_trks.shape == (3,)
    assert isinstance(vmax_trks, np.ndarray) and vmax_trks.shape == (3,)
    assert isinstance(u850_trks, np.ndarray) and u850_trks.shape == (3,)
    assert isinstance(v850_trks, np.ndarray) and v850_trks.shape == (3,)
    assert isinstance(tc_month, np.ndarray) and tc_month.shape == (3,)
    assert isinstance(tc_years, np.ndarray) and tc_years.shape == (3,)
    assert isinstance(tc_time, np.ndarray) and tc_time.shape == (3,)


@patch("glob.glob", return_value=[])
def test_file_not_found(mock_glob):
    """Test that FileNotFoundError is raised when no NetCDF file is found."""
    with pytest.raises(FileNotFoundError, match="No files found for the given path"):
        load_netcdf_track_data("/path")


@patch("glob.glob", return_value=["/path/tracks_NA_E3SM-1-0_001.nc"])
@patch("xarray.open_dataset")
def test_missing_variable(mock_open_dataset, mock_glob):
    """Test that KeyError is raised when a required variable is missing."""
    mock_dataset = xr.Dataset(
        {"lat_trks": ("tracks", np.array([10.0, 20.0, 30.0]))})
    mock_open_dataset.return_value = mock_dataset

    with pytest.raises(KeyError, match="Missing expected data in the dataset"):
        load_netcdf_track_data("/path")
