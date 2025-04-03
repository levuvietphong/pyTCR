"""
Test suite for physics-related functions in PyTCR.
"""

import numpy as np
import pytest
from tcr.physics import estimate_track_density, estimate_pdi


@pytest.fixture
def test_data():
    """
    Fixture to setup test variables for the functions.
    """
    lat_trks = np.array([[10, 11, 12], [20, 21, 22], [30, 31, 32]])
    lon_trks = np.array([[100, 101, 102], [110, 111, 112], [120, 121, 122]])
    vmax_trks = np.array([[10, 20, 30], [20, 30, 40], [30, 40, 50]])
    num_trks = 2
    threshold = 15
    cellsize = 1
    interval = 1
    dt = 1

    return lat_trks, lon_trks, vmax_trks, num_trks, threshold, cellsize, interval, dt

@pytest.fixture(autouse=True)
def set_seed():
    """Ensure consistent test outputs by setting a random seed."""
    np.random.seed(42)


def test_estimate_track_density(test_data):
    """Test that track_density runs and returns correct shape and values"""
    lat_trks, lon_trks, vmax_trks, num_trks, threshold, cellsize, interval, _ = test_data

    latitude, longitude, density_all = estimate_track_density(
        lat_trks, lon_trks, vmax_trks, num_trks, threshold, cellsize, interval
    )

    # Validate the shape of the output
    assert latitude.shape == (180, )
    assert longitude.shape == (360, )
    assert density_all.shape == (180, 360)

    # Test if density_all produces expected values
    assert np.sum(density_all) == 5


def test_estimate_pdi(test_data):
    """Test that estimate_pdi runs and returns correct shape and values"""
    lat_trks, lon_trks, vmax_trks, num_trks, _, cellsize, _, dt = test_data

    latitude, longitude, pdi_all = estimate_pdi(
        lat_trks, lon_trks, vmax_trks, num_trks, cellsize, dt
    )

    # Validate the shape of the output
    assert latitude.shape == (180, )
    assert longitude.shape == (360, )
    assert pdi_all.shape == (180, 360)

    # Test if pdi_all produces expected values
    assert np.sum(pdi_all) == 486e6
