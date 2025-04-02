"""
Test suite for rainfall-related functions in PyTCR.
"""

import pytest
import numpy as np
from tcr.rainfall import calculate_rainfall_rate, calculate_etr_swath


@pytest.fixture
def storm_test_data():
    """Returns sample storm track data"""
    nt = 0
    time_steps = 10
    num_storms = 1

    latitude = np.random.uniform(10, 30, (num_storms, time_steps))
    longitude = np.random.uniform(-80, -60, (num_storms, time_steps))
    radius_storm = np.random.uniform(10, 40, (num_storms, time_steps))
    velocity = np.random.uniform(20, 50, (num_storms, time_steps))
    radius_storm_secondary = np.random.uniform(5, 20, (num_storms, time_steps))
    velocity_secondary = np.random.uniform(10, 40, (num_storms, time_steps))
    ut = np.random.uniform(-5, 5, (num_storms, time_steps))
    vt = np.random.uniform(-5, 5, (num_storms, time_steps))
    u850 = np.random.uniform(-10, 10, (num_storms, time_steps))
    v850 = np.random.uniform(-10, 10, (num_storms, time_steps))

    months = np.full((num_storms, time_steps), 8)  # August
    days = np.full((num_storms, time_steps), 15)  # 15th
    hours = np.full((num_storms, time_steps), 12)  # Noon

    return (
        nt, latitude, longitude, radius_storm, velocity,
        radius_storm_secondary, velocity_secondary, ut, vt,
        u850, v850, months, days, hours
    )


def test_calculate_rainfall_rate(storm_test_data):
    """Test that calculate_rainfall_rate runs and returns correct shape"""
    nt, lat, lon, r_storm, vel, r_sec, v_sec, ut, vt, u850, v850, months, days, hours = storm_test_data

    rainrate, x, y = calculate_rainfall_rate(
        nt, lat, lon, r_storm, vel, r_sec, v_sec, ut, vt,
        u850, v850, months, days, hours,
        monthplot=8, dayplot=15, hourplot=12
    )

    assert isinstance(rainrate, np.ndarray)
    assert rainrate.shape[0] > 0 and rainrate.shape[1] > 0  # Non-empty output
    assert isinstance(x, np.ndarray)
    assert isinstance(y, np.ndarray)


def test_calculate_etr_swath(storm_test_data):
    """Test that calculate_etr_swath runs and returns correct shape"""
    nt, lat, lon, r_storm, vel, r_sec, v_sec, ut, vt, u850, v850, _, _, _ = storm_test_data

    x, y, rain_swath = calculate_etr_swath(
        nt, lat, lon, r_storm, vel, r_sec, v_sec, ut, vt, u850, v850
    )

    assert isinstance(rain_swath, np.ndarray)
    # Non-empty output
    assert isinstance(x, np.ndarray)
    assert isinstance(y, np.ndarray)
    assert rain_swath.shape[0] > 0 and rain_swath.shape[1] > 0
    assert rain_swath.shape[0] == len(y)
    assert rain_swath.shape[1] == len(x)
