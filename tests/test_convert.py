"""
Test suite for unit conversion in PyTCR.
"""

import numpy as np
import pytest
from tcr.iodata import convert_to_mps
from tcr.convert import to_datetime_from_netcdf


def test_knots_to_mps():
    """Test the distance unit conversion."""
    assert convert_to_mps(1) == pytest.approx(0.5144444)
    assert convert_to_mps(0) == 0
    assert convert_to_mps(10) == pytest.approx(5.144444)


def test_datetime():
    """Test datetime conversion from netcdf to python format."""
    years = np.array([1990, 1991])
    months = np.array([9, 10])
    times = np.array([3600, 7200])
    outputs_obs = np.array([[6.521508e08, 6.521544e08], [6.862788e08, 6.862824e08]])
    outputs_mod = to_datetime_from_netcdf(years, months, times)
    np.testing.assert_allclose(outputs_mod, outputs_obs, rtol=1e-8)
