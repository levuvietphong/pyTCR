import pytest
import numpy as np
from tcr.wind import windprofiles


@pytest.fixture
def test_data():
    """Fixture to setup test variables for the functions."""
    vm = np.array([50, 60, 70])
    rm = np.array([30, 40, 50])
    r = np.array([10, 20, 30])
    wp = 3
    vm2 = np.array([10, 20, 30])
    rm2 = np.array([10, 15, 20])

    return vm, rm, r, wp, vm2, rm2


def test_windprofiles_basic(test_data):
    """Test windprofiles with basic input values."""
    vm, rm, r, wp, _, _ = test_data
    result = windprofiles(vm, rm, r, wp)
    assert result.shape == r.shape
    assert np.all(result >= 0)


def test_windprofiles_secondary_eyewall(test_data):
    """Test windprofiles with secondary eyewall parameters."""
    vm, rm, r, wp, vm2, rm2 = test_data
    result = windprofiles(vm, rm, r, wp, vm2, rm2)
    assert result.shape == r.shape
    assert np.all(result >= 0)


def test_windprofiles_different_profiles(test_data):
    """Test windprofiles with different wind profile models."""
    vm, rm, r, _, _, _ = test_data
    for wp in [1, 2, 3]:
        result = windprofiles(vm, rm, r, wp)
        assert result.shape == r.shape
        assert np.all(result >= 0)


def test_windprofiles_opt_false(test_data):
    """Test windprofiles with opt=False."""
    vm, rm, r, wp, _, _ = test_data
    result = windprofiles(vm, rm, r, wp, opt=False)
    assert result.shape == r.shape
    assert np.all(result >= 0)


def test_windprofiles_large_radii(test_data):
    """Test windprofiles with large radii values."""
    vm, rm, _, wp, _, _ = test_data
    r = np.array([100, 200, 300])
    result = windprofiles(vm, rm, r, wp)
    assert result.shape == r.shape
    assert np.all(result >= 0)
