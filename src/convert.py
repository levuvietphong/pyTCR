"""
Functions for time and unit conversion in PyTCR
"""

import numpy as np
import pandas as pd


def to_datetime_from_Matlab(years1d, months, days, hours):
    """
    Convert year, month, day, hour in integer to python datetime format
    """
    _, n = months.shape
    years = np.tile(years1d.T, (1, n))

    # Flatten the 2D arrays to work with pandas Dataframe
    years_flat = years.flatten()
    months_flat = months.flatten()
    days_flat = days.flatten()
    hours_flat = hours.flatten()

    df = pd.DataFrame({
        'year': years_flat,
        'month': months_flat,
        'day': days_flat,
        'hour': hours_flat
    })

    # Convert to datetime, handling NaN values
    datetimes_flat = pd.to_datetime(df, errors='coerce')
    datetime_numeric = datetimes_flat.astype(int) / 10**9  # Convert to seconds
    datetime_array = datetime_numeric.to_numpy()

    # Reshape back to the original 2D structure
    datetimes_2d = datetime_array.reshape(years.shape)

    return datetimes_2d


def to_datetime_from_netcdf(years, months, times):
    """
    Convert year, month, time to python datetime format

    Inputs:
    -------
        - years : 1d-array of year info for each storm [1 x nstorms]
        - months : 1d-array of month info for each storm [1 x nstorms]
        - times : 1d-array of time duration (in second) for all storms.
                  The interval is 3600 s. The length of all storm is 15 days.
                  We assume day always starts at the first date of the month

    Returns:
    --------
        - datetime:
    """

    m = len(years)
    n = len(times)

    # Create days and hours from times
    hours = (times / 3600) % 24
    days = (times / 3600) // 24

    years2d = np.tile(years[:, np.newaxis], (1, n))
    months2d = np.tile(months[:, np.newaxis], (1, n))
    days2d = np.tile(days.T, (m, 1))
    hours2d = np.tile(hours.T, (m, 1))

    # Flatten the 2D arrays to work with pandas Dataframe
    years_flat = years2d.flatten()
    months_flat = months2d.flatten()
    days_flat = days2d.flatten()
    hours_flat = hours2d.flatten()

    df = pd.DataFrame({
        'year': years_flat,
        'month': months_flat,
        'day': days_flat,
        'hour': hours_flat
    })

    # Convert to datetime, handling NaN values
    datetimes_flat = pd.to_datetime(df, errors='coerce')
    datetime_numeric = datetimes_flat.astype(int) / 10**9  # Convert to seconds
    datetime_array = datetime_numeric.to_numpy()

    # Reshape back to the original 2D structure
    datetimes_2d = datetime_array.reshape(years2d.shape)

    return datetimes_2d
