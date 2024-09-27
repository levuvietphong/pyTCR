"""
Functions for time and unit conversion in PyTCR
"""

import numpy as np
import pandas as pd


def to_datetime_from_Matlab(years1d, months, days, hours):
    """
    Convert year, month, day, and hour from integer format to Python datetime

    Parameters:
    -----------
    years1d : numpy.ndarray
        1D array containing the year information.
    months : numpy.ndarray
        2D array containing the month information.
    days : numpy.ndarray
        2D array containing the day information.
    hours : numpy.ndarray
        2D array containing the hour information.

    Returns:
    --------
    numpy.ndarray
        2D array of datetime values in seconds since the epoch.
    """
    # Ensure the input arrays have compatible shapes
    if months.shape != days.shape or months.shape != hours.shape:
        raise ValueError(
            "Input arrays 'months', 'days', and 'hours' must have the same shape."
        )

    # Create a 2D array of years by repeating the 1D array across columns
    _, n = months.shape
    years = np.tile(years1d.T, (1, n))

    # Flatten the 2D arrays to work with pandas DataFrame
    years_flat = years.flatten()
    months_flat = months.flatten()
    days_flat = days.flatten()
    hours_flat = hours.flatten()

    # Create a DataFrame from the flattened arrays
    df = pd.DataFrame(
        {"year": years_flat, "month": months_flat, "day": days_flat, "hour": hours_flat}
    )

    # Convert to datetime, handling NaN values
    datetimes_flat = pd.to_datetime(df, errors="coerce")
    datetime_numeric = (
        datetimes_flat.view("int64") // 10**9
    )  # Convert to seconds since the epoch
    datetime_array = datetime_numeric.to_numpy()

    # Reshape back to the original 2D structure
    datetimes_2d = datetime_array.reshape(years.shape)

    return datetimes_2d


def to_datetime_from_netcdf(years, months, times):
    """
    Convert year, month, and time duration to Python datetime format.

    Parameters:
    -----------
    years : numpy.ndarray
        1D array containing the year information for each storm [1 x nstorms].
    months : numpy.ndarray
        1D array containing the month information for each storm [1 x nstorms].
    times : numpy.ndarray
        1D array containing the time duration (in seconds) for all storms.
        The interval is 3600 seconds, and the length of all storms is 15 days.
        It is assumed that the day always starts at the first date of the month.

    Returns:
    --------
    numpy.ndarray
        2D array of datetime values corresponding to the input years, months, and times.
    """

    num_storms = len(years)
    num_times = len(times)

    # Calculate days and hours from times
    hours = (times / 3600) % 24
    days = (times / 3600) // 24 + 1

    # Create 2D arrays for years, months, days, and hours
    years_2d = np.tile(years[:, np.newaxis], (1, num_times))
    months_2d = np.tile(months[:, np.newaxis], (1, num_times))
    days_2d = np.tile(days.T, (num_storms, 1))
    hours_2d = np.tile(hours.T, (num_storms, 1))

    # Flatten the 2D arrays to work with pandas DataFrame
    years_flat = years_2d.flatten()
    months_flat = months_2d.flatten()
    days_flat = days_2d.flatten()
    hours_flat = hours_2d.flatten()

    # Create a DataFrame and convert to datetime, handling NaN values
    df = pd.DataFrame(
        {"year": years_flat, "month": months_flat, "day": days_flat, "hour": hours_flat}
    )
    datetimes_flat = pd.to_datetime(df, errors="coerce")
    datetime_numeric = datetimes_flat.astype(int) / 10**9  # Convert to seconds
    datetime_array = datetime_numeric.to_numpy()

    # Reshape back to the original 2D structure
    datetimes_2d = datetime_array.reshape(years_2d.shape)

    return datetimes_2d
