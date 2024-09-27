"""
Functions for plottings in PyTCR
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import ticker
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm, LinearSegmentedColormap

from cartopy.feature import ShapelyFeature
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Point, Polygon


def format_mapping(ax, extent=(-180, 180, -90, 90, 10, 10),
                   shapefile=None, show_gridlabel=False, show_coastlines=True,
                   add_features=True, utm=False):
    """
    Format and decorate the map by setting extent, grid spacing, and optional
    shapefile overlay.

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axis object on which to format the map.
    extent : tuple, optional
        Bounding box and spacing in data coordinates (left, right, bottom, 
        top, dx, dy). Defines the spatial extent of the map. Default is 
        (-180, 180, -90, 90, 10, 10).
    shapefile : str or shapefile-like object, optional
        A shapefile to overlay on the map. Provides additional geographic 
        context.
    show_gridlabel : bool, optional
        Whether to display grid labels on the map. Controls the visibility of 
        gridline labels. Default is False.
    """

    xmin, xmax, ymin, ymax, dx, dy = extent

    n = 60
    aoi = mpath.Path(
        list(zip(np.linspace(xmin, xmax, n), np.full(n, ymax))) +
        list(zip(np.full(n, xmax), np.linspace(ymax, ymin, n))) +
        list(zip(np.linspace(xmax, xmin, n), np.full(n, ymin))) +
        list(zip(np.full(n, xmin), np.linspace(ymin, ymax, n)))
    )

    if not utm:
        ax.set_extent((xmin, xmax, ymin, ymax))
        ax.set_boundary(aoi, transform=ccrs.PlateCarree())

    if show_coastlines:
        ax.coastlines(lw=1)
    if add_features:
        ax.add_feature(cfeature.LAND, zorder=1, edgecolor='k', lw=0.5,
                       alpha=0.75)

    if show_gridlabel:
        if utm:
            ax.gridlines(draw_labels=True, rotate_labels=False,
                         x_inline=False, y_inline=False, lw=0.25, ls='--',
                         crs=ccrs.PlateCarree())
        else:  # lat-lon grid
            # Find the smallest multiple of 10 greater than or equal to xmin
            xstart = np.ceil(xmin / 10) * 10
            ystart = np.ceil(ymin / 10) * 10

            # Find the largest multiple of 10 less than or equal to xmax
            xend = np.floor(xmax / 10) * 10 + 1
            yend = np.floor(ymax / 10) * 10 + 1

            ax.gridlines(
                draw_labels=True,
                xlocs=np.arange(xstart, xend, dx),
                ylocs=np.arange(ystart, yend, dy),
                rotate_labels=False,
                x_inline=False, y_inline=False,
                lw=0.25, ls='--'
            )

    if shapefile is not None:
        shapereg = ShapelyFeature(shapefile.geometries(), ccrs.PlateCarree(),
                                  facecolor='none', edgecolor='k')
        ax.add_feature(shapereg, facecolor='none', zorder=4, lw=0.75)


def plot_density(ax, lat, lon, density, levels, extent=None, alpha=1,
                 cmap='viridis', logscale=False, show_gridlabel=False,
                 show_coastlines=True, add_features=True,
                 shapefile=None, title=None, title_ypos=1,
                 title_fontcolor='k', title_fontstyle='regular', 
                 method='contourf', utm=False):
    """
    Plot the density of tropical cyclone (TC) tracks.

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axis on which to plot the density map.
    lat : array-like
        Latitude coordinates of the tropical cyclone tracks.
    lon : array-like
        Longitude coordinates of the tropical cyclone tracks.
    density : 2D array-like
        Density values corresponding to the latitude and longitude grid.
    levels : array-like
        Contour levels for plotting density.
    extent : tuple, optional
        Bounding box and spacing in data coordinates (left, right, bottom, 
        top, dx, dy). Defines the spatial extent of the map. Default is None.
    alpha : float, optional
        Transparency level of the density plot. Default is 1.0.
    cmap : str or Colormap, optional
        Colormap to use for visualizing density values. Default is 'viridis'.
    logscale : bool, optional
        Whether to use a logarithmic scale for density values. Default is 
        False.
    show_gridlabel : bool, optional
        Whether to display grid labels on the map. Default is False.
    shapefile : str or shapefile-like object, optional
        A shapefile to overlay on the map. Default is None.
    title : str, optional
        Title of the plot. Default is None.
    title_ypos : float, optional
        Y-coordinate for the title position. Default is 1.0.
    title_fontcolor : str, optional
        Color of the plot title. Default is 'k'.
    title_fontstyle : str, optional
        Font style of the plot title. Default is 'regular'.
    """

    # Format the map
    format_mapping(
        ax, extent=extent, shapefile=shapefile, show_gridlabel=show_gridlabel,
        show_coastlines=show_coastlines, add_features=add_features, utm=utm,
    )

    # Contour fill map
    if logscale:
        norm = LogNorm(vmin=levels[0], vmax=levels[-1])
        locator = ticker.LogLocator()
    else:
        norm = None
        locator = None

    if method == 'contourf':
        im = ax.contourf(lon, lat, density, alpha=alpha, cmap=cmap,
                         levels=levels, norm=norm, locator=locator,
                         extend='both', transform=ccrs.PlateCarree())
    elif method == 'pcolormesh':
        im = ax.pcolormesh(lon, lat, density, alpha=alpha, cmap=cmap,
                           vmin=levels[0], vmax=levels[-1], norm=norm,
                           transform=ccrs.PlateCarree())
    else:
        raise ValueError(f"Invalid method: {method}")

    # Set plot title if provided
    if title:
        ax.set_title(title, y=title_ypos, color=title_fontcolor,
                     fontweight=title_fontstyle, fontsize=16)

    return im


def plot_tracks(ax, lats, lons, vmaxs, track_inds, interval=1,
                wind_speed_threshold=0, extent=None, alpha=1,
                cmap='viridis', show_gridlabel=False, show_coastlines=True,
                add_features=True, shapefile=None, title=None,
                title_ypos=1, wind_color=False, title_fontcolor='k',
                title_fontstyle='regular',
                norm=plt.Normalize(0, 100)):
    """
    Plot the tracks of tropical cyclones (TCs).

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axis on which to plot the tracks.
    lats : array-like
        Latitude data for the tracks.
    lons : array-like
        Longitude data for the tracks.
    vmaxs : array-like
        Maximum wind speeds associated with the tracks.
    track_inds : array-like
        Indices of the tropical cyclone tracks to plot.
    interval : int, optional
        Interval of track indices to plot. Default is 1.
    wind_speed_threshold : float, optional
        Minimum wind speed required to plot a track. Default is 0.
    extent : tuple, optional
        Bounding box and spacing in data coordinates (left, right, bottom, top,
        dx, dy). Defines the spatial extent of the map. Default is None.
    alpha : float, optional
        Transparency level of the plot. Default is 1.0.
    cmap : str or Colormap, optional
        Colormap to use for the plot. Default is 'viridis'.
    show_gridlabel : bool, optional
        Whether to display grid labels on the map. Default is False.
    shapefile : str or shapefile-like object, optional
        A shapefile to overlay on the map. Default is None.
    title : str, optional
        Title of the plot. Default is None.
    title_ypos : float, optional
        Y-coordinate for the title position. Default is 1.0.
    wind_color : bool, optional
        Whether to color the tracks based on wind speed. Default is False.
    title_fontcolor : str, optional
        Color of the plot title. Default is 'k'.
    title_fontstyle : str, optional
        Font style of the plot title. Default is 'regular'.
    norm : matplotlib.colors.Normalize, optional
        Normalization for the colormap scaling.

    Returns:
    --------
    line : matplotlib.collections.LineCollection
        The LineCollection object representing the plotted tracks.
    """

    # Format the map
    format_mapping(
        ax, extent=extent, shapefile=shapefile, show_gridlabel=show_gridlabel,
        show_coastlines=show_coastlines, add_features=add_features
    )

    # Plot hurricane tracks
    for i in track_inds:
        lat = lats[i, :]
        indmax = np.where(lat > -90)[0][-1] + 1
        lat = lat[:indmax:interval]
        lon = lons[i, :indmax:interval]
        vs = vmaxs[i, :indmax:interval]
        ind = np.where(vs >= wind_speed_threshold)[0]
        if len(ind) > 0:
            lat = lat[ind]
            lon = lon[ind]
            vs = vs[ind]
            # points = np.array([lon, lat]).T.reshape(-1, 1, 2)
            points = np.column_stack((lon, lat)).reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = LineCollection(
                segments,
                cmap=cmap if wind_color else plt.cm.gray,
                norm=norm,
                alpha=alpha if wind_color else 0.25,
                transform=ccrs.PlateCarree(),
            )
            lc.set_array(vs)
            lc.set_linewidth(1)
            line = ax.add_collection(lc)
            ax.plot(
                lon[0], lat[0], "ok", ms=3,
                alpha=alpha if wind_color else 0.25,
                transform=ccrs.PlateCarree(), zorder=1,
            )

    if title:
        ax.set_title(title, y=title_ypos, color=title_fontcolor,
                     fontweight=title_fontstyle, fontsize=16)

    return line


def plot_exceedance_probability(
    data, ax, ylabel="Value", xlabel="Exceedance Probability", title=None,
    fontweight='regular', fontsize=12
):
    """
    Plots the exceedance probability of the given data.

    Parameters:
    - data: array-like, the data for which to plot the exceedance probability
    - xlabel: str, label for the x-axis (default is 'Value')
    - ylabel: str, label for the y-axis (default is 'Exceedance Probability')
    - title: str, title of the plot
    """
    # Ensure data is a NumPy array
    data = np.array(data)
    num_storms = len(data)

    # Sort the data in descending order
    sorted_data = np.sort(data)[::-1]

    # Calculate exceedance probability
    n = len(sorted_data)
    exceedance_prob = np.arange(1, n + 1) / n  # Ranks divided by total count

    ax.plot(exceedance_prob, sorted_data, "o", color="tab:blue",
            ms=3, label="Exceedance Probability")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f'{title} (# storms: {num_storms})', fontweight=fontweight,
                 fontsize=fontsize)
    ax.set_xscale('log')
    ax.set_xlim([1e-3, 1e0])
    ax.set_ylim(bottom=0)
    ax.grid(True, linestyle='--', alpha=0.5, color='k')


def get_tracks_landfall_region(lat_tracks, lon_tracks, polygon, num_tracks=None):
    """
    Get the indices of tracks that make landfall within a specified polygon.

    Parameters:
    -----------
    lat_tracks : numpy.ndarray
        2D array of latitudes for each track [num_tracks, num_points].
    lon_tracks : numpy.ndarray
        2D array of longitudes for each track [num_tracks, num_points].
    num_tracks : int
        Number of tracks to sample and check for landfall.
    polygon : shapely.geometry.Polygon
        Polygon representing the region to check for landfall.

    Returns:
    --------
    ind_poly : list
        List of indices of tracks that make landfall within the polygon.
    """

    num_tcs = lat_tracks.shape[0]
    if num_tracks is not None:
        ind_tracks = np.sort(
            np.random.choice(np.arange(num_tcs), num_tracks, replace=False)
        )
    else:
        ind_tracks = np.arange(num_tcs)

    ind_poly = []

    for tc_id in ind_tracks:
        temp = lat_tracks[tc_id, :]
        indmax = np.where(temp > -90)[0][-1] + 1

        lat = lat_tracks[tc_id, :indmax]
        lon = lon_tracks[tc_id, :indmax]

        # Convert longitudes to [-180, 180] range
        lon[lon > 180] -= 360

        for i in range(indmax):
            point = Point(lon[i], lat[i])
            if point.within(polygon):
                ind_poly.append(tc_id)
                break

    return ind_poly


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    """
    Truncate a colormap to a specific range and return a new colormap.

    Parameters:
    -----------
    cmap : matplotlib.colors.Colormap
        The colormap to be truncated.
    minval : float, optional
        Minimum value for the colormap range. Default is 0.0.
    maxval : float, optional
        Maximum value for the colormap range. Default is 1.0.
    n : int, optional
        Number of colors in the truncated colormap. Default is 100.

    Returns:
    --------
    new_cmap : matplotlib.colors.LinearSegmentedColormap
        A new colormap truncated to the specified range.
    """
    # Ensure minval and maxval are within the range [0, 1]
    minval = max(0.0, minval)
    maxval = min(1.0, maxval)

    # Create new colormap by sampling the original colormap
    new_cmap = LinearSegmentedColormap.from_list(
        f"trunc({cmap.name},{minval:.2f},{maxval:.2f})",
        cmap(np.linspace(minval, maxval, n))
    )

    return new_cmap


def create_buffer_around_POI(lat_poi, lon_poi, radius=1.0):
    """
    Creates a buffer zone around a Point of Interest (POI) defined by its
    latitude and longitude.

    Parameters:
    -----------
    lat_poi : float
        Latitude of the Point of Interest.
    lon_poi : float
        Longitude of the Point of Interest.
    radius : float, optional
        Radius of the buffer zone in degrees. Default is 1.0.

    Returns:
    --------
    Polygon
        A Polygon object representing the buffer zone around the POI.

    Raises:
    ------
    ValueError
        If radius is not a positive number.
    """
    if radius <= 0:
        raise ValueError("Radius must be a positive number.")
    point = Point(lon_poi, lat_poi)
    buffer = point.buffer(radius)
    return Polygon(buffer.exterior.coords)
