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
from shapely.geometry import Point


def format_mapping(ax, extent=(-180, 180, -90, 90, 10, 10),
                   shapefile=None, show_gridlabel=False):
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

    ax.coastlines(lw=1)
    ax.set_extent((xmin, xmax, ymin, ymax))
    ax.add_feature(cfeature.LAND, zorder=1, edgecolor='k', lw=0.5, alpha=0.75)
    ax.set_boundary(aoi, transform=ccrs.PlateCarree())

    if show_gridlabel:
        # Find the smallest multiple of 10 greater than or equal to xmin
        xstart = np.ceil(xmin / 10) * 10
        ystart = np.ceil(ymin / 10) * 10

        # Find the largest multiple of 10 less than or equal to xmax
        xend = np.floor(xmax / 10) * 10 + 1
        yend = np.floor(ymax / 10) * 10 + 1

        ax.gridlines(draw_labels=True, xlocs=np.arange(xstart, xend, dx),
                     ylocs=np.arange(ystart, yend, dy), rotate_labels=False,
                     x_inline=False, y_inline=False, lw=0.25, ls='--')

    if shapefile is not None:
        shapereg = ShapelyFeature(shapefile.geometries(), ccrs.PlateCarree(),
                                  facecolor='none', edgecolor='k')
        ax.add_feature(shapereg, facecolor='none', zorder=4, lw=0.75)


def plot_density(ax, lat, lon, density, levels, extent=None, alpha=1,
                 cmap='viridis', logscale=False, show_gridlabel=False,
                 shapefile=None, title=None, title_ypos=1,
                 title_fontcolor='k', title_fontstyle='regular'):
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
        ax, extent=extent, shapefile=shapefile, show_gridlabel=show_gridlabel
    )

    # Contour fill map
    if logscale:
        norm = LogNorm(vmin=levels[0], vmax=levels[-1])
        locator = ticker.LogLocator()
    else:
        norm = None
        locator = None

    im = ax.contourf(lon, lat, density, alpha=alpha, cmap=cmap,
                     levels=levels, norm=norm, locator=locator,
                     extend='both', transform=ccrs.PlateCarree())

    # Set plot title if provided
    if title:
        ax.set_title(title, y=title_ypos, color=title_fontcolor,
                     fontweight=title_fontstyle, fontsize=16)

    return im


def plot_tracks(ax, lats, lons, vmaxs, track_inds, interval=1,
                wind_speed_threshold=0, extent=None, alpha=1,
                cmap='viridis', show_gridlabel=False, shapefile=None,
                title=None, title_ypos=1, wind_color=False,
                title_fontcolor='k', title_fontstyle='regular',
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
        ax, extent=extent, shapefile=shapefile, show_gridlabel=show_gridlabel
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


def get_tracks_landfall_region(lat_tracks, lon_tracks, num_tracks, polygon):
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
    ind_tracks = np.random.choice(
        np.arange(num_tcs), num_tracks, replace=False
    )
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
