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


def format_mapping(ax, extent=None, shapefile=None, show_gridlabel=False):
    """
    Format and decorate the map by setting extent, grid spacing, and optional shapefile overlay.

    Parameters:
    -----------
    - ax: matplotlib.axes.Axes
        The axis object on which to format the map.
    - extent: tuple, optional
        Bounding box and spacing in data coordinates (left, right, bottom, top, dx, dy).
        Defines the spatial extent of the map.
    - shapefile: str or shapefile-like object, optional
        A shapefile to overlay on the map. Provides additional geographic context.
    - show_gridlabels: bool, optional (default=False)
        Whether to display grid labels on the map. Controls the visibility of gridline labels.

    """
    if extent is None:
        extent = (-180, 180, -90, 90, 10, 10)

    xmin = extent[0]
    xmax = extent[1]
    ymin = extent[2]
    ymax = extent[3]
    dx = extent[4]
    dy = extent[5]

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

    # PLOT US REGIONS
    if shapefile is not None:
        shapereg = ShapelyFeature(shapefile.geometries(), ccrs.PlateCarree(),
                                  facecolor='none', edgecolor='k')
        ax.add_feature(shapereg, facecolor='none', zorder=4, lw=0.75)


def plot_density(ax, lat, lon, density, levels, extent=None, alpha=1.0,
                 cmap='viridis', logscale=False, show_gridlabel=False, shapefile=None,
                 title=None, title_ypos=0.95, title_fontcolor='k', title_fontstyle='regular'):
    """
    This function plots the density of TC tracks

    Parameters:
    -----------
    - ax: matplotlib.axes.Axes
        The axis on which to plot the density map
    - lat: array-like
        Latitude coordinates of the tropical cyclone tracks
    - lon: array-like
        Longitude coordinates of the tropical cyclone tracks
    - density: 2D array-like
        Density values corresponding to the latitude and longitude grid
    - levels: array-like
        Contour levels for plotting density
    - extent: tuple, optional
        Bounding box and spacing in data coordinates (left, right, bottom, top, dx, dy).
        Defines the spatial extent of the map.
    - alpha: float, optional (default=10)
        Transparency level of the density plot
    - cmap: str or Colormap, optional (default='viridis')
        Colormap to use for visualizing density values
    - logscale: bool, optional (default=False)
        Whether to use a logarithmic scale for density values
    - show_gridlabels: bool, optional (default=False)
        Whether to display grid labels on the map
    - shapefile: str or shapefile-like object, optional
        A shapefile to overlay on the map
    - title: str, optional
        Title of the plot
    - title_ypos: float, optional (default=095)
        Y-coordinate for the title position
    - title_fontsize: int, optional (default=12)
        Font size of the title
    - title_fontcolor: str, optional (default='k')
        Color of the plot title
    - title_fontstyle: str, optional (default='regular')
        Font style of the plot title
    """

    # Map format
    format_mapping(ax, extent=extent, shapefile=shapefile, show_gridlabel=show_gridlabel)

    # Contour fill map
    if logscale:
        im = ax.contourf(lon, lat, density, alpha=alpha, cmap=cmap,
                         levels=levels, norm=LogNorm(vmin=levels[0], vmax=levels[-1]),
                         locator=ticker.LogLocator(), extend='both', transform=ccrs.PlateCarree())
    else:
        im = ax.contourf(lon, lat, density, alpha=alpha, cmap=cmap,
                         levels=levels, extend='both', transform=ccrs.PlateCarree())

    if title is not None:
        ax.set_title(title, y=title_ypos, color=title_fontcolor,
                     fontweight=title_fontstyle, fontsize=16)

    return im


def plot_tracks(ax, lats, lons, vmaxs, track_inds, interval=1, wind_speed_threshold=0, extent=None,
                alpha=1.0, cmap='viridis', show_gridlabel=False, shapefile=None, title=None,
                title_ypos=0.95, wind_color=False, title_fontcolor='k', title_fontstyle='regular',
                norm=plt.Normalize(0, 100)):
    """
    This function plots the tracks of TCs

    Parameters:
    -----------
    - ax: matplotlib.axes.Axes
        The axis on which to plot the tracks.
    - lats: array-like
        Latitude data for the tracks.
    - lons: array-like
        Longitude data for the tracks.
    - vmaxs: array-like
        Maximum wind speeds associated with the tracks.
    - track_inds: array-like
        Indices of the tropical cyclone tracks to plot.
    - interval: int, optional (default=1)
        Interval of track indices to plot.        
    - wind_speed_threshold: float, optional (default=0)
        Minimum wind speed required to plot a track.
    - extent: tuple, optional
        Bounding box and spacing in data coordinates (left, right, bottom, top, dx, dy).
        Defines the spatial extent of the map.
    - alpha: float, optional (default=1.0)
        Transparency level of the plot.
    - cmap: str or Colormap, optional (default='viridis')
        Colormap to use for the plot.
    - show_gridlabels: bool, optional (default=False)
        Whether to display grid labels on the map.        
    - shapefile: str or shapefile-like object, optional
        A shapefile to overlay on the map.
    - title: str, optional
        Title of the plot.
    - title_ypos: float, optional (default=0.95)
        Y-coordinate for the title position.
    - wind_color: bool, optional (default=False)
        Whether to color the tracks based on wind speed.
    - title_fontcolor: str, optional (default='black')
        Color of the plot title.
    - title_fontstyle: str, optional (default='regular')
        Font style of the plot title.
    - norm: matplotlib.colors.Normalize, optional (default=plt.Normalize(0, 100))
        Normalization for the colormap scaling.
    """

    # Map format
    format_mapping(ax, extent=extent, shapefile=shapefile, show_gridlabel=show_gridlabel)

    # PLOT HURRICANE TRACKS
    for i in track_inds:
        lat = lats[i, :].data
        indmax = np.where(lat > -90)[0][-1]+1
        lat = lat[:indmax:interval]
        lon = lons[i, :indmax:interval].data
        vs = vmaxs[i, :indmax:interval].data
        ind = np.where(vs >= wind_speed_threshold)[0]
        if len(ind) > 0:
            lat = lat[ind]
            lon = lon[ind]
            vs = vs[ind]
            points = np.array([lon, lat]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            if wind_color:
                lc = LineCollection(segments, cmap=cmap, norm=norm,
                                    alpha=alpha, transform=ccrs.PlateCarree())
            else:
                lc = LineCollection(segments, cmap=plt.cm.gray, norm=norm,
                                    alpha=0.25, transform=ccrs.PlateCarree())
            lc.set_array(vs)
            lc.set_linewidth(1)
            line = ax.add_collection(lc)
            if wind_color:
                ax.plot(lon[0], lat[0], 'ok', ms=3, alpha=1.0,
                        transform=ccrs.PlateCarree(), zorder=1)
            else:
                ax.plot(lon[0], lat[0], 'ok', ms=3, alpha=0.25,
                        transform=ccrs.PlateCarree(), zorder=1)

    if title is not None:
        ax.set_title(title, y=title_ypos, color=title_fontcolor,
                     fontweight=title_fontstyle, fontsize=16)

    return line


def get_tracks_landfall_region(lat_tracks, lon_tracks, num_tracks, polygon):
    """
    Get the indices of tracks that make landfall within a specified polygon.

    Parameters:
    -----------
    - lat_tracks: 2D array-like
        Array of latitudes for each track. Shape is (num_tracks, num_points).
    - lon_tracks: 2D array-like
        Array of longitudes for each track. Shape is (num_tracks, num_points).
    - num_tracks: int
        Number of tracks to sample and check for landfall.
    - polygon: shapely.geometry.Polygon
        Polygon representing the region to check for landfall.

    Returns:
    --------
    - ind_poly: list
        List of indices of tracks that make landfall within the polygon.

    """

    num_tcs = np.shape(lon_tracks)[0]
    ind_tracks = np.random.choice(np.arange(num_tcs), num_tracks, replace=False)
    num_tcs = len(ind_tracks)
    ind_poly = []

    for tc_id in ind_tracks:
        temp = lat_tracks[tc_id, :]
        indmax = np.where(temp > -90)[0][-1]+1

        lat = lat_tracks[tc_id, :indmax]
        lon = lon_tracks[tc_id, :indmax]

        # Convert longitudes to [-180, 180] range
        lon[lon > 180] -= 360

        for i in range(indmax):
            p1 = Point(lon[i], lat[i])
            ans = p1.within(polygon)
            if np.array(ans)[0]:
                ind_poly.append(tc_id)
                break

    return ind_poly


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    """
    Truncate a colormap to a specific range and return a new colormap.

    Parameters:
    -----------
    - cmap: matplotlib.colors.Colormap
        The colormap to be truncated.
    - minval: float, optional (default=0.0)
        Minimum value for the colormap range.
    - maxval: float, optional (default=1.0)
        Maximum value for the colormap range.
    - n: int, optional (default=100)
        Number of colors in the truncated colormap.

    Returns:
    --------
    - LinearSegmentedColormap
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
