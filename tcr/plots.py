"""
Functions for plottings in PyTCR
"""

import os
import colorsys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import ticker
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm, LinearSegmentedColormap
import matplotlib.colors as mcolors

import geopandas as gpd
from cartopy.feature import ShapelyFeature
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Point, Polygon


def format_mapping(ax, extent=(-180, 180, -90, 90), projection=None,
                   shapefile=None, gridlabel=False, gridspace=10,
                   coastlines=True, add_features=True):
    """
    Customizes the map by setting its spatial extent, grid spacing, and
    optionally overlays a shapefile.

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axis object where the map will be formatted.
    extent : tuple, optional
        Specifies the map's spatial extent in data coordinates (left, right,
        bottom, top). Defaults to (-180, 180, -90, 90).
    shapefile : str or shapefile-like object, optional
        Path to a shapefile or shapefile object to overlay on the map.
    gridlabel : bool, optional
        Toggles the display of grid labels on the map.
        Default is False, hiding gridline labels.
    projection : cartopy.crs, optional
        The projection to use for the map. Defaults to None, which implies a
        PlateCarree projection.
    coastlines : bool, optional
        Toggles the display of coastlines on the map.
        Default is True, showing coastlines.
    add_features : bool, optional
        Toggles the display of additional geographical features.
        Default is True, showing these features.
    gridspace : int, optional
        Specifies the spacing between gridlines. Defaults to 10 degrees.
    """

    xmin, xmax, ymin, ymax = extent
    if isinstance(projection, ccrs.UTM):
        utm = True
    else:
        utm = False

    if utm:
        utm_extent = projection.transform_points(
            ccrs.PlateCarree(), np.array([xmin, xmax]), np.array([ymin, ymax])
        )
        ax.set_extent(
            (
                utm_extent[0, 0],
                utm_extent[1, 0],
                utm_extent[0, 1],
                utm_extent[1, 1],
            ),
            crs=projection,
        )
    else:
        ax.set_extent((xmin, xmax, ymin, ymax))
        n = 60
        aoi = mpath.Path(
            list(zip(np.linspace(xmin, xmax, n), np.full(n, ymax))) +
            list(zip(np.full(n, xmax), np.linspace(ymax, ymin, n))) +
            list(zip(np.linspace(xmax, xmin, n), np.full(n, ymin))) +
            list(zip(np.full(n, xmin), np.linspace(ymin, ymax, n)))
        )
        ax.set_boundary(aoi, transform=ccrs.PlateCarree())

    if coastlines:
        ax.coastlines(lw=1)

    if add_features:
        ax.add_feature(cfeature.LAND, zorder=1, edgecolor='k', lw=0.5,
                       alpha=0.75)

    if gridlabel:
        if utm:
            ax.gridlines(draw_labels=True, rotate_labels=False,
                         x_inline=False, y_inline=False, lw=0.25, ls='--',
                         crs=ccrs.PlateCarree())
        else:  # lat-lon grid
            # Find the smallest multiple of 10 greater than or equal to xmin
            xstart = np.floor(xmin / gridspace) * gridspace
            ystart = np.floor(ymin / gridspace) * gridspace

            # Find the largest multiple of gridspace less than or equal to xmax
            xend = np.ceil(xmax / gridspace) * gridspace + gridspace/2.
            yend = np.ceil(ymax / gridspace) * gridspace + gridspace/2.

            ax.gridlines(
                draw_labels=True,
                xlocs=np.arange(xstart, xend, gridspace),
                ylocs=np.arange(ystart, yend, gridspace),
                rotate_labels=False, x_inline=False, y_inline=False,
                lw=0.25, ls='--', transform=ccrs.PlateCarree()
            )

    if shapefile is not None:
        shapereg = ShapelyFeature(shapefile.geometries(), ccrs.PlateCarree(),
                                  facecolor='none', edgecolor='k')
        ax.add_feature(shapereg, facecolor='none', zorder=4, lw=0.75)


def plot_density(ax, lat, lon, density, levels, extent=None, projection=None,
                 alpha=1, cmap='viridis', logscale=False, gridlabel=False,
                 gridspace=10, coastlines=True, add_features=True,
                 shapefile=None, title=None, title_ypos=1,
                 title_fontcolor='k', title_fontstyle='regular',
                 method='contourf'):
    """
    Plot the density of tropical cyclone (TC) tracks.

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axis on which to plot the density map
    lat : array-like
        Latitude coordinates of the tropical cyclone tracks (degree)
    lon : array-like
        Longitude coordinates of the tropical cyclone tracks (degree)
    density : 2D array-like
        Density values corresponding to the latitude and longitude grid
    levels : array-like
        Contour levels for plotting density
    extent : tuple, optional
        Bounding box and spacing in data coordinates (left, right, bottom, top)
        Defines the spatial extent of the map, default is None
    projection : str or Cartopy projection, optional
        The projection to use for the map
        Default is None, which implies a PlateCarree projection.
    alpha : float, optional
        Transparency level of the density plot, default is 1.0
    cmap : str or Colormap, optional
        Colormap to use for visualizing density values, default is 'viridis'
    logscale : bool, optional
        Whether to use a logarithmic scale for density values, default is False
    gridlabel : bool, optional
        Whether to display grid labels on the map, default is False
    gridspace : int, optional
        The spacing between gridlines, default is 10
    coastlines : bool, optional
        Whether to display coastlines on the map, default is True
    add_features : bool, optional
        Whether to add additional features to the map, default is True
    shapefile : str or shapefile-like object, optional
        A shapefile to overlay on the map, default is None
    title : str, optional
        Title of the plot, default is None
    title_ypos : float, optional
        Y-coordinate for the title position, default is 1.0
    title_fontcolor : str, optional
        Color of the plot title, default is 'k'
    title_fontstyle : str, optional
        Font style of the plot title, default is 'regular'
    """

    # Format the map
    format_mapping(
        ax, extent=extent, projection=projection, shapefile=shapefile,
        gridlabel=gridlabel, gridspace=gridspace, coastlines=coastlines,
        add_features=add_features
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
        raise ValueError(f"Invalid plotting method: {method}")

    # Set plot title if provided
    if title:
        ax.set_title(title, y=title_ypos, color=title_fontcolor,
                     fontweight=title_fontstyle, fontsize=16)

    return im


def plot_tracks(ax, lats, lons, vmaxs, track_inds, interval=1,
                wind_speed_threshold=0, extent=None, projection=None,
                alpha=1, cmap='viridis', gridlabel=False, gridspace=10,
                coastlines=True, add_features=True, shapefile=None,
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
        Latitude data for the tracks (degrees).
    lons : array-like
        Longitude data for the tracks (degrees).
    vmaxs : array-like
        Maximum wind speeds associated with the tracks (m/s).
    track_inds : array-like
        Indices of the tropical cyclone tracks to plot.
    interval : int, optional
        Interval of track indices to plot. Default is 1.
    wind_speed_threshold : float, optional
        Minimum wind speed required to plot a track (m/s).
    extent : tuple, optional
        Bounding box and spacing (degree) in data coordinates
        (left, right, bottom, top, dx, dy). Defines the spatial extent of the map.
        Default is None.
    projection : str, optional
        Projection to use for the map.  Default is None, which implies a
        PlateCarree projection.
    alpha : float, optional
        Transparency level of the plot. Default is 1.0.
    cmap : str or Colormap, optional
        Colormap to use for the plot. Default is 'viridis'.
    gridlabel : bool, optional
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
    coastlines : bool, optional
        Whether to display coastlines on the map. Default is True.
    add_features : bool, optional
        Whether to add additional features such as country borders and lakes.
        Default is True.
    gridlines : bool, optional
        Whether to display gridlines on the map. Default is False.
    gridspaces : list, optional
        List of grid spacings for the x and y axes. Default is [10, 10].

    Returns:
    --------
    line : matplotlib.collections.LineCollection
        The LineCollection object representing the plotted tracks.
    """

    # Format the map
    format_mapping(
        ax, extent=extent, projection=projection, shapefile=shapefile,
        gridlabel=gridlabel, gridspace=gridspace, coastlines=coastlines,
        add_features=add_features
    )

    # Plot hurricane tracks
    line = None  # Initialize line variable
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
    fontweight='regular', fontsize=12, color="tab:blue"
):
    """
    Plots the exceedance probability of the given data.

    Parameters:
    - data: array-like, the data for which to plot the exceedance probability
    - ax: matplotlib axes object, the axes on which to plot
    - xlabel: str, label for the x-axis (default is 'Exceedance Probability')
    - ylabel: str, label for the y-axis (default is 'Value')
    - title: str, title of the plot
    - fontweight: str, font weight for the title (default is 'regular')
    - fontsize: int, font size for the title (default is 12)
    """
    # Ensure data is a NumPy array
    data = np.array(data)
    num_storms = len(data)

    # Sort the data in descending order
    sorted_data = np.sort(data)[::-1]

    # Calculate exceedance probability
    n = len(sorted_data)
    exceedance_prob = np.arange(1, n + 1) / n  # Ranks divided by total count

    ax.plot(exceedance_prob, sorted_data, "o", color=color,
            ms=3, label="Exceedance Probability")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f'{title} (# storms: {num_storms})', fontweight=fontweight,
                 fontsize=fontsize)
    ax.set_xscale('log')
    ax.set_xlim([1e-3, 1e0])
    ax.set_ylim(bottom=0)
    ax.grid(True, linestyle='--', alpha=0.5, color='k')


def get_tracks_landfall_region(lat_tracks, lon_tracks, polygon, buffer=0,
                               num_tracks=None):
    """
    Get the indices of tracks that make landfall within a specified polygon.

    Parameters:
    -----------
    lat_tracks : numpy.ndarray
        2D array of latitudes (degrees) for each track [num_tracks, num_points]
    lon_tracks : numpy.ndarray
        2D array of longitudes (degrees) for each track [num_tracks, num_points]
    num_tracks : int
        Number of tracks to sample and check for landfall
    buffer : float, optional
        Buffer distance around the polygon (degrees)
    polygon : shapely.geometry.Polygon
        Polygon representing the region to check for landfall

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

    polygon_buffer = polygon.buffer(buffer)
    ind_poly = []

    for tc_id in ind_tracks:
        temp = lat_tracks[tc_id, :]
        indmax = np.where(temp > -90)[0][-1] + 1

        lat = np.array(lat_tracks[tc_id, :indmax])
        lon = np.array(lon_tracks[tc_id, :indmax])

        # Convert longitudes to [-180, 180] range
        lon[lon > 180] -= 360

        for i in range(indmax):
            point = Point(lon[i], lat[i])
            if point.within(polygon_buffer):
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
        Latitude of the Point of Interest (degree)
    lon_poi : float
        Longitude of the Point of Interest (degree)
    radius : float, optional
        Radius of the buffer zone in degrees, default is 1.0

    Returns:
    --------
    Polygon
        A Polygon object representing the buffer zone around the POI
    """
    if radius <= 0:
        raise ValueError("Radius must be a positive number.")

    if lon_poi > 180:
        lon_poi -= 360
    elif lon_poi < -180:
        raise ValueError("Longitude must be greater than or equal to -180 degrees.")

    point = Point(lon_poi, lat_poi)
    buffer = point.buffer(radius)
    return Polygon(buffer.exterior.coords)


def gmtColormap_openfile(cptf, name=None):
    """Read a GMT color map from an OPEN cpt file

    Parameters
    ----------
    cptf : open file or url handle
        path to .cpt file
    name : str, optional
        name for color map
        if not provided, the file name will be used
    """
    # generate cmap name
    if name is None:
        name = '_'.join(os.path.basename(cptf.name).split('.')[:-1])

    # process file
    x = []
    r = []
    g = []
    b = []
    lastls = None
    colorModel = "RGB"  # Set default value
    for line in cptf.readlines():
        ls = line.split()

        # skip empty lines
        if not ls:
            continue

        # parse header info
        if ls[0] in ["#", b"#"]:
            if ls[-1] in ["HSV", b"HSV"]:
                colorModel = "HSV"
            else:
                colorModel = "RGB"
            continue

        # skip BFN info
        if ls[0] in ["B", b"B", "F", b"F", "N", b"N"]:
            continue

        # parse color vectors
        x.append(float(ls[0]))
        r.append(float(ls[1]))
        g.append(float(ls[2]))
        b.append(float(ls[3]))

        # save last row
        lastls = ls

    x.append(float(lastls[4]))
    r.append(float(lastls[5]))
    g.append(float(lastls[6]))
    b.append(float(lastls[7]))

    x = np.array(x)
    r = np.array(r)
    g = np.array(g)
    b = np.array(b)

    if colorModel == "HSV":
        for i in range(r.shape[0]):
            # convert HSV to RGB
            rr, gg, bb = colorsys.hsv_to_rgb(r[i]/360., g[i], b[i])
            r[i] = rr
            g[i] = gg
            b[i] = bb
    elif colorModel == "RGB":
        r /= 255.
        g /= 255.
        b /= 255.

    red = []
    blue = []
    green = []
    xNorm = (x - x[0])/(x[-1] - x[0])
    for i in range(len(x)):
        red.append([xNorm[i], r[i], r[i]])
        green.append([xNorm[i], g[i], g[i]])
        blue.append([xNorm[i], b[i], b[i]])

    # return colormap
    cdict = dict(red=red, green=green, blue=blue)
    return mcolors.LinearSegmentedColormap(name=name, segmentdata=cdict)


def load_cpt_colormap(cptfile, name=None):
    """Read a color map from a cpt file

    Parameters
    ----------
    cptfile : str or open file-like object
        path to .cpt file
    name : str, optional
        name for color map
        if not provided, the file name will be used
    """
    with open(cptfile, 'r', encoding='utf-8') as cptf:
        return gmtColormap_openfile(cptf, name=name)


def shape2grid(x, y, shapefile):
    """
    Convert shapefile to a grid mask.

    Parameters:
    -----------
    x : array-like
        X-coordinates of the grid.
    y : array-like
        Y-coordinates of the grid.
    shapefile : str
        Path to the shapefile.

    Returns:
    --------
    numpy.ndarray
        Grid mask where 1 indicates inside the shape, NaN outside.
    """
    # Load the shapefile
    gdf = gpd.read_file(shapefile)
    # Create a meshgrid
    X, Y = np.meshgrid(x-360, y)

    # Check if each point is inside the shapefile
    points = [Point(x, y) for x, y in zip(X.ravel(), Y.ravel())]
    inside = [gdf.contains(point).any() for point in points]

    # Reshape the result to the shape of the meshgrid
    inside = np.array(inside).reshape(X.shape)
    inside = np.where(inside, 1, np.nan)
    return inside
