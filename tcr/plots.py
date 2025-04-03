"""
This module provides a collection of functions for advanced plotting and visualization 
in PyTCR, including map formatting, density plotting, track visualization, and more.
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
    Configures the map's appearance by defining its spatial extent, grid spacing, 
    and optionally overlaying a shapefile.

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The Matplotlib axis object on which the map will be rendered.
    extent : tuple, optional
        Defines the spatial bounds of the map in data coordinates 
        (left, right, bottom, top). Default is (-180, 180, -90, 90).
    shapefile : str or shapefile-like object, optional
        Path to a shapefile or a shapefile-like object to overlay on the map.
    gridlabel : bool, optional
        Enables or disables the display of gridline labels. Default is False.
    projection : cartopy.crs.Projection, optional
        The map projection to use. Default is None, which implies a 
        PlateCarree projection.
    coastlines : bool, optional
        Determines whether coastlines are displayed. Default is True.
    add_features : bool, optional
        Specifies whether additional geographical features (e.g., land, borders) 
        are displayed. Default is True.
    gridspace : int, optional
        Sets the spacing between gridlines in degrees. Default is 10.

    Notes:
    ------
    - If a shapefile is provided, its geometries will be overlaid on the map.
    - Gridline labels and spacing can be customized for better visualization.
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

    Notes:
    ------
    - The function supports both contour and pcolormesh methods for plotting density.
    - If `logscale` is True, the density values are plotted on a logarithmic scale.
    - The `extent` parameter defines the spatial bounds of the map.

    Example:
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import cartopy.crs as ccrs
    >>> fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    >>> lat = np.linspace(-30, 30, 100)
    >>> lon = np.linspace(120, 180, 100)
    >>> lon, lat = np.meshgrid(lon, lat)
    >>> density = np.exp(-((lat - 0)**2 + (lon - 150)**2) / (2 * 10**2))
    >>> levels = np.linspace(0, 1, 11)
    >>> plot_density(ax, lat, lon, density, levels, extent=(120, 180, -30, 30),
    ...              cmap='viridis', title="TC Density")
    >>> plt.show()
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

    Notes:
    ------
    - Tracks with wind speeds below the `wind_speed_threshold` will not be plotted.
    - The function uses `LineCollection` to efficiently plot multiple line segments.
    - If `wind_color` is True, the tracks will be colored based on wind speed using the provided colormap.

    Example:
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import cartopy.crs as ccrs
    >>> from matplotlib.colors import Normalize
    >>> fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    >>> lats = np.random.uniform(-30, 30, (5, 100))
    >>> lons = np.random.uniform(120, 180, (5, 100))
    >>> vmaxs = np.random.uniform(10, 50, (5, 100))
    >>> track_inds = [0, 1, 2, 3, 4]
    >>> plot_tracks(ax, lats, lons, vmaxs, track_inds, wind_color=True, 
    ...             cmap='viridis', norm=Normalize(10, 50), title="TC Tracks")
    >>> plt.show()
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
    Generate a plot illustrating the exceedance probability of the provided dataset.

    Parameters:
    -----------
    data : array-like
        The dataset for which the exceedance probability plot will be generated.
    ax : matplotlib.axes.Axes
        The Matplotlib Axes object on which the plot will be rendered.
    ylabel : str, optional
        Label for the y-axis. Default is 'Value'.
    xlabel : str, optional
        Label for the x-axis. Default is 'Exceedance Probability'.
    title : str, optional
        Title of the plot. If provided, it will be displayed above the plot. Default is None.
    fontweight : str, optional
        Font weight for the plot title. Default is 'regular'.
    fontsize : int, optional
        Font size for the plot title. Default is 12.
    color : str, optional
        Color of the data points in the plot. Default is 'tab:blue'.

    Notes:
    ------
    - The exceedance probability is calculated as the rank of each data point divided by the total number of points.
    - The x-axis is displayed on a logarithmic scale to emphasize the lower probability range.
    - The function ensures that the data is converted to a NumPy array for consistent processing.
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

    Notes:
    ------
    - The function applies a buffer to the input polygon if a non-zero buffer distance is provided.
    - Longitude values are normalized to the range [-180, 180] to ensure compatibility with the
      polygon's coordinates.
    - Tracks are considered to make landfall if any point along the track lies within the buffered
      polygon.

    Example:
    --------
    >>> from shapely.geometry import Polygon
    >>> import numpy as np
    >>> lat_tracks = np.random.uniform(-90, 90, (100, 50))
    >>> lon_tracks = np.random.uniform(-180, 180, (100, 50))
    >>> polygon = Polygon([(-10, 30), (-10, 40), (10, 40), (10, 30)])
    >>> landfall_indices = get_tracks_landfall_region(lat_tracks, lon_tracks, polygon, buffer=0.5)
    >>> print(landfall_indices)
    [0, 5, 12, 45]
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
        The original colormap to be truncated.
    minval : float, optional
        The lower bound of the colormap range, normalized to [0, 1]. Default is 0.0.
    maxval : float, optional
        The upper bound of the colormap range, normalized to [0, 1]. Default is 1.0.
    n : int, optional
        The number of discrete colors in the truncated colormap. Default is 100.

    Returns:
    --------
    matplotlib.colors.LinearSegmentedColormap
        A new colormap truncated to the specified range.

    Notes:
    ------
    - The `minval` and `maxval` parameters must be within the range [0, 1].
    - The function creates a new colormap by sampling the original colormap
      within the specified range.

    Example:
    --------
    >>> import matplotlib.pyplot as plt
    >>> from matplotlib import cm
    >>> truncated_cmap = truncate_colormap(cm.viridis, minval=0.2, maxval=0.8)
    >>> plt.imshow([[0, 1]], cmap=truncated_cmap, aspect='auto')
    >>> plt.colorbar()
    >>> plt.show()
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
    Generate a buffer zone around a specified Point of Interest (POI).

    This function creates a circular buffer zone around a given latitude and
    longitude, representing the area within a specified radius.

    Parameters:
    -----------
    lat_poi : float
        Latitude of the Point of Interest (in degrees). Must be within the range [-90, 90].
    lon_poi : float
        Longitude of the Point of Interest (in degrees). Must be within the range [-180, 180].
    radius : float, optional
        Radius of the buffer zone (in degrees). Must be a positive value. Default is 1.0.

    Returns:
    --------
    shapely.geometry.Polygon
        A Polygon object representing the buffer zone around the POI.

    Raises:
    -------
    ValueError
        If the radius is not positive or if the longitude is outside the valid range.

    Notes:
    ------
    - The buffer zone is created using a simple Euclidean distance approximation
      and may not account for distortions in certain map projections.
    - If the longitude exceeds 180 degrees, it is normalized to the range [-180, 180].

    Examples:
    ---------
    >>> from shapely.geometry import Polygon
    >>> buffer = create_buffer_around_POI(37.7749, -122.4194, radius=0.5)
    >>> isinstance(buffer, Polygon)
    True
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
    """
    Reads a GMT color map from an open .cpt file and generates a Matplotlib colormap.

    Parameters
    ----------
    cptf : file-like object
        An open file or URL handle pointing to a .cpt file containing the colormap data.
    name : str, optional
        A custom name for the generated colormap. If not provided, the name will be 
        derived from the file name.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
        A Matplotlib colormap object created from the .cpt file.

    Notes
    -----
    - The .cpt file should follow the GMT colormap format.
    - Supports both RGB and HSV color models, automatically converting HSV to RGB.
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
    """
    Load a colormap from a GMT .cpt file.

    Parameters
    ----------
    cptfile : str or file-like object
        Path to the .cpt file or an open file-like object containing the colormap data.
    name : str, optional
        Custom name for the colormap. If not provided, the colormap name will be derived
        from the file name.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
        A colormap object that can be used with Matplotlib plotting functions.

    Notes
    -----
    - The .cpt file should follow the GMT colormap format.
    - Supports both RGB and HSV color models, automatically converting HSV to RGB.
    - If the file is provided as a path, it will be opened and read internally.
    """
    with open(cptfile, 'r', encoding='utf-8') as cptf:
        return gmtColormap_openfile(cptf, name=name)


def shape2grid(x, y, shapefile):
    """
    Generates a grid mask from a shapefile, marking grid points that fall
    inside the shape with 1 and points outside with NaN.

    Parameters:
    -----------
    x : array-like
        1D array of X-coordinates representing the grid's longitude values.
    y : array-like
        1D array of Y-coordinates representing the grid's latitude values.
    shapefile : str
        File path to the shapefile containing the geometric shape(s) to be
        used for masking.

    Returns:
    --------
    numpy.ndarray
        A 2D array (grid mask) with the same shape as the meshgrid created 
        from `x` and `y`. Grid cells inside the shape are marked with 1, 
        while cells outside are assigned NaN.

    Notes
    -----
    - The function assumes the shapefile is in a compatible coordinate 
      reference system (CRS) with the provided `x` and `y` coordinates.
    - The X-coordinates are adjusted by subtracting 360 to handle 
      longitude wrapping, if necessary.
    - The function uses the `shapely` library to determine whether 
      grid points are inside the shape(s) defined in the shapefile.

    Examples
    >>> import numpy as np
    >>> x = np.linspace(-180, 180, 360)
    >>> y = np.linspace(-90, 90, 180)
    >>> shapefile = "path/to/shapefile.shp"
    >>> grid_mask = shape2grid(x, y, shapefile)
    >>> print(grid_mask.shape)
    (180, 360)
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
