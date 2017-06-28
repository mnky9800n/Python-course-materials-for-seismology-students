# plotting utilities

import numpy as np
import matplotlib.pyplot as plt
from utilities import stats
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.patches import Polygon as pg
from utilities.util import *


def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.
 
    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")
 
    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)
 
    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):
 
        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu
 
    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
 
    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
 
    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi
 
    return (glon2, glat2, baz)
 
def equi(m, centerlon, centerlat, radius, *args, **kwargs):
    # TODO : pick a better name
    """
    plots circle on matplotlib basemap map

    :param m:
    :type m:
    :param centerlon:
    :type centerlon:
    :param centerlat:
    :type centerlat:
    :param radius:
    :type radius:
    :param args:
    :type args:
    :param kwargs:
    :type kwargs:
    :return:
    :rtype:
    """
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])
 
    X,Y = m(X,Y)
    plt.plot(X,Y,**kwargs)

def plot_circle_on_map(m, centerlon, centerlat, radius, **kwargs):
    """
    Wrapper for equi

    :param m:
    :type m:
    :param centerlon:
    :type centerlon:
    :param centerlat:
    :type centerlat:
    :param radius:
    :type radius:
    :param kwargs:
    :type kwargs:
    :return:
    :rtype:
    """
    return equi(m, centerlon, centerlat, radius, **kwargs)

def plot_line_on_map(m, point_1, point_2, s=5, l=5, color='red'):
    """
    plots line on matplotlib basemap map

    :param m:
    :type m:
    :param point_1:
    :type point_1:
    :param point_2:
    :type point_2:
    :param s:
    :type s:
    :param l:
    :type l:
    :param color:
    :type color:
    :return:
    :rtype:
    """
    x = (point_1[0], point_2[0])
    y = (point_1[1], point_2[1])
    x, y = m(x, y)
    m.plot(x, y, marker='D', color=color, markersize=s, linewidth=l)
    
def plot_text_on_map(m, ax, lat, lon, text, fontsize=15):
    """
    plots text on matplotlib basemap map

    :param m:
    :type m:
    :param ax:
    :type ax:
    :param lat:
    :type lat:
    :param lon:
    :type lon:
    :param text:
    :type text:
    :param fontsize:
    :type fontsize:
    :return:
    :rtype:
    """
    x, y = m(lon, lat)
    ax.text(s=text, x=x, y=y, fontsize=fontsize)

def draw_screen_poly( lats, lons, m, ax):
    """
    draws polygon on matplotlib basemap map

    :param lats:
    :type lats:
    :param lons:
    :type lons:
    :param m:
    :type m:
    :param ax:
    :type ax:
    :return:
    :rtype:
    """
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = pg( xy, facecolor='None', edgecolor='red', linestyle='--', linewidth=5 )
    ax.add_patch(poly)
    
def plot_seismicity_map(dataframe, lon_lat_min_max=None, **kwargs):
    """
    Plots seismicity map for provided catalog
    
    Plots depth as color, magnitude proportional
    to size. Plots down to 100km. >100km depths
    are reduced to 100km in color chart.
    
    Fails on catalogs with >1 million events
    (perhaps smaller).
    
    
    """
    df = dataframe.copy()
    if lon_lat_min_max is None:
        lat_min = np.floor(df.lat.min())
        lat_max = np.ceil(df.lat.max())
        lon_min = np.floor(df.lon.min())
        lon_max = np.ceil(df.lon.max())
    else:
        lon_min, lon_max, lat_min, lat_max = lon_lat_min_max

    fig, ax = plt.subplots(1, figsize=(8,8))
    m = Basemap(projection='merc'
           ,llcrnrlat=lat_min
           ,urcrnrlat=lat_max
           ,llcrnrlon=lon_min
           ,urcrnrlon=lon_max
           ,resolution='i'
           ,area_thresh=1000
           ,ax=ax)

    lat_labels = np.arange(lat_min, lat_max, int((lat_max - lat_min) / 4) + 1)
    lon_labels = np.arange(lon_min, lon_max, int((lat_max - lat_min) / 4) + 1)

    m.drawparallels(lat_labels, labels=lat_labels)
    m.drawmeridians(lon_labels, labels=lon_labels)
    
    m.drawcoastlines()
    m.fillcontinents(color='0.72', zorder=0)
    
    x, y = m(df.lon.values, df.lat.values)
    # TODO : make the color user changeable
    cbar = ax.scatter(x, y, c=df.depth.values, s=1*np.exp(df.mag.values/2.), edgecolor='None'
                  , cmap='rainbow', alpha=0.5, vmin=0, vmax=100)
    c1 = fig.colorbar(cbar, label='depth (km)',fraction=0.0346, pad=0.084)
    c1.ax.invert_yaxis()
    
    return m, fig, ax

def plot_seismicity_rate(dataframe, fig, ax, **kwargs):
    """
    plots seismicity rate for earthquake catalog
    
    accepts kwargs for matplotlib axes object.
    
    assumes dataframe has timestamp column
    """
    df = dataframe.copy()
#     df = df.set_index('timestamp')
    df['seismicity_rate'] = np.arange(0, df.shape[0], 1)
    
    # fig, ax = plt.subplots(figsize=(10, 4))
    df.seismicity_rate.plot(ax=ax, **kwargs)
    ax.set_xlabel('')
    ax.set_ylabel('cumulative eq count')
    if df.shape[0] > 10000:
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        
    return fig, ax
    
def plot_radius_time_sweep(dataframe, vertical_axis, value, fig, ax, colorbar=True, **kwargs):

    """
    expects a dataframe with the following columns
    
    start_time : numpy.datetime64
    radius : float64 or int
    
    ::params::
    ----------------------------------
    dataframe : pandas.dataframe
    value : dataframe column header to be plotted as color map
    kwargs : any values that can be used with matplotlib.pyplot.pcolormesh
    
    ::dependencies::
    ----------------------------------
    matplotlib
    numpy
    pandas
    
    ::notes::
    ----------------------------------
    this could probably, and should probably, be generalized for any axis
    case. or at least such that the vertical axis is specified by the user
    """

    # TODO: provide option to for time axis to be "years before" instead of explicit date
    zi = dataframe.pivot(index='start_time', columns=vertical_axis, values=value)
    xi_label = [np.datetime64(z, 'Y').astype(str) for z in zi.index]
    xi = np.arange(len(xi_label))
    yi = zi.columns
    xi, yi = np.meshgrid(xi, yi)
    zi = np.ma.masked_invalid(zi).transpose()
    cbar = ax.pcolormesh(xi, yi, zi, **kwargs)
    if colorbar is True:
        fig.colorbar(cbar, label=str(value))
    
    xi_ticks = [i for i, j in  enumerate(replace_unique_items(xi_label)) if j is not None]
    xi_labels = [j for i, j in  enumerate(replace_unique_items(xi_label)) if j is not None]
    
    ax.set_xticks(xi_ticks)
    ax.set_xticklabels(xi_labels, rotation=90)
    
    ax.set_ylabel(vertical_axis)
    
    return fig, ax


def plot_fmd_diagram(df, fig, ax, bins=100, range=[0, 10], **kwargs):
    """
    Plots fmd diagram with fit line for given magnitudes.
    
    :param df: 
    :param fig: 
    :param ax: 
    :param bins: 
    :param range: 
    :param kwargs: 
    :return: 
    """
    mags = df.mag
    hist, edges = np.histogram(a=mags, bins=bins, range=range)
    chist = np.cumsum(hist[::-1])[::-1]

    a, b, bstd, n, mc = stats.calc_fmd_stats_with_mc(mags)

    ax.scatter(edges[:-1], chist, marker='s', color='black', **kwargs)

    x = np.arange(range[0], range[1], 0.1)
    y = 10 ** (a - b * x)

    ax.plot(x, y, color='red')

    ax.set_yscale('log')
    ax.set_ylim(1e0, 10 ** a)
    ax.set_xlim(0, mags.max() + 1)
    ax.set_title('b={b}$\pm${bstd}, n={n}, mc={mc}'.format(b=round(b, 2), bstd=round(bstd, 2), n=n, mc=mc))
