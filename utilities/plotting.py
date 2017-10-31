# plotting utilities

import numpy as np
import matplotlib.pyplot as plt
from utilities import stats
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.patches import Polygon as pg
from utilities.util import *



def plot_seismicity_rate(dataframe, fig, ax, **kwargs):
    """
    plots seismicity rate for earthquake catalog
    
    accepts kwargs for matplotlib axes object.
    
    assumes dataframe has timestamp column

    dataframe : pandas.DataFrame
    fig : mpl Figure
    ax : mpl Axes
    kwargs : axes kwargs
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
    
    df : pandas.DataFrame
    fig : mpl Figure
    ax : mpl Axes
    bins : int
    range : list
    kwargs : axes kwargs
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
