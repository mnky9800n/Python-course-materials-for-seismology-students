# -*- encoding: utf-8 -*-
"""
functions useful for statistical seismology.

"""

import numpy as np
import pandas as pd
from utilities import polygon_selection

def mc_maximum_curvature(magnitudes):
    """
    calculates magnitude of completeness using maximum
    curvature method.

    citation: ???

    :param catalog : pandas Series
    :param method : string
    """

    minimum = round(magnitudes.min(), 2)
    bins = np.arange(start=minimum, stop=10, step=0.1)
    hist, edges = np.histogram(a=magnitudes, range=(minimum, 10), bins=bins)
    hist_maximum_index = np.argmax(hist)

    return round(edges[hist_maximum_index], 2)


def fmd_values(magnitudes, bin_width=0.1):
    """
    returns a,b,bstd, n-values

    :param magnitudes:
    :type magnitudes:
    :param bin_width:
    :type bin_width:
    :return:
    :rtype:
    """


    length = magnitudes.shape[0]
    minimum = magnitudes.min()
    average = magnitudes.mean()
    b_value = (1 / (average - (minimum - (bin_width / 2)))) * np.log10(np.exp(1))

    sigma_mag = np.sum([((m - average) ** 2) / (length * (length - 1)) for m in magnitudes])
    b_error = 2.3 * b_value ** 2 * np.sqrt(sigma_mag)

    a_value = np.log10(length) + b_value * minimum

    return a_value, b_value, b_error, length

def calc_fmd_stats_with_mc(magnitudes):
    """
    calculates fmd statistics (a, b, bstd, n, mc) using maximum curvature
    method to calculate magnitude of completeness

    :param magnitudes: magnitudes array
    :type magnitudes: pandas.Series
    :return: calculated fmd statistics
    :rtype: list
    """
    if len(magnitudes) > 0:
        mc = mc_maximum_curvature(magnitudes) + 0.2
        magnitudes = magnitudes[magnitudes >= mc]
        if len(magnitudes) > 0:
            fmd_stats = fmd_values(magnitudes)
            return fmd_stats + (mc,)
        else: return (np.nan, np.nan, np.nan, np.nan, np.nan)
    else:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)

def get_cumdist(data):
    hist, edges = np.histogram(a=data, bins=100, range=(0,10))
    chist = np.cumsum(hist[::-1])
    return edges, hist, chist

def generate_b_based_on_location_mag_normal_error(df, location, radius):
    """
    Assumes a dataframe with columns:
    
    lon        : longitude
    lat        : latitude
    mag        : magnitude
    hz_err_deg : error in the horizontal plane
    
    assumes a 0.1 standard deviation of magnitude error
    assumes all data comes from a normal distribution    
    """
    err_df = df.copy()
    err_df['lon'] = np.random.normal(df['lon'].values, df['hz_err_deg'].values+0.001)
    err_df['lat'] = np.random.normal(df['lat'].values, df['hz_err_deg'].values+0.001)
    err_df['mag'] = np.random.normal(df['mag'], 0.1)
    node_df = get_node_data(node=location, radius=radius, data=err_df)
    a, b, bstd, n, mc = calc_fmd_stats_with_mc(node_df.mag)
    return a, b, bstd, n, mc

def calc_bootstrapped_fmd_values(df, n_calculations):
    """
    calculates bootstrapped fmd values

    :param df: input dataframe
    :type df: pandas.dataframe
    :param n_calculations: number of times to bootstrap input dataframe
    :type n_calculations: int
    :return: fmd statistics (a,b,bstd,n,mc) calculated
    :rtype: list
    """
    fmd_values = []
    for n in range(n_calculations):
        fmd_values.append(calc_fmd_stats_with_mc(df.ix[np.random.choice(df.index, size=(1, df.shape[0]))[0]].mag))
    return fmd_values

def get_catalog_shifted_by_location_normal_error(df):
    """
    shifts catalog dataframe locations by errors assuming a
    normal distribution

    :param df: catalog dataframe
    :type df: pandas.dataframe
    :return: shifted catalog dataframe
    :rtype: pandas.dataframe
    """
    # TODO: this should be moved to polygon_selection
    err_df = df.copy()
    err_df['hz_err_deg'] = err_df['horizontal_error'] / 111.113
    err_df['lon'] = np.random.normal(err_df['lon'].values, err_df['hz_err_deg'].values+0.001)
    err_df['lat'] = np.random.normal(err_df['lat'].values, err_df['hz_err_deg'].values+0.001)
    return err_df


def calculate_b_value_parameter_sweep(dataframe, location, n_iterations, parameters):
    """
    calculates grid search data for fmd statistics
    """
    # TODO: this is probably too specific for here
    rows = []
    for r, t in parameters:
        raw_df = dataframe.loc[dataframe.index >= t].copy()
        raw_df = polygon_selection.get_node_data(data=raw_df, radius=r, node=location, m=1)

        try:
            b = calc_bootstrapped_fmd_values(get_catalog_shifted_by_location_normal_error(raw_df)
                                                   , n_iterations)
            bdf = pd.DataFrame(np.array(b), columns=['a', 'b', 'bstd', 'n', 'mc'])
            row = [(r,), (t,), bdf.mean().values, bdf.std().values]
            row = np.concatenate(np.array(row))

        except ValueError:
            row = (r,) + (t,) + (np.nan,) * 10

        rows.append(row)

    bdf = pd.DataFrame(rows, columns=['radius', 'start_time', 'a_avg', 'b_avg', 'bstd_avg', 'n_avg', 'mc_avg'
        , 'a_std', 'b_std', 'bstd_std', 'n_std', 'mc_std'])
    bdf[['radius', 'a_avg', 'b_avg', 'bstd_avg', 'n_avg', 'mc_avg'
        , 'a_std', 'b_std', 'bstd_std', 'n_std', 'mc_std']] = bdf[['radius', 'a_avg'
        , 'b_avg', 'bstd_avg'
        , 'n_avg', 'mc_avg'
        , 'a_std', 'b_std'
        , 'bstd_std', 'n_std'
        , 'mc_std']].apply(pd.to_numeric)
    bdf['start_time'] = pd.to_datetime(bdf['start_time'])

    return bdf