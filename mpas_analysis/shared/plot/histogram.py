# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
"""
Functions for plotting histograms (and comparing with reference data sets)
"""
# Authors
# -------
# Carolyn Begeman, Adrian Turner, Xylar Asay-Davis

import matplotlib
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np

from mpas_analysis.shared.timekeeping.utility import date_to_days

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.plot.ticks import plot_xtick_format
from mpas_analysis.shared.plot.title import limit_title

def histogram_analysis_plot(config, dsvalues, calendar, title, xlabel, ylabel,
                            bins=20, range=None, density=True, weights=None, lineColors=None,
                            lineStyles=None, markers=None, lineWidths=None,
                            legendText=None,
                            titleFontSize=None, axisFontSize=None, defaultFontSize=None,
                            figsize=(12, 6), dpi=None,
                            legendLocation='upper right',
                            maxTitleLength=90):

    """
    Plots the list of histogram data sets.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    dsvalues : list of xarray DataSets
        the data set(s) to be plotted. For area or volume quantities, multiply by the
        respective area or volume before inputting. Datasets should already be sliced
        within the time range specified in the config file.

    title : str
        the title of the plot

    xlabel, ylabel : str
        axis labels

    calendar : str
        the calendar to use for formatting the time axis

    density : logical
        if True, normalize the histogram so that the area under the curve is 1

    weights: list of numpy data arrays or NoneType's of length dsvalues
        the weights corresponding to each entry in dsvalues

    lineColors, lineStyles, legendText : list of str, optional
        control line color, style, and corresponding legend
        text.  Default is black, solid line, and no legend.

    lineWidths : list of float, optional
        control line width.  Default is 1.0.

    titleFontSize : int, optional
        the size of the title font

    defaultFontSize : int, optional
        the size of text other than the title

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure, taken from section ``plot``
        option ``dpi`` in the config file by default

    legendLocation : str, optional
        The location of the legend (see ``pyplot.legend()`` for details)

    maxTitleLength : int, optional
        the maximum number of characters in the title and legend, beyond which
        they are truncated with a trailing ellipsis

    Returns
    -------
    fig : ``matplotlib.figure.Figure``
        The resulting figure
    """
    # Authors
    # -------
    # Carolyn Begeman, Adrian Turner, Xylar Asay-Davis

    if defaultFontSize is None:
        defaultFontSize = config.getint('plot', 'defaultFontSize')
    matplotlib.rc('font', size=defaultFontSize)

    if dpi is None:
        dpi = config.getint('plot', 'dpi')

    fig = plt.figure(figsize=figsize, dpi=dpi)
    if title is not None:
        if titleFontSize is None:
            titleFontSize = config.get('plot', 'titleFontSize')
        title_font = {'size': titleFontSize,
                      'color': config.get('plot', 'titleFontColor'),
                      'weight': config.get('plot', 'titleFontWeight')}
    if axisFontSize is None:
        axisFontSize = config.get('plot', 'axisFontSize')
    axis_font = {'size': axisFontSize}

    ax = plt.gca()
    labelCount = 0
    for dsIndex, dsvalue in enumerate(dsvalues):
        if dsvalue is None:
            continue
        if legendText is None:
            label = None
        else:
            label = legendText[dsIndex]
            if label is not None:
                label = limit_title(label, maxTitleLength)
            labelCount += 1
        if lineColors is None:
            color = 'k'
        else:
            color = lineColors[dsIndex]
        if lineStyles is None:
            linestyle = '-'
        else:
            linestyle = lineStyles[dsIndex]
        if markers is None:
            marker = None
        else:
            marker = markers[dsIndex]
        if lineWidths is None:
            linewidth = 1.
        else:
            linewidth = lineWidths[dsIndex]

        hist_values = dsvalue.values.ravel()
        weight = weights[dsIndex]
        hist_type = 'step'
        ax.hist(hist_values, range=range, bins=bins, weights=weight,
                linestyle=linestyle, linewidth=linewidth,
                histtype=hist_type, label=label, density=density)
        if labelCount > 1:
            plt.legend(loc=legendLocation)

    if title is not None:
        title = limit_title(title, maxTitleLength)
        plt.title(title, **title_font)
    if xlabel is not None:
        plt.xlabel(xlabel, **axis_font)
    if ylabel is not None:
        plt.ylabel(ylabel, **axis_font)

    return fig
