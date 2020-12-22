# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2020 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2020 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2020 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
"""
Plotting utilities, including routines for plotting:
    * time series (and comparing with reference data sets)
    * remapped horizontal fields (and comparing with reference data sets)
    * vertical sections on native grid
    * NINO34 time series and spectra
"""
# Authors
# -------
# Xylar Asay-Davis, Milena Veneziani, Luke Van Roekel, Greg Streletz

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import matplotlib.pyplot as plt

from mpas_analysis.shared.plot.title import limit_title


def plot_1D(config, xArrays, fieldArrays, errArrays,
            lineColors=None, lineStyles=None, markers=None, lineWidths=None,
            legendText=None, title=None, xlabel=None, ylabel=None,
            fileout='plot_1D.png',
            figsize=(10, 4), dpi=None,
            xLim=None,
            yLim=None,
            invertYAxis=False,
            maxTitleLength=80):  # {{{
    """
    Plots a 1D line plot with error bars if available.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    xArrays : list of float arrays
        x array (latitude, or any other x axis except time)

    fieldArrays : list of float arrays
        y array (any field as function of x)

    errArrays : list of float arrays
        error array (y errors)

    lineColors, lineStyles, markers, legendText : list of str, optional
        control line color, style, marker, and corresponding legend
        text.  Default is black, solid line with no marker, and no legend.

    lineWidths : list of float, optional
        control line width.  Default is 1.0.

    title : str, optional
        title of plot

    xlabel, ylabel : str, optional
        label of x- and y-axis

    fileout : str, optional
        the file name to be written

    figsize : tuple of float, optional
        size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure, taken from section ``plot``
        option ``dpi`` in the config file by default

    xLim : float array, optional
        x range of plot

    yLim : float array, optional
        y range of plot

    invertYAxis : logical, optional
        if True, invert Y axis

    maxTitleLength : int, optional
        the maximum number of characters in the title and legend, beyond which
        they are truncated with a trailing ellipsis
    """
    # Authors
    # -------
    # Mark Petersen, Milena Veneziani

    # set up figure
    if dpi is None:
        dpi = config.getint('plot', 'dpi')
    plt.figure(figsize=figsize, dpi=dpi)

    plotLegend = False
    for dsIndex in range(len(xArrays)):
        xArray = xArrays[dsIndex]
        fieldArray = fieldArrays[dsIndex]
        errArray = errArrays[dsIndex]
        if xArray is None:
            continue

        if legendText is None:
            label = None
        else:
            label = limit_title(legendText[dsIndex], maxTitleLength)
            plotLegend = True
        if lineColors is None:
            color = 'k'
        else:
            color = lineColors[dsIndex]
        if markers is None:
            marker = None
        else:
            marker = markers[dsIndex]
        if lineStyles is None:
            linestyle = '-'
        else:
            linestyle = lineStyles[dsIndex]
        if lineWidths is None:
            linewidth = 1.
        else:
            linewidth = lineWidths[dsIndex]

        plt.plot(xArray, fieldArray, color=color, linestyle=linestyle,
                 marker=marker, linewidth=linewidth, label=label)
        if errArray is not None:
            plt.fill_between(xArray, fieldArray, fieldArray + errArray,
                             facecolor=color, alpha=0.2)
            plt.fill_between(xArray, fieldArray, fieldArray - errArray,
                             facecolor=color, alpha=0.2)
    plt.grid()
    plt.axhline(0.0, linestyle='-', color='k')  # horizontal lines
    if plotLegend and len(xArrays) > 1:
        plt.legend()

    axis_font = {'size': config.get('plot', 'axisFontSize')}
    title_font = {'size': config.get('plot', 'titleFontSize'),
                  'color': config.get('plot', 'titleFontColor'),
                  'weight': config.get('plot', 'titleFontWeight')}
    if title is not None:
        title = limit_title(title, max_title_length=maxTitleLength)
        plt.title(title, **title_font)
    if xlabel is not None:
        plt.xlabel(xlabel, **axis_font)
    if ylabel is not None:
        plt.ylabel(ylabel, **axis_font)

    if invertYAxis:
        plt.gca().invert_yaxis()

    if xLim:
        plt.xlim(xLim)
    if yLim:
        plt.ylim(yLim)

    if fileout is not None:
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    plt.close()  # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
