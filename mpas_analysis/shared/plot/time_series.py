# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/main/LICENSE
"""
Functions for plotting time series (and comparing with reference data sets)
"""
# Authors
# -------
# Xylar Asay-Davis, Milena Veneziani, Luke Van Roekel, Greg Streletz

import matplotlib
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np

from mpas_analysis.shared.timekeeping.utility import date_to_days

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.plot.ticks import plot_xtick_format
from mpas_analysis.shared.plot.title import limit_title


def timeseries_analysis_plot(config, dsvalues, calendar, title, xlabel, ylabel,
                             movingAveragePoints=None, lineColors=None,
                             lineStyles=None, markers=None, lineWidths=None,
                             legendText=None, maxPoints=None,
                             titleFontSize=None, defaultFontSize=None,
                             figsize=(12, 6), dpi=None,
                             firstYearXTicks=None, yearStrideXTicks=None,
                             maxXTicks=20, obsMean=None, obsUncertainty=None,
                             obsLegend=None, legendLocation='lower left',
                             maxTitleLength=None):
    """
    Plots the list of time series data sets.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    dsvalues : list of xarray DataSets
        the data set(s) to be plotted

    title : str
        the title of the plot

    xlabel, ylabel : str
        axis labels

    calendar : str
        the calendar to use for formatting the time axis

    movingAveragePoints : int, optional
        the number of time points over which to perform a moving average

    lineColors, lineStyles, markers, legendText : list of str, optional
        control line color, style, marker, and corresponding legend
        text.  Default is black, solid line with no marker, and no legend.

    lineWidths : list of float, optional
        control line width.  Default is 1.0.

    maxPoints : list of {None, int}, optional
        the approximate maximum number of time points to use in a time series.
        This can be helpful for reducing the number of symbols plotted if
        plotting with markers.  Otherwise the markers become indistinguishable
        from each other.

    titleFontSize : int, optional
        the size of the title font

    defaultFontSize : int, optional
        the size of text other than the title

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure, taken from section ``plot``
        option ``dpi`` in the config file by default

    firstYearXTicks : int, optional
        The year of the first tick on the x axis.  By default, the first time
        entry is the first tick.

    yearStrideXTicks : int, optional
        The number of years between x ticks. By default, the stride is chosen
        automatically to have ``maxXTicks`` tick marks or fewer.

    maxXTicks : int, optional
        the maximum number of tick marks that will be allowed along the x axis.
        This may need to be adjusted depending on the figure size and aspect
        ratio.

    obsMean, obsUncertainty : list of float, optional
        Mean values and uncertainties for observations to be plotted as error
        bars. The two lists must have the same number of elements.

    obsLegend : list of str, optional
        The label in the legend for each element in ``obsMean`` (and
        ``obsUncertainty``)

    legendLocation : str, optional
        The location of the legend (see ``pyplot.legend()`` for details)

    maxTitleLength : int or None, optional
        the maximum number of characters in the title, beyond which it is
        truncated with a trailing ellipsis.  The default is from the
        ``maxTitleLength`` config option.

    Returns
    -------
    fig : ``matplotlib.figure.Figure``
        The resulting figure
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani, Stephen Price

    if maxTitleLength is None:
        maxTitleLength = config.getint('plot', 'maxTitleLength')

    if defaultFontSize is None:
        defaultFontSize = config.getint('plot', 'defaultFontSize')
    matplotlib.rc('font', size=defaultFontSize)

    if dpi is None:
        dpi = config.getint('plot', 'dpi')
    fig = plt.figure(figsize=figsize, dpi=dpi)

    minDays = []
    maxDays = []
    labelCount = 0
    for dsIndex in range(len(dsvalues)):
        dsvalue = dsvalues[dsIndex]
        if dsvalue is None:
            continue
        if movingAveragePoints == 1 or movingAveragePoints is None:
            mean = dsvalue
        else:
            mean = pd.Series.rolling(dsvalue.to_pandas(), movingAveragePoints,
                                     center=True).mean()
            mean = xr.DataArray.from_series(mean)
        minDays.append(mean.Time.min())
        maxDays.append(mean.Time.max())

        if maxPoints is not None and maxPoints[dsIndex] is not None:
            nTime = mean.sizes['Time']
            if maxPoints[dsIndex] < nTime:
                stride = int(round(nTime / float(maxPoints[dsIndex])))
                mean = mean.isel(Time=slice(0, None, stride))

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

        plt.plot(mean['Time'].values, mean.values, color=color,
                 linestyle=linestyle, marker=marker, linewidth=linewidth,
                 label=label)

    if obsMean is not None:
        obsCount = len(obsMean)
        assert(len(obsUncertainty) == obsCount)

        # space the observations along the time line, leaving gaps at either
        # end
        start = np.amin(minDays)
        end = np.amax(maxDays)
        obsTimes = np.linspace(start, end, obsCount + 2)[1:-1]
        obsSymbols = ['o', '^', 's', 'D', '*']
        obsColors = [config.get('timeSeries', 'obsColor{}'.format(index+1))
                     for index in range(5)]
        for iObs in range(obsCount):
            if obsMean[iObs] is not None:
                symbol = obsSymbols[np.mod(iObs, len(obsSymbols))]
                color = obsColors[np.mod(iObs, len(obsColors))]
                plt.errorbar(obsTimes[iObs],
                             obsMean[iObs],
                             yerr=obsUncertainty[iObs],
                             fmt=symbol,
                             color=color,
                             ecolor=color,
                             capsize=0,
                             label=obsLegend[iObs])
                # plot a box around the error bar to make it more visible
                boxHalfWidth = 0.01 * (end - start)
                boxHalfHeight = obsUncertainty[iObs]
                boxX = obsTimes[iObs] + \
                    boxHalfWidth * np.array([-1, 1, 1, -1, -1])
                boxY = obsMean[iObs] + \
                    boxHalfHeight * np.array([-1, -1, 1, 1, -1])

                plt.plot(boxX, boxY, '-', color=color, linewidth=3)
                labelCount += 1

    if labelCount > 1:
        plt.legend(loc=legendLocation)

    ax = plt.gca()

    if titleFontSize is None:
        titleFontSize = config.get('plot', 'titleFontSize')
    axis_font = {'size': config.get('plot', 'axisFontSize')}
    title_font = {'size': titleFontSize,
                  'color': config.get('plot', 'titleFontColor'),
                  'weight': config.get('plot', 'titleFontWeight')}

    if firstYearXTicks is not None:
        minDays = date_to_days(year=firstYearXTicks, calendar=calendar)

    plot_xtick_format(calendar, minDays, maxDays, maxXTicks,
                      yearStride=yearStrideXTicks)

    # Add a y=0 line if y ranges between positive and negative values
    yaxLimits = ax.get_ylim()
    if yaxLimits[0] * yaxLimits[1] < 0:
        x = ax.get_xlim()
        plt.plot(x, np.zeros(np.size(x)), 'k-', linewidth=1.2, zorder=1)

    if title is not None:
        title = limit_title(title, maxTitleLength)
        plt.title(title, **title_font)
    if xlabel is not None:
        plt.xlabel(xlabel, **axis_font)
    if ylabel is not None:
        plt.ylabel(ylabel, **axis_font)

    return fig


def timeseries_analysis_plot_polar(config, dsvalues, title,
                                   movingAveragePoints=None, lineColors=None,
                                   lineStyles=None, markers=None,
                                   lineWidths=None, legendText=None,
                                   titleFontSize=None, defaultFontSize=None,
                                   figsize=(15, 6), dpi=None,
                                   maxTitleLength=None):
    """
    Plots the list of time series data sets on a polar plot.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    dsvalues : list of xarray DataSets
        the data set(s) to be plotted

    movingAveragePoints : int
        the numer of time points over which to perform a moving average

    title : str
        the title of the plot

    lineColors, lineStyles, markers, legendText : list of str, optional
        control line color, style, marker, and corresponding legend
        text.  Default is black, solid line with no marker, and no legend.

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

    maxTitleLength : int or None, optional
        the maximum number of characters in the title, beyond which it is
        truncated with a trailing ellipsis.  The default is from the
        ``maxTitleLength`` config option.

    Returns
    -------
    fig : ``matplotlib.figure.Figure``
        The resulting figure
    """
    # Authors
    # -------
    # Adrian K. Turner, Xylar Asay-Davis

    if maxTitleLength is None:
        maxTitleLength = config.getint('plot', 'maxTitleLength')

    if defaultFontSize is None:
        defaultFontSize = config.getint('plot', 'defaultFontSize')
    matplotlib.rc('font', size=defaultFontSize)

    if dpi is None:
        dpi = config.getint('plot', 'dpi')
    fig = plt.figure(figsize=figsize, dpi=dpi)

    minDays = []
    maxDays = []
    labelCount = 0
    for dsIndex in range(len(dsvalues)):
        dsvalue = dsvalues[dsIndex]
        if dsvalue is None:
            continue
        mean = pd.Series.rolling(dsvalue.to_pandas(), movingAveragePoints,
                                 center=True).mean()
        mean = xr.DataArray.from_series(mean)
        minDays.append(mean.Time.min())
        maxDays.append(mean.Time.max())

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

        plt.polar((mean['Time'] / 365.0) * np.pi * 2.0, mean, color=color,
                  linestyle=linestyle, marker=marker, linewidth=linewidth,
                  label=label)

    if labelCount > 1:
        plt.legend(loc='lower left')

    ax = plt.gca()

    # set azimuthal axis formatting
    majorTickLocs = np.zeros(12)
    minorTickLocs = np.zeros(12)
    majorTickLocs[0] = 0.0
    minorTickLocs[0] = (constants.daysInMonth[0] * np.pi) / 365.0
    for month in range(1, 12):
        majorTickLocs[month] = majorTickLocs[month - 1] + \
            ((constants.daysInMonth[month - 1] * np.pi * 2.0) / 365.0)
        minorTickLocs[month] = minorTickLocs[month - 1] + \
            (((constants.daysInMonth[month - 1] +
               constants.daysInMonth[month]) * np.pi) / 365.0)

    ax.set_xticks(majorTickLocs)
    ax.set_xticklabels([])

    ax.set_xticks(minorTickLocs, minor=True)
    ax.set_xticklabels(constants.abrevMonthNames, minor=True)

    if titleFontSize is None:
        title = limit_title(title, maxTitleLength)
        titleFontSize = config.get('plot', 'titleFontSize')

    title_font = {'size': titleFontSize,
                  'color': config.get('plot', 'titleFontColor'),
                  'weight': config.get('plot', 'titleFontWeight')}
    if title is not None:
        plt.title(title, **title_font)

    return fig
