# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2018 Los Alamos National Security, LLC. All rights reserved.
# Copyright (c) 2018 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2018 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
"""
Plotting utilities, including routines for plotting:
    * time series (and comparing with control data sets)
    * remapped horizontal fields (and comparing with control data sets)
    * vertical sections on native grid
    * NINO34 time series and spectra
"""
# Authors
# -------
# Xylar Asay-Davis, Milena Veneziani, Luke Van Roekel, Greg Streletz

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import xarray as xr
import pandas as pd
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import addcyclic
from matplotlib.ticker import FuncFormatter, FixedLocator
import numpy as np
from functools import partial
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import xml.etree.ElementTree as ET
from six.moves import configparser
import cmocean
import pkg_resources

from mpas_analysis.shared.timekeeping.utility import days_to_datetime, \
    date_to_days

from mpas_analysis.shared.constants import constants


def timeseries_analysis_plot(config, dsvalues, N, title, xlabel, ylabel,
                             fileout, calendar, lineColors=None,
                             lineStyles=None, markers=None, lineWidths=None,
                             legendText=None, maxPoints=None,
                             titleFontSize=None, figsize=(15, 6), dpi=None,
                             firstYearXTicks=None, yearStrideXTicks=None,
                             maxXTicks=20, obsMean=None, obsUncertainty=None,
                             obsLegend=None, legendLocation='lower left'):

    """
    Plots the list of time series data sets and stores the result in an image
    file.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    dsvalues : list of xarray DataSets
        the data set(s) to be plotted

    N : int
        the numer of time points over which to perform a moving average

    title : str
        the title of the plot

    xlabel, ylabel : str
        axis labels

    fileout : str
        the file name to be written

    calendar : str
        the calendar to use for formatting the time axis

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
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani, Stephen Price

    if dpi is None:
        dpi = config.getint('plot', 'dpi')
    plt.figure(figsize=figsize, dpi=dpi)

    minDays = []
    maxDays = []
    labelCount = 0
    for dsIndex in range(len(dsvalues)):
        dsvalue = dsvalues[dsIndex]
        if dsvalue is None:
            continue
        if N == 1 or N is None:
            mean = dsvalue
        else:
            mean = pd.Series.rolling(dsvalue.to_pandas(), N,
                                     center=True).mean()
            mean = xr.DataArray.from_series(mean)
        minDays.append(mean.Time.min())
        maxDays.append(mean.Time.max())

        if maxPoints is not None and maxPoints[dsIndex] is not None:
            nTime = mean.sizes['Time']
            if maxPoints[dsIndex] < nTime:
                stride = int(round(nTime/float(maxPoints[dsIndex])))
                mean = mean.isel(Time=slice(0, None, stride))

        if legendText is None:
            label = None
        else:
            label = legendText[dsIndex]
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
        obsTimes = np.linspace(start, end, obsCount+2)[1:-1]
        obsSymbols = ['o', '^', 's', 'D', '*']
        obsColors = ['b', 'g', 'c', 'm', 'r']
        for iObs in range(obsCount):
            if obsMean[iObs] is not None:
                symbol = obsSymbols[np.mod(iObs, len(obsSymbols))]
                color = obsColors[np.mod(iObs, len(obsColors))]
                fmt = '{}{}'.format(color, symbol)
                plt.errorbar(obsTimes[iObs],
                             obsMean[iObs],
                             yerr=obsUncertainty[iObs],
                             fmt=fmt,
                             ecolor=color,
                             capsize=0,
                             label=obsLegend[iObs])
                # plot a box around the error bar to make it more visible
                boxHalfWidth = 0.01*(end - start)
                boxHalfHeight = obsUncertainty[iObs]
                boxX = obsTimes[iObs] + \
                    boxHalfWidth*np.array([-1, 1, 1, -1, -1])
                boxY = obsMean[iObs] + \
                    boxHalfHeight*np.array([-1, -1, 1, 1, -1])

                plt.plot(boxX, boxY, '{}-'.format(color), linewidth=3)
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
    if yaxLimits[0]*yaxLimits[1] < 0:
        x = ax.get_xlim()
        plt.plot(x, np.zeros(np.size(x)), 'k-', linewidth=1.2, zorder=1)

    if title is not None:
        plt.title(title, **title_font)
    if xlabel is not None:
        plt.xlabel(xlabel, **axis_font)
    if ylabel is not None:
        plt.ylabel(ylabel, **axis_font)
    if fileout is not None:
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    plt.close()


def timeseries_analysis_plot_polar(config, dsvalues, N, title,
                                   fileout, lineColors=None, lineStyles=None,
                                   markers=None, lineWidths=None,
                                   legendText=None, titleFontSize=None,
                                   figsize=(15, 6), dpi=None):

    """
    Plots the list of time series data sets on a polar plot and stores
    the result in an image file.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    dsvalues : list of xarray DataSets
        the data set(s) to be plotted

    N : int
        the numer of time points over which to perform a moving average

    title : str
        the title of the plot

    fileout : str
        the file name to be written

    lineColors, lineStyles, markers, legendText : list of str, optional
        control line color, style, marker, and corresponding legend
        text.  Default is black, solid line with no marker, and no legend.

    lineWidths : list of float, optional
        control line width.  Default is 1.0.

    titleFontSize : int, optional
        the size of the title font

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure, taken from section ``plot``
        option ``dpi`` in the config file by default
    """
    # Authors
    # -------
    # Adrian K. Turner, Xylar Asay-Davis

    if dpi is None:
        dpi = config.getint('plot', 'dpi')
    plt.figure(figsize=figsize, dpi=dpi)

    minDays = []
    maxDays = []
    labelCount = 0
    for dsIndex in range(len(dsvalues)):
        dsvalue = dsvalues[dsIndex]
        if dsvalue is None:
            continue
        mean = pd.Series.rolling(dsvalue.to_pandas(), N, center=True).mean()
        mean = xr.DataArray.from_series(mean)
        minDays.append(mean.Time.min())
        maxDays.append(mean.Time.max())

        if legendText is None:
            label = None
        else:
            label = legendText[dsIndex]
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

        plt.polar((mean['Time']/365.0)*np.pi*2.0, mean, color=color,
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
        majorTickLocs[month] = majorTickLocs[month-1] + \
            ((constants.daysInMonth[month-1] * np.pi * 2.0) / 365.0)
        minorTickLocs[month] = minorTickLocs[month-1] + \
            (((constants.daysInMonth[month-1] +
               constants.daysInMonth[month]) * np.pi) / 365.0)

    ax.set_xticks(majorTickLocs)
    ax.set_xticklabels([])

    ax.set_xticks(minorTickLocs, minor=True)
    ax.set_xticklabels(constants.abrevMonthNames, minor=True)

    if titleFontSize is None:
        titleFontSize = config.get('plot', 'titleFontSize')

    title_font = {'size': titleFontSize,
                  'color': config.get('plot', 'titleFontColor'),
                  'weight': config.get('plot', 'titleFontWeight')}
    if title is not None:
        plt.title(title, **title_font)

    if fileout is not None:
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    plt.close()


def plot_polar_comparison(
        config,
        Lons,
        Lats,
        modelArray,
        controlArray,
        diffArray,
        colorMapSectionName,
        fileout,
        title=None,
        plotProjection='npstere',
        latmin=50.0,
        lon0=0,
        modelTitle='Model',
        controlTitle='Observations',
        diffTitle='Model-Observations',
        cbarlabel='units',
        titleFontSize=None,
        figsize=None,
        dpi=None,
        vertical=False):

    """
    Plots a data set around either the north or south pole.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    Lons, Lats : float arrays
        longitude and latitude arrays

    modelArray, controlArray : float arrays
        model and observational or reference run data sets

    diffArray : float array
        difference between modelArray and controlArray

    colorMapSectionName : str
        section name in ``config`` where color map info can be found.

    fileout : str
        the file name to be written

    title : str, optional
        the subtitle of the plot

    plotProjection : str, optional
        Basemap projection for the plot

    modelTitle : str, optional
        title of the model panel

    controlTitle : str, optional
        title of the observations or reference run panel

    diffTitle : str, optional
        title of the difference (bias) panel

    cbarlabel : str, optional
        label on the colorbar

    titleFontSize : int, optional
        size of the title font

    figsize : tuple of float, optional
        the size of the figure in inches.  If ``None``, the figure size is
        ``(8, 22)`` if ``vertical == True`` and ``(22, 8)`` otherwise.

    dpi : int, optional
        the number of dots per inch of the figure, taken from section ``plot``
        option ``dpi`` in the config file by default

    vertical : bool, optional
        whether the subplots should be stacked vertically rather than
        horizontally
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani

    def do_subplot(ax, field, title, colormap, norm, levels, ticks, contours,
                   lineWidth, lineColor):
        """
        Make a subplot within the figure.
        """

        m = Basemap(projection=plotProjection, boundinglat=latmin,
                    lon_0=lon0, resolution='l', ax=ax)

        fieldPeriodic, LatsPeriodic, LonsPeriodic = addcyclic(field, Lats,
                                                              Lons)

        x, y = m(LonsPeriodic, LatsPeriodic)  # compute map proj coordinates

        ax.set_title(title, y=1.06, **plottitle_font)

        m.drawcoastlines()
        m.fillcontinents(color='grey', lake_color='white')
        m.drawparallels(np.arange(-80., 81., 10.))
        m.drawmeridians(np.arange(-180., 181., 20.),
                        labels=[True, False, True, True])

        if levels is None:
            plotHandle = m.pcolormesh(x, y, fieldPeriodic, cmap=colormap,
                                      norm=norm)
        else:
            plotHandle = m.contourf(x, y, fieldPeriodic, cmap=colormap,
                                    norm=norm, levels=levels)

        if contours is not None:
            matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
            m.contour(x, y, field, levels=contours, colors=lineColor,
                      linewidths=lineWidth)

        cbar = m.colorbar(plotHandle, location='right', pad="3%",
                          spacing='uniform', ticks=ticks, boundaries=levels)
        cbar.set_label(cbarlabel)

    if dpi is None:
        dpi = config.getint('plot', 'dpi')

    dictModelRef = setup_colormap(config, colorMapSectionName, suffix='Result')
    dictDiff = setup_colormap(config, colorMapSectionName, suffix='Difference')

    if controlArray is None:
        if figsize is None:
            figsize = (8, 8.5)
        subplots = [111]
    elif vertical:
        if figsize is None:
            figsize = (8, 22)
        subplots = [311, 312, 313]
    else:
        if figsize is None:
            figsize = (22, 8.5)
        subplots = [131, 132, 133]

    fig = plt.figure(figsize=figsize, dpi=dpi)

    if (title is not None):
        if titleFontSize is None:
            titleFontSize = config.get('plot', 'titleFontSize')
        title_font = {'size': titleFontSize,
                      'color': config.get('plot', 'titleFontColor'),
                      'weight': config.get('plot', 'titleFontWeight')}
        fig.suptitle(title, y=0.95, **title_font)

    plottitle_font = {'size': config.get('plot',
                                         'threePanelPlotTitleFontSize')}

    ax = plt.subplot(subplots[0])
    do_subplot(ax=ax, field=modelArray, title=modelTitle, **dictModelRef)

    if controlArray is not None:
        ax = plt.subplot(subplots[1])
        do_subplot(ax=ax, field=controlArray, title=controlTitle,
                   **dictModelRef)

        ax = plt.subplot(subplots[2])
        do_subplot(ax=ax, field=diffArray, title=diffTitle, **dictDiff)

    plt.tight_layout(pad=4.)
    if vertical:
        plt.subplots_adjust(top=0.9)

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    plt.close()


def plot_global_comparison(
        config,
        Lons,
        Lats,
        modelArray,
        controlArray,
        diffArray,
        colorMapSectionName,
        fileout,
        title=None,
        modelTitle='Model',
        controlTitle='Observations',
        diffTitle='Model-Observations',
        cbarlabel='units',
        titleFontSize=None,
        figsize=None,
        dpi=None,
        lineWidth=1,
        lineColor='black'):

    """
    Plots a data set as a longitude/latitude map.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    Lons, Lats : float arrays
        longitude and latitude arrays

    modelArray, controlArray : float arrays
        model and observational or reference run data sets

    diffArray : float array
        difference between modelArray and controlArray

    colorMapSectionName : str
        section name in ``config`` where color map info can be found.

    fileout : str
        the file name to be written

    title : str, optional
        the subtitle of the plot

    modelTitle : str, optional
        title of the model panel

    controlTitle : str, optional
        title of the observations or reference run panel

    diffTitle : str, optional
        title of the difference (bias) panel

    cbarlabel : str, optional
        label on the colorbar

    titleFontSize : int, optional
        size of the title font

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure, taken from section ``plot``
        option ``dpi`` in the config file by default

    lineWidth : int, optional
        the line width of contour lines (if specified)

    lineColor : str, optional
        the color of contour lines (if specified)
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani

    def plot_panel(title, array, colormap, norm, levels, ticks, contours,
                   lineWidth, lineColor):

        plt.title(title, y=1.06, **plottitle_font)

        m.drawcoastlines()
        m.fillcontinents(color='grey', lake_color='white')
        m.drawparallels(np.arange(-80., 80., 20.),
                        labels=[True, False, False, False])
        m.drawmeridians(np.arange(-180., 180., 60.),
                        labels=[False, False, False, True])

        if levels is None:
            plotHandle = m.pcolormesh(x, y, array, cmap=colormap, norm=norm)
        else:
            plotHandle = m.contourf(x, y, array, cmap=colormap, norm=norm,
                                    levels=levels, extend='both')

        if contours is not None:
            matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
            m.contour(x, y, array, levels=contours, colors=lineColor,
                      linewidths=lineWidth)

        cbar = m.colorbar(plotHandle, location='right', pad="5%",
                          spacing='uniform', ticks=ticks, boundaries=ticks)
        cbar.set_label(cbarlabel)

    # set up figure
    if dpi is None:
        dpi = config.getint('plot', 'dpi')
    if figsize is None:
        # set the defaults, depending on if we have 1 or 3 panels
        if controlArray is None:
            figsize = (8, 5)
        else:
            figsize = (8, 13)
    fig = plt.figure(figsize=figsize, dpi=dpi)
    if (title is not None):
        if titleFontSize is None:
            titleFontSize = config.get('plot', 'titleFontSize')
        title_font = {'size': titleFontSize,
                      'color': config.get('plot', 'titleFontColor'),
                      'weight': config.get('plot', 'titleFontWeight')}
        fig.suptitle(title, y=0.95, **title_font)

    plottitle_font = {'size': config.get('plot',
                                         'threePanelPlotTitleFontSize')}

    m = Basemap(projection='cyl', llcrnrlat=-85, urcrnrlat=86, llcrnrlon=-180,
                urcrnrlon=181, resolution='l')
    x, y = m(Lons, Lats)  # compute map proj coordinates

    dictModelRef = setup_colormap(config, colorMapSectionName, suffix='Result')
    dictDiff = setup_colormap(config, colorMapSectionName, suffix='Difference')

    if controlArray is not None:
        plt.subplot(3, 1, 1)

    plot_panel(modelTitle, modelArray, **dictModelRef)

    if controlArray is not None:
        plt.subplot(3, 1, 2)
        plot_panel(controlTitle, controlArray, **dictModelRef)

        plt.subplot(3, 1, 3)
        plot_panel(diffTitle, diffArray, **dictDiff)

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    plt.close()


def plot_polar_projection_comparison(
        config,
        x,
        y,
        invalidMask,
        modelArray,
        controlArray,
        diffArray,
        fileout,
        colorMapSectionName,
        title=None,
        modelTitle='Model',
        controlTitle='Observations',
        diffTitle='Model-Observations',
        cbarlabel='units',
        titleFontSize=None,
        figsize=None,
        dpi=None,
        lineWidth=0.5,
        lineColor='black',
        vertical=False):

    """
    Plots a data set as a longitude/latitude map.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    x, y : numpy ndarrays
        1D x and y arrays defining the projection grid

    ivalidMask : numpy ndarray
        mask of where there is no data (e.g. land, for ocean or sea-ice output)

    modelArray, controlArray : numpy ndarrays
        model and observational or reference run data sets

    diffArray : float array
        difference between modelArray and controlArray

    fileout : str
        the file name to be written

    colorMapSectionName : str
        section name in ``config`` where color map info can be found.

    title : str, optional
        the subtitle of the plot

    plotProjection : str, optional
        Basemap projection for the plot

    modelTitle : str, optional
        title of the model panel

    controlTitle : str, optional
        title of the observations or reference run panel

    diffTitle : str, optional
        title of the difference (bias) panel

    cbarlabel : str, optional
        label on the colorbar

    titleFontSize : int, optional
        size of the title font

    figsize : tuple of float, optional
        the size of the figure in inches.  If ``None``, the figure size is
        ``(8, 22)`` if ``vertical == True`` and ``(22, 8)`` otherwise.

    dpi : int, optional
        the number of dots per inch of the figure, taken from section ``plot``
        option ``dpi`` in the config file by default

    lineWidth : int, optional
        the line width of contour lines (if specified)

    lineColor : str, optional
        the color of contour lines (if specified)

    vertical : bool, optional
        whether the subplots should be stacked vertically rather than
        horizontally
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def plot_panel(ax, title, array, colormap, norm, levels, ticks, contours,
                   lineWidth, lineColor):

        plt.title(title, y=1.06, **plottitle_font)

        if levels is None:
            plotHandle = plt.pcolormesh(x, y, array, cmap=colormap, norm=norm)
        else:
            plotHandle = plt.contourf(xCenter, yCenter, array, cmap=colormap,
                                      norm=norm, levels=levels, extend='both')

        plt.pcolormesh(x, y, invalidMask, cmap=landColorMap)
        plt.contour(xCenter, yCenter, invalidMask.mask, (0.5,), colors='k',
                    linewidths=0.5)

        if contours is not None:
            matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
            plt.contour(x, y, array, levels=contours, colors=lineColor,
                        linewidths=lineWidth)

        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        cbar = plt.colorbar(plotHandle, cax=cax)
        cbar.set_label(cbarlabel)
        if ticks is not None:
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(['{}'.format(tick) for tick in ticks])

        ax.axis('off')
        ax.set_aspect('equal')
        ax.autoscale(tight=True)

    # set up figure
    if dpi is None:
        dpi = config.getint('plot', 'dpi')

    if controlArray is None:
        if figsize is None:
            figsize = (8, 7.5)
        subplots = [111]
    elif vertical:
        if figsize is None:
            figsize = (8, 22)
        subplots = [311, 312, 313]
    else:
        if figsize is None:
            figsize = (22, 7.5)
        subplots = [131, 132, 133]

    dictModelRef = setup_colormap(config, colorMapSectionName, suffix='Result')
    dictDiff = setup_colormap(config, colorMapSectionName, suffix='Difference')

    fig = plt.figure(figsize=figsize, dpi=dpi)

    if (title is not None):
        if titleFontSize is None:
            titleFontSize = config.get('plot', 'titleFontSize')
        title_font = {'size': titleFontSize,
                      'color': config.get('plot', 'titleFontColor'),
                      'weight': config.get('plot', 'titleFontWeight')}
        fig.suptitle(title, y=0.95, **title_font)

    plottitle_font = {'size': config.get('plot',
                                         'threePanelPlotTitleFontSize')}

    # set up land colormap
    colorList = [(0.8, 0.8, 0.8), (0.8, 0.8, 0.8)]
    landColorMap = cols.LinearSegmentedColormap.from_list('land', colorList)

    # locations of centers for contour plots
    xCenter = 0.5*(x[1:] + x[0:-1])
    yCenter = 0.5*(y[1:] + y[0:-1])

    ax = plt.subplot(subplots[0])
    plot_panel(ax, modelTitle, modelArray, **dictModelRef)

    if controlArray is not None:
        ax = plt.subplot(subplots[1])
        plot_panel(ax, controlTitle, controlArray, **dictModelRef)

        ax = plt.subplot(subplots[2])
        plot_panel(ax, diffTitle, diffArray, **dictDiff)

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    plt.close()


def plot_vertical_section_comparison(
        config,
        xArray,
        depthArray,
        modelArray,
        controlArray,
        diffArray,
        fileout,
        colorMapSectionName,
        cbarLabel=None,
        xlabel=None,
        ylabel=None,
        title=None,
        modelTitle='Model',
        controlTitle='Observations',
        diffTitle='Model-Observations',
        titleFontSize=None,
        plotTitleFontSize=None,
        axisFontSize=None,
        figsize=None,
        dpi=None,
        lineWidth=2,
        lineStyle='solid',
        lineColor='black',
        backgroundColor='grey',
        xLim=None,
        yLim=None,
        secondXAxisData=None,
        secondXAxisLabel=None,
        thirdXAxisData=None,
        thirdXAxisLabel=None,
        numUpperTicks=None,
        upperXAxisTickLabelPrecision=None,
        invertYAxis=True,
        xArrayIsTime=False,
        N=None,
        firstYearXTicks=None,
        yearStrideXTicks=None,
        maxXTicks=20,
        calendar='gregorian',
        compareAsContours=False,
        contourLineStyle=None,
        comparisonContourLineStyle=None,
        comparisonContourLineColor=None,
        labelContours=False,
        contourLabelPrecision=1):

    """
    Plots vertical section plots in a three-panel format, comparing model data
    (in modelArray) to some control dataset (in controlArray), which can be
    either observations or an alternative model, and also presenting the
    difference plot of the two.  If controlArray is None, then only one panel
    is plotted, displaying the model data.

    If compareAsContours is true, the contours of modelArray and controlArray
    are plotted on a single plot.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    xArray : float array
        x array (latitude, longitude, spherical distance, or distance along
        a transect;  or, time for Hovmoller plots)

    depthArray : float array
        depth array [m]

    modelArray, controlArray : float arrays
        model and observational or reference run data sets

    diffArray : float array
        difference between modelArray and controlArray

    fileout : str
        the file name to be written

    colorMapSectionName : str
        section name in ``config`` where color map info can be found.

    cbarLabel : str, optional
        the label for the colorbar.  If compareAsContours and labelContours are
        both True, colorbarLabel is used as follows (typically in order to
        indicate the units that are associated with the contour labels):
        if controlArray is None, the colorbarLabel string is parenthetically
        appended to the plot title;  if controlArray is not None, it is
        parenthetically appended to the legend entries of the contour
        comparison plot.

    xlabel, ylabel : str, optional
        label of x- and y-axis

    title : str, optional
        the subtitle of the plot

    modelTitle : str, optional
        title of the model panel

    controlTitle : str, optional
        title of the observations or reference run panel

    diffTitle : str, optional
        title of the difference (bias) panel

    titleFontSize : int, optional
        size of the title font

    plotTitleFontSize : int, optional
        size of the title font for the individual plots

    axisFontSize : int, optional
        size of the axis font

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure, taken from section ``plot``
        option ``dpi`` in the config file by default

    lineWidth : int, optional
        the line width of contour lines (if specified)

    lineStyle : str, optional
        the line style of contour lines (if specified); this applies to the
        contour lines on heatmaps and to the contour lines of the model field
        on contour comparison plots (the line style of the contour lines of
        the control field on contour comparison plots is set using the
        contourComparisonLineStyle argument).

    lineColor : str, optional
        the color of contour lines (if specified); this applies to the
        contour lines on heatmaps and to the contour lines of the model field
        on contour comparison plots (the line color of the contour lines of
        the control field on contour comparison plots is set using the
        contourComparisonLineColor argument).

    backgroundColor : str, optional
        the background color for the plot (NaNs and masked areas will be
        shown in this color)

    xLim : float array, optional
        x range of plot

    yLim : float array, optional
        y range of plot

    secondXAxisData : the data to use to display a second x axis (which will be
        placed above the plot).  This array must have the same number of values
        as xArray, and it is assumed that the values in this array define
        locations along the x axis that are the same as those defined by the
        corresponding values in xArray, but in some different unit system.

    secondXAxisLabel : the label for the second x axis, if requested

    thirdXAxisData : the data to use to display a third x axis (which will be
        placed above the plot and above the second x axis, which must be
        specified if a third x axis is to be specified).  This array must have
        the same number of values as xArray, and it is assumed that the values
        in this array define locations along the x axis that are the same as
        those defined by the corresponding values in xArray, but in some
        different unit system (which is presumably also different from the unit
        system used for the values in the secondXAxisData array).  The typical
        use for this third axis is for transects, for which the primary x axis
        represents distance along a transect, and the second and third x axes
        are used to display the corresponding latitudes and longitudes.

    thirdXAxisLabel : the label for the third x axis, if requested

    numUpperTicks : the approximate number of ticks to use on the upper x axis
        or axes (these are the second and third x axes, which are placed above
        the plot if they have been requested by specifying the secondXAxisData
        or thirdXAxisData arrays above)

    upperXAxisTickLabelPrecision : the number of decimal places (to the right
        of the decimal point) to use for values at upper axis ticks.  This
        value can be adjusted (in concert with numUpperTicks) to avoid problems
        with overlapping numbers along the upper axis.

    invertYAxis : logical, optional
        if True, invert Y axis

    xArrayIsTime : logical, optional
        if True, format the x axis for time (this applies only to the primary
        x axis, not to the optional second or third x axes)

    N : int, optional
        the number of points over which to perform a moving average
        NOTE: this option is mostly intended for use when xArrayIsTime is True,
        although it will work with other data as well.  Also, the moving
        average calculation is based on number of points, not actual x axis
        values, so for best results, the values in the xArray should be equally
        spaced.

    firstYearXTicks : int, optional
        The year of the first tick on the x axis.  By default, the first time
        entry is the first tick.

    yearStrideXTicks : int, optional
        The number of years between x ticks. By default, the stride is chosen
        automatically to have ``maxXTicks`` tick marks or fewer.

    maxXTicks : int, optional
        the maximum number of tick marks that will be allowed along the primary
        x axis.  This may need to be adjusted depending on the figure size and
        aspect ratio.  NOTE:  maxXTicks is only used if xArrayIsTime is True

    calendar : str, optional
        the calendar to use for formatting the time axis
        NOTE:  calendar is only used if xArrayIsTime is True

    compareAsContours : bool, optional
       if compareAsContours is True, instead of creating a three panel plot
       showing modelArray, controlArray, and their difference, the function
       will plot the contours of modelArray and controlArray on a single plot
       (unless controlArray is None, in which case only the contours of
       modelArray will be plotted on the single panel plot).

    comparisonContourLineStyle : str, optional
        the line style of contour lines of the control field on a contour
        comparison plot

    comparisonContourLineColor : str, optional
        the line color of contour lines of the control field on a contour
        comparison plot

    labelContours : bool, optional
        whether or not to label contour lines (if specified) with their values

    contourLabelPrecision : int, optional
        the precision (in terms of number of figures to the right of the
        decimal point) of contour labels

    """
    # Authors
    # -------
    # Greg Streletz, Xylar Asay-Davis, Milena Veneziani

    if controlArray is None or compareAsContours:
        singlePanel = True
    else:
        singlePanel = False

    # set up figure
    if dpi is None:
        dpi = config.getint('plot', 'dpi')
    if figsize is None:
        # set the defaults, depending on if we have 1 or 3 panels, and
        # depending on how many x axes are to be displayed on the plots
        if singlePanel:
            if compareAsContours and controlArray is not None:
                if thirdXAxisData is not None:
                    figsize = (8, 8)
                else:
                    figsize = (8, 7)
            else:
                figsize = (8, 5)
        elif thirdXAxisData is not None:
            figsize = (8, 17)
        else:
            figsize = (8, 13)

    fig = plt.figure(figsize=figsize, dpi=dpi)

    if (title is not None):
        if titleFontSize is None:
            titleFontSize = config.get('plot', 'threePanelTitleFontSize')
        title_font = {'size': titleFontSize,
                      'color': config.get('plot', 'threePanelTitleFontColor'),
                      'weight': config.get('plot',
                                           'threePanelTitleFontWeight')}
        fig.suptitle(title, y=0.95, **title_font)

    if plotTitleFontSize is None:
        plotTitleFontSize = config.get('plot', 'threePanelPlotTitleFontSize')

    if thirdXAxisData is not None:
        if singlePanel:
            titleY = 1.64
        else:
            titleY = 1.34
    elif secondXAxisData is not None:
        titleY = 1.20
    else:
        titleY = 1.06

    if axisFontSize is None:
        axisFontSize = config.get('plot', 'threePanelAxisFontSize')

    if not singlePanel:
        plt.subplot(3, 1, 1)

    if not compareAsContours or controlArray is None:
        title = modelTitle
        contourComparisonFieldArray = None
        comparisonFieldName = None
        originalFieldName = None
    else:
        title = None
        contourComparisonFieldArray = controlArray
        comparisonFieldName = controlTitle
        originalFieldName = modelTitle

    plot_vertical_section(
            config,
            xArray,
            depthArray,
            modelArray,
            colorMapSectionName,
            suffix='Result',
            colorbarLabel=cbarLabel,
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
            fileout=None,
            titleFontSize=plotTitleFontSize,
            titleY=titleY,
            axisFontSize=axisFontSize,
            xLim=xLim,
            yLim=yLim,
            lineWidth=lineWidth,
            lineStyle=lineStyle,
            lineColor=lineColor,
            secondXAxisData=secondXAxisData,
            secondXAxisLabel=secondXAxisLabel,
            thirdXAxisData=thirdXAxisData,
            thirdXAxisLabel=thirdXAxisLabel,
            numUpperTicks=numUpperTicks,
            upperXAxisTickLabelPrecision=upperXAxisTickLabelPrecision,
            invertYAxis=invertYAxis,
            xArrayIsTime=xArrayIsTime,
            N=None,
            firstYearXTicks=firstYearXTicks,
            yearStrideXTicks=yearStrideXTicks,
            maxXTicks=maxXTicks, calendar=calendar,
            backgroundColor=backgroundColor,
            plotAsContours=compareAsContours,
            contourComparisonFieldArray=contourComparisonFieldArray,
            comparisonFieldName=comparisonFieldName,
            originalFieldName=originalFieldName,
            comparisonContourLineStyle=comparisonContourLineStyle,
            comparisonContourLineColor=comparisonContourLineColor,
            labelContours=labelContours,
            contourLabelPrecision=contourLabelPrecision)

    if not singlePanel:
        plt.subplot(3, 1, 2)
        plot_vertical_section(
                config,
                xArray,
                depthArray,
                controlArray,
                colorMapSectionName,
                suffix='Result',
                colorbarLabel=cbarLabel,
                title=controlTitle,
                xlabel=xlabel,
                ylabel=ylabel,
                fileout=None,
                titleFontSize=plotTitleFontSize,
                titleY=titleY,
                axisFontSize=axisFontSize,
                xLim=xLim,
                yLim=yLim,
                lineWidth=lineWidth,
                lineStyle=lineStyle,
                lineColor=lineColor,
                secondXAxisData=secondXAxisData,
                secondXAxisLabel=secondXAxisLabel,
                thirdXAxisData=thirdXAxisData,
                thirdXAxisLabel=thirdXAxisLabel,
                upperXAxisTickLabelPrecision=upperXAxisTickLabelPrecision,
                numUpperTicks=numUpperTicks,
                invertYAxis=invertYAxis,
                xArrayIsTime=xArrayIsTime,
                N=None,
                firstYearXTicks=firstYearXTicks,
                yearStrideXTicks=yearStrideXTicks,
                maxXTicks=maxXTicks,
                calendar=calendar,
                backgroundColor=backgroundColor,
                labelContours=labelContours,
                contourLabelPrecision=contourLabelPrecision)

        plt.subplot(3, 1, 3)
        plot_vertical_section(
                config,
                xArray,
                depthArray,
                diffArray,
                colorMapSectionName,
                suffix='Difference',
                colorbarLabel=cbarLabel,
                title=diffTitle,
                xlabel=xlabel,
                ylabel=ylabel,
                fileout=None,
                titleFontSize=plotTitleFontSize,
                titleY=titleY,
                axisFontSize=axisFontSize,
                xLim=xLim,
                yLim=yLim,
                lineWidth=lineWidth,
                lineStyle=lineStyle,
                lineColor=lineColor,
                secondXAxisData=secondXAxisData,
                secondXAxisLabel=secondXAxisLabel,
                thirdXAxisData=thirdXAxisData,
                thirdXAxisLabel=thirdXAxisLabel,
                upperXAxisTickLabelPrecision=upperXAxisTickLabelPrecision,
                numUpperTicks=numUpperTicks,
                invertYAxis=invertYAxis,
                xArrayIsTime=xArrayIsTime,
                N=None,
                firstYearXTicks=firstYearXTicks,
                yearStrideXTicks=yearStrideXTicks,
                maxXTicks=maxXTicks,
                calendar=calendar,
                backgroundColor=backgroundColor,
                labelContours=labelContours,
                contourLabelPrecision=contourLabelPrecision)

    if singlePanel:
        if thirdXAxisData is not None and controlArray is None:
            plt.tight_layout(pad=0.0, h_pad=2.0, rect=[0.0, 0.0, 1.0, 0.98])
        else:
            plt.tight_layout(pad=0.0, h_pad=2.0, rect=[0.0, 0.0, 1.0, 0.80])
    else:
        plt.tight_layout(pad=0.0, h_pad=2.0, rect=[0.0, 0.0, 1.0, 0.88])

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    plt.close()


def plot_1D(config, xArrays, fieldArrays, errArrays,
            lineColors=None, lineStyles=None, markers=None, lineWidths=None,
            legendText=None, title=None, xlabel=None, ylabel=None,
            fileout='plot_1D.png',
            figsize=(10, 4), dpi=None,
            xLim=None,
            yLim=None,
            invertYAxis=False):  # {{{

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
            label = legendText[dsIndex]
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
            plt.fill_between(xArray, fieldArray, fieldArray+errArray,
                             facecolor=color, alpha=0.2)
            plt.fill_between(xArray, fieldArray, fieldArray-errArray,
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

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    plt.close()

    return  # }}}


def plot_vertical_section(
        config,
        xArray,
        depthArray,
        fieldArray,
        colorMapSectionName,
        suffix='',
        colorbarLabel=None,
        title=None,
        xlabel=None,
        ylabel=None,
        fileout=None,
        figsize=(10, 4),
        dpi=None,
        titleFontSize=None,
        titleY=None,
        axisFontSize=None,
        xLim=None,
        yLim=None,
        lineWidth=2,
        lineStyle='solid',
        lineColor='black',
        backgroundColor='grey',
        secondXAxisData=None,
        secondXAxisLabel=None,
        thirdXAxisData=None,
        thirdXAxisLabel=None,
        numUpperTicks=None,
        upperXAxisTickLabelPrecision=None,
        invertYAxis=True,
        xArrayIsTime=False,
        N=None,
        firstYearXTicks=None,
        yearStrideXTicks=None,
        maxXTicks=20,
        calendar='gregorian',
        plotAsContours=False,
        contourComparisonFieldArray=None,
        comparisonFieldName=None,
        originalFieldName=None,
        comparisonContourLineStyle=None,
        comparisonContourLineColor=None,
        labelContours=False,
        contourLabelPrecision=1):  # {{{

    """
    Plots a data set as a x distance (latitude, longitude,
    or spherical distance) vs depth map (vertical section).

    Or, if xArrayIsTime is True, plots data set on a vertical
    Hovmoller plot (depth vs. time).

    Typically, the fieldArray data are plotted using a heatmap, but if
    contourComparisonFieldArray is not None, then contours of both
    fieldArray and contourComparisonFieldArray are plotted instead.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    xArray : float array
        x array (latitude, longitude, or spherical distance; or, time for a
        Hovmoller plot)

    depthArray : float array
        depth array [m]

    fieldArray : float array
        field array to plot

    colorMapSectionName : str
        section name in ``config`` where color map info can be found.

    suffix : str, optional
        the suffix used for colorbar config options

    colorbarLabel : str, optional
        the label for the colorbar.  If plotAsContours and labelContours are
        both True, colorbarLabel is used as follows (typically in order to
        indicate the units that are associated with the contour labels):
        if contourComparisonFieldArray is None, the colorbarLabel string is
        parenthetically appended to the plot title;  if
        contourComparisonFieldArray is not None, it is parenthetically appended
        to the legend entries of the contour comparison plot.

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

    titleFontSize : int, optional
        size of the title font

    titleY : float, optional
        the y value to use for placing the plot title

    axisFontSize : int, optional
        size of the axis font

    xLim : float array, optional
        x range of plot

    yLim : float array, optional
        y range of plot

    lineWidth : int, optional
        the line width of contour lines (if specified)

    lineStyle : str, optional
        the line style of contour lines (if specified); this applies to the
        style of contour lines of fieldArray (the style of the contour lines
        of contourComparisonFieldArray is set using
        contourComparisonLineStyle).

    lineColor : str, optional
        the color of contour lines (if specified); this applies to the
        contour lines of fieldArray (the color of the contour lines of
        contourComparisonFieldArray is set using contourComparisonLineColor

    backgroundColor : str, optional
        the background color for the plot (NaNs will be shown in this color)

    secondXAxisData : the data to use to display a second x axis (which will be
        placed above the plot).  This array must have the same number of values
        as xArray, and it is assumed that the values in this array define
        locations along the x axis that are the same as those defined by the
        corresponding values in xArray, but in some different unit system.

    secondXAxisLabel : the label for the second x axis, if requested

    thirdXAxisData : the data to use to display a third x axis (which will be
        placed above the plot and above the second x axis, which must be
        specified if a third x axis is to be specified).  This array must have
        the same number of values as xArray, and it is assumed that the values
        in this array define locations along the x axis that are the same as
        those defined by the corresponding values in xArray, but in some
        different unit system (which is presumably also different from the unit
        system used for the values in the secondXAxisData array).  The typical
        use for this third axis is for transects, for which the primary x axis
        represents distance along a transect, and the second and third x axes
        are used to display the corresponding latitudes and longitudes.

    thirdXAxisLabel : the label for the third x axis, if requested

    numUpperTicks : the approximate number of ticks to use on the upper x axis
        or axes (these are the second and third x axes, which are placed above
        the plot if they have been requested by specifying the secondXAxisData
        or thirdXAxisData arrays above)

    upperXAxisTickLabelPrecision : the number of decimal places (to the right
        of the decimal point) to use for values at upper axis ticks.  This
        value can be adjusted (in concert with numUpperTicks) to avoid problems
        with overlapping numbers along the upper axis.

    invertYAxis : logical, optional
        if True, invert Y axis

    xArrayIsTime : logical, optional
        if True, format the x axis for time (this applies only to the primary
        x axis, not to the optional second or third x axes)

    N : int, optional
        the number of points over which to perform a moving average
        NOTE: this option is mostly intended for use when xArrayIsTime is True,
        although it will work with other data as well.  Also, the moving
        average calculation is based on number of points, not actual x axis
        values, so for best results, the values in the xArray should be equally
        spaced.

    firstYearXTicks : int, optional
        The year of the first tick on the x axis.  By default, the first time
        entry is the first tick.

    yearStrideXTicks : int, optional
        The number of years between x ticks. By default, the stride is chosen
        automatically to have ``maxXTicks`` tick marks or fewer.

    maxXTicks : int, optional
        the maximum number of tick marks that will be allowed along the primary
        x axis.  This may need to be adjusted depending on the figure size and
        aspect ratio.  NOTE:  maxXTicks is only used if xArrayIsTime is True

    calendar : str, optional
        the calendar to use for formatting the time axis
        NOTE:  calendar is only used if xArrayIsTime is True

    plotAsContours : bool, optional
        if plotAsContours is True, instead of plotting fieldArray as a
        heatmap, the function will plot only the contours of fieldArray.  In
        addition, if contourComparisonFieldArray is not None, the contours
        of this field will be plotted on the same plot.  The selection of
        contour levels is still determined as for the contours on the heatmap
        plots, via the 'contours' entry in colorMapSectionName.

    contourComparisonFieldArray : float array, optional
        a comparison field array (typically observational data or results from
        another simulation run), assumed to be of the same shape as fieldArray,
        and related to xArray and depthArray in the same way fieldArray is.
        If contourComparisonFieldArray is None, then fieldArray will be plotted
        as a heatmap.  However, if countourComparisonFieldArray is not None,
        then contours of both fieldArray and contourComparisonFieldArray will
        be plotted in order to enable a comparison of the two fields on the
        same plot.  If plotAsContours is False, this parameter is ignored.

    comparisonFieldName : str, optional
       the name for the comparison field.  If contourComparisonFieldArray is
       None, this parameter is ignored.

    originalFieldName : str, optional
       the name for the fieldArray field (for the purposes of labeling the
       contours on a contour comparison plot).  If contourComparisonFieldArray
       is None, this parameter is ignored.

    comparisonContourLineStyle : str, optional
        the line style of contour lines of the comparisonFieldName field on
        a contour comparison plot

    comparisonContourLineColor : str, optional
        the line color of contour lines of the comparisonFieldName field on
        a contour comparison plot

    labelContours : bool, optional
        whether or not to label contour lines (if specified) with their values

    contourLabelPrecision : int, optional
        the precision (in terms of number of figures to the right of the
        decimal point) of contour labels

    """
    # Authors
    # -------
    # Milena Veneziani, Mark Petersen, Xylar Asay-Davis, Greg Streletz

    dimX = xArray.shape
    dimZ = depthArray.shape
    dimF = fieldArray.shape
    if contourComparisonFieldArray is not None:
        dimC = contourComparisonFieldArray.shape

    if len(dimX) != 1 and len(dimX) != 2:
        raise ValueError('xArray must have either one or two dimensions '
                         '(has %d)' % dimX)

    if len(dimZ) != 1 and len(dimZ) != 2:
        raise ValueError('depthArray must have either one or two dimensions '
                         '(has %d)' % dimZ)

    if len(dimF) != 2:
        raise ValueError('fieldArray must have two dimensions (has %d)' % dimF)

    if contourComparisonFieldArray is not None:
        if len(dimC) != 2:
            raise ValueError('contourComparisonFieldArray must have two '
                             'dimensions (has %d)' % dimC)
        elif (fieldArray.shape[0] != contourComparisonFieldArray.shape[0]) or \
             (fieldArray.shape[1] != contourComparisonFieldArray.shape[1]):
            raise ValueError('size mismatch between fieldArray (%d x %d) and '
                             'contourComparisonFieldArray (%d x %d)' %
                             (fieldArray.shape[0], fieldArray.shape[1],
                              contourComparisonFieldArray.shape[0],
                              contourComparisonFieldArray.shape[1]))

    # verify that the dimensions of fieldArray are consistent with those of
    # xArray and depthArray
    if len(dimX) == 1 and len(dimZ) == 1:
        num_x = dimX[0]
        num_z = dimZ[0]
        if num_x != fieldArray.shape[1] or num_z != fieldArray.shape[0]:
            raise ValueError('size mismatch between xArray (%d), '
                             'depthArray (%d), and fieldArray (%d x %d)' %
                             (num_x, num_z, fieldArray.shape[0],
                              fieldArray.shape[1]))
    elif len(dimX) == 1:
        num_x = dimX[0]
        num_x_Z = dimZ[1]
        num_z = dimZ[0]
        if num_x != fieldArray.shape[1] or num_z != fieldArray.shape[0] or \
                num_x != num_x_Z:
            raise ValueError('size mismatch between xArray (%d), '
                             'depthArray (%d x %d), and fieldArray (%d x %d)' %
                             (num_x, num_z, num_x_Z,
                              fieldArray.shape[0],
                              fieldArray.shape[1]))
    elif len(dimZ) == 1:
        num_x = dimX[1]
        num_z_X = dimX[0]
        num_z = dimZ[0]
        if num_x != fieldArray.shape[1] or num_z != fieldArray.shape[0] or \
                num_z != num_z_X:
            raise ValueError('size mismatch between xArray (%d x %d), '
                             'depthArray (%d), and fieldArray (%d x %d)' %
                             (num_z_X, num_x, num_z,
                              fieldArray.shape[0],
                              fieldArray.shape[1]))
    else:
        num_x = dimX[1]
        num_z_X = dimX[0]
        num_x_Z = dimZ[1]
        num_z = dimZ[0]
        if num_x != fieldArray.shape[1] or num_z != fieldArray.shape[0] \
                or num_x != num_x_Z or num_z != num_z_X:
            raise ValueError('size mismatch between xArray (%d x %d), '
                             'depthArray (%d x %d), and fieldArray (%d x %d)' %
                             (num_z_X, num_x, num_z, num_x_Z,
                              fieldArray.shape[0],
                              fieldArray.shape[1]))

    # Verify that the upper x-axis parameters are consistent with each other
    # and with xArray
    if (secondXAxisData is None and thirdXAxisData is not None):
        raise ValueError('secondXAxisData cannot be None if thirdXAxisData '
                         'is not None')
    if (secondXAxisData is not None):
        arrayShape = secondXAxisData.shape
        if len(arrayShape) == 1 and arrayShape[0] != num_x:
            raise ValueError('secondXAxisData has %d x values, '
                             'but should have num_x = %d x values' %
                             (arrayShape[0], num_x))
        elif len(arrayShape) == 2 and arrayShape[1] != num_x:
            raise ValueError('secondXAxisData has %d x values, '
                             'but should have num_x = %d x values' %
                             (arrayShape[1], num_x))
        elif len(arrayShape) > 2:
            raise ValueError('secondXAxisData must be a 1D or 2D array, '
                             'but is of dimension %d' %
                             (len(arrayShape)))
    if (thirdXAxisData is not None):
        arrayShape = thirdXAxisData.shape
        if len(arrayShape) == 1 and arrayShape[0] != num_x:
            raise ValueError('thirdXAxisData has %d x values, '
                             'but should have num_x = %d x values' %
                             (arrayShape[0], num_x))
        elif len(arrayShape) == 2 and arrayShape[1] != num_x:
            raise ValueError('thirdXAxisData has %d x values, '
                             'but should have num_x = %d x values' %
                             (arrayShape[1], num_x))
        elif len(arrayShape) > 2:
            raise ValueError('thirdXAxisData must be a 1D or 2D array, '
                             'but is of dimension %d' %
                             (len(arrayShape)))

    # define x and y as the appropriate 2D arrays for plotting
    if len(dimX) == 1 and len(dimZ) == 1:
        x, y = np.meshgrid(xArray, depthArray)
    elif len(dimX) == 1:
        x, y = np.meshgrid(xArray, np.zeros(num_z))
        y = depthArray
    elif len(dimZ) == 1:
        x, y = np.meshgrid(np.zeros(num_x), depthArray)
        x = xArray
    else:
        x = xArray
        y = depthArray

    # set up figure
    if dpi is None:
        dpi = config.getint('plot', 'dpi')
    if fileout is not None:
        plt.figure(figsize=figsize, dpi=dpi)

    # compute moving averages with respect to the x dimension
    if N is not None and N != 1:
        movingAverageDepthSlices = []
        for nVertLevel in range(len(depthArray)):
            depthSlice = fieldArray[[nVertLevel]][0]
            # in case it's not an xarray already
            depthSlice = xr.DataArray(depthSlice)
            mean = pd.Series.rolling(depthSlice.to_series(), N,
                                     center=True).mean()
            mean = xr.DataArray.from_series(mean)
            mean = mean[int(N/2.0):-int(round(N/2.0)-1)]
            movingAverageDepthSlices.append(mean)
        xArray = xArray[int(N/2.0):-int(round(N/2.0)-1)]
        fieldArray = xr.DataArray(movingAverageDepthSlices)

    colormapDict = setup_colormap(config, colorMapSectionName, suffix=suffix)

    if not plotAsContours:    # display a heatmap of fieldArray

        if colormapDict['levels'] is None:
            # interpFieldArray contains the values at centers of grid cells,
            # for pcolormesh plots (using bilinear interpolation)
            interpFieldArray = \
                    0.5 * (0.5*(fieldArray[1:, 1:] + fieldArray[0:-1, 1:]) +
                           0.5*(fieldArray[1:, 0:-1] + fieldArray[0:-1, 0:-1]))

            plotHandle = plt.pcolormesh(x, y, interpFieldArray,
                                        cmap=colormapDict['colormap'],
                                        norm=colormapDict['norm'])
        else:
            plotHandle = plt.contourf(x, y, fieldArray,
                                      cmap=colormapDict['colormap'],
                                      norm=colormapDict['norm'],
                                      levels=colormapDict['levels'],
                                      extend='both')

        cbar = plt.colorbar(plotHandle,
                            orientation='vertical',
                            spacing='uniform',
                            aspect=9,
                            ticks=colormapDict['ticks'],
                            boundaries=colormapDict['ticks'])

        if colorbarLabel is not None:
            cbar.set_label(colorbarLabel)

    else:     # display a white heatmap to get a white background for non-land
        zeroArray = np.ma.where(fieldArray != np.nan, 0.0, fieldArray)
        plotHandle = plt.contourf(x, y, zeroArray, colors='white')

    # set the color for NaN or masked regions, and draw a black
    # outline around them; technically, the contour level used should
    # be 1.0, but the contours don't show up when using 1.0, so 0.999
    # is used instead
    ax = plt.gca()
    ax.set_facecolor(backgroundColor)
    landArray = np.ma.where(fieldArray != np.nan, 1.0, fieldArray)
    landArray = np.ma.masked_where(landArray == np.nan, landArray, copy=True)
    landArray = landArray.filled(0.0)
    plt.contour(x, y, landArray, levels=[0.999], colors='black', linewidths=1)

    # plot contours, if they were requested
    contourLevels = colormapDict['contours']
    if contourLevels is not None:
        if len(contourLevels) == 0:
            # automatic calculation of contour levels
            contourLevels = None
        cs1 = plt.contour(x, y, fieldArray,
                          levels=contourLevels,
                          colors=lineColor,
                          linestyles=lineStyle,
                          linewidths=lineWidth)
        if labelContours:
            fmt_string = "%%1.%df" % int(contourLabelPrecision)
            plt.clabel(cs1, fmt=fmt_string)
        if plotAsContours and contourComparisonFieldArray is not None:
            cs2 = plt.contour(x, y, contourComparisonFieldArray,
                              levels=contourLevels,
                              colors=comparisonContourLineColor,
                              linestyles=comparisonContourLineStyle,
                              linewidths=lineWidth)
            if labelContours:
                plt.clabel(cs2, fmt=fmt_string)

    if plotAsContours and contourComparisonFieldArray is not None:
        h1, _ = cs1.legend_elements()
        h2, _ = cs2.legend_elements()
        if labelContours:
            originalFieldName = originalFieldName + " (" + colorbarLabel + ")"
            comparisonFieldName = comparisonFieldName + " (" + \
                colorbarLabel + ")"
        ax.legend([h1[0], h2[0]], [originalFieldName, comparisonFieldName],
                  loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=1)

    if (title is not None):
        if plotAsContours and labelContours \
           and contourComparisonFieldArray is None:
            title = title + " (" + colorbarLabel + ")"
        if titleFontSize is None:
            titleFontSize = config.get('plot', 'titleFontSize')
        title_font = {'size': titleFontSize,
                      'color': config.get('plot', 'titleFontColor'),
                      'weight': config.get('plot', 'titleFontWeight')}
        if titleY is not None:
            plt.title(title, y=titleY, **title_font)
        else:
            plt.title(title, **title_font)

    if (xlabel is not None) or (ylabel is not None):
        if axisFontSize is None:
            axisFontSize = config.get('plot', 'axisFontSize')
        axis_font = {'size': axisFontSize}

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

    if xArrayIsTime:
        if firstYearXTicks is None:
            minDays = [xArray[0]]
        else:
            minDays = date_to_days(year=firstYearXTicks, calendar=calendar)
        maxDays = [xArray[-1]]

        plot_xtick_format(calendar, minDays, maxDays, maxXTicks,
                          yearStride=yearStrideXTicks)

    # add a second x-axis scale, if it was requested
    if secondXAxisData is not None:
        ax2 = ax.twiny()
        ax2.set_facecolor(backgroundColor)
        ax2.set_xlabel(secondXAxisLabel, **axis_font)
        ax2.set_xlim(ax.get_xlim())
        stride = int(round(float(num_x) / (float(numUpperTicks) - 1.0))) - 1
        if stride <= 0:
            stride = 1
        elif stride == 1:
            stride = 2
        ax2.set_xticks(x.flatten()[:num_x:stride])
        formatString = "{{0:.{0:d}f}}$\degree$".format(
            upperXAxisTickLabelPrecision)
        ax2.set_xticklabels([formatString.format(member)
                             for member in secondXAxisData[::stride]])

    # add a third x-axis scale, if it was requested
    if thirdXAxisData is not None:
        ax3 = ax.twiny()
        ax3.set_facecolor(backgroundColor)
        ax3.set_xlabel(thirdXAxisLabel, **axis_font)
        ax3.set_xlim(ax.get_xlim())
        ax3.set_xticks(x.flatten()[:num_x:stride])
        ax3.set_xticklabels([formatString.format(member)
                             for member in thirdXAxisData[::stride]])
        ax3.spines['top'].set_position(('outward', 36))

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    if fileout is not None:
        plt.close()

    return  # }}}


def setup_colormap(config, configSectionName, suffix=''):

    '''
    Set up a colormap from the registry

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    configSectionName : str
        name of config section

    suffix: str, optional
        suffix of colormap related options

    Returns
    -------
    colormapDict : dict
        A dictionary of colormap information.

        'colormap' specifies the name of the new colormap

        'norm' is a matplotlib norm object used to normalize the colormap

        'levels' is an array of contour levels or ``None`` if not using indexed
        color map

        'ticks' is an array of values where ticks should be placed

        'contours' is an array of contour values to plot or ``None`` if none
        have been specified

        'lineWidth' is the width of contour lines or ``None`` if not specified

        'lineColor' is the color of contour lines or ``None`` if not specified
    '''
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani, Greg Streletz

    _register_custom_colormaps()

    colormapType = config.get(configSectionName,
                              'colormapType{}'.format(suffix))
    if colormapType == 'indexed':
        (colormap, norm, levels, ticks) = _setup_indexed_colormap(
            config, configSectionName, suffix=suffix)
    elif colormapType == 'continuous':
        (colormap, norm, ticks) = _setup_colormap_and_norm(
            config, configSectionName, suffix=suffix)
        levels = None
    else:
        raise ValueError('config section {} option colormapType{} is not '
                         '"indexed" or "continuous"'.format(
                                 configSectionName, suffix))

    option = 'contourLevels{}'.format(suffix)
    if config.has_option(configSectionName, option):
        contours = config.getExpression(configSectionName,
                                        option,
                                        usenumpyfunc=True)
        if contours == 'none':
            contours = None
    else:
        contours = None

    option = 'contourThickness{}'.format(suffix)
    if config.has_option(configSectionName, option):
        lineWidth = config.getfloat(configSectionName, option)
    else:
        lineWidth = None

    option = 'contourColor{}'.format(suffix)
    if config.has_option(configSectionName, option):
        lineColor = config.get(configSectionName, option)
    else:
        lineColor = None

    return {'colormap': colormap, 'norm': norm, 'levels': levels,
            'ticks': ticks, 'contours': contours, 'lineWidth': lineWidth,
            'lineColor': lineColor}


def plot_xtick_format(calendar, minDays, maxDays, maxXTicks, yearStride=None):
    '''
    Formats tick labels and positions along the x-axis for time series
    / index plots

    Parameters
    ----------
    calendar : str
        the calendar to use for formatting the time axis

    minDays : float
        start time for labels

    maxDays : float
        end time for labels

    maxXTicks : int
        the maximum number of tick marks to display, used to sub-sample ticks
        if there are too many

    yearStride : int, optional
        the number of years to skip over between ticks
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    ax = plt.gca()

    start = days_to_datetime(np.amin(minDays), calendar=calendar)
    end = days_to_datetime(np.amax(maxDays), calendar=calendar)

    if yearStride is not None or end.year - start.year > maxXTicks/2:
        if yearStride is None:
            yearStride = 1
        else:
            maxXTicks = None
        major = [date_to_days(year=year, calendar=calendar)
                 for year in np.arange(start.year, end.year+1, yearStride)]
        formatterFun = partial(_date_tick, calendar=calendar,
                               includeMonth=False)
    else:
        # add ticks for months
        major = []
        for year in range(start.year, end.year+1):
            for month in range(1, 13):
                major.append(date_to_days(year=year, month=month,
                                          calendar=calendar))
        formatterFun = partial(_date_tick, calendar=calendar,
                               includeMonth=True)

    ax.xaxis.set_major_locator(FixedLocator(major, maxXTicks))
    ax.xaxis.set_major_formatter(FuncFormatter(formatterFun))

    plt.setp(ax.get_xticklabels(), rotation=30)

    plt.autoscale(enable=True, axis='x', tight=True)


def _setup_colormap_and_norm(config, configSectionName, suffix=''):

    '''
    Set up a colormap from the registry

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    configSectionName : str
        name of config section

    suffix: str, optional
        suffix of colormap related options

    Returns
    -------
    colormap : srt
        new colormap

    norm : ``mapplotlib.colors.Normalize``
        the norm used to normalize the colormap

    ticks : array of float
        the tick marks on the colormap
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    _register_custom_colormaps()

    colormap = plt.get_cmap(config.get(configSectionName,
                                       'colormapName{}'.format(suffix)))

    normType = config.get(configSectionName, 'normType{}'.format(suffix))

    kwargs = config.getExpression(configSectionName,
                                  'normArgs{}'.format(suffix))

    if normType == 'symLog':
        norm = cols.SymLogNorm(**kwargs)
    elif normType == 'log':
        norm = cols.LogNorm(**kwargs)
    elif normType == 'linear':
        norm = cols.Normalize(**kwargs)
    else:
        raise ValueError('Unsupported norm type {} in section {}'.format(
            normType, configSectionName))

    try:
        ticks = config.getExpression(
                configSectionName, 'colorbarTicks{}'.format(suffix),
                usenumpyfunc=True)
    except(configparser.NoOptionError):
        ticks = None

    return (colormap, norm, ticks)


def _setup_indexed_colormap(config, configSectionName, suffix=''):

    '''
    Set up a colormap from the registry

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    configSectionName : str
        name of config section

    suffix: str, optional
        suffix of colormap related options

    colorMapType

    Returns
    -------
    colormap : srt
        new colormap

    norm : ``mapplotlib.colors.Normalize``
        the norm used to normalize the colormap

    ticks : array of float
        the tick marks on the colormap
    '''
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani, Greg Streletz

    colormap = plt.get_cmap(config.get(configSectionName,
                                       'colormapName{}'.format(suffix)))

    indices = config.getExpression(configSectionName,
                                   'colormapIndices{}'.format(suffix),
                                   usenumpyfunc=True)

    try:
        levels = config.getExpression(
                configSectionName, 'colorbarLevels{}'.format(suffix),
                usenumpyfunc=True)
    except(configparser.NoOptionError):
        levels = None

    if levels is not None:
        # set under/over values based on the first/last indices in the colormap
        underColor = colormap(indices[0])
        overColor = colormap(indices[-1])
        if len(levels)+1 == len(indices):
            # we have 2 extra values for the under/over so make the colormap
            # without these values
            indices = indices[1:-1]
        elif len(levels)-1 != len(indices):
            # indices list must be either one element shorter
            # or one element longer than colorbarLevels list
            raise ValueError('length mismatch between indices and '
                             'colorbarLevels')
        colormap = cols.ListedColormap(colormap(indices),
                                       'colormapName{}'.format(suffix))
        colormap.set_under(underColor)
        colormap.set_over(overColor)

    norm = cols.BoundaryNorm(levels, colormap.N)

    try:
        ticks = config.getExpression(
                configSectionName, 'colorbarTicks{}'.format(suffix),
                usenumpyfunc=True)
    except(configparser.NoOptionError):
        ticks = levels

    return (colormap, norm, levels, ticks)


def _date_tick(days, pos, calendar='gregorian', includeMonth=True):
    days = np.maximum(days, 0.)
    date = days_to_datetime(days, calendar)
    if includeMonth:
        return '{:04d}-{:02d}'.format(date.year, date.month)
    else:
        return '{:04d}'.format(date.year)


def _register_custom_colormaps():
    name = 'ferret'
    backgroundColor = (0.9, 0.9, 0.9)

    red = np.array([[0, 0.6],
                    [0.15, 1],
                    [0.35, 1],
                    [0.65, 0],
                    [0.8, 0],
                    [1, 0.75]])

    green = np.array([[0, 0],
                      [0.1, 0],
                      [0.35, 1],
                      [1, 0]])

    blue = np.array([[0, 0],
                     [0.5, 0],
                     [0.9, 0.9],
                     [1, 0.9]])

    colorCount = 21
    colorList = np.ones((colorCount, 4), float)
    colorList[:, 0] = np.interp(np.linspace(0, 1, colorCount),
                                red[:, 0], red[:, 1])
    colorList[:, 1] = np.interp(np.linspace(0, 1, colorCount),
                                green[:, 0], green[:, 1])
    colorList[:, 2] = np.interp(np.linspace(0, 1, colorCount),
                                blue[:, 0], blue[:, 1])
    colorList = colorList[::-1, :]

    colorMap = cols.LinearSegmentedColormap.from_list(
        name, colorList, N=255)

    colorMap.set_bad(backgroundColor)
    plt.register_cmap(name, colorMap)

    name = 'erdc_iceFire_H'

    colorArray = np.array([
        [-1, 4.05432e-07, 0, 5.90122e-06],
        [-0.87451, 0, 0.120401, 0.302675],
        [-0.74902, 0, 0.216583, 0.524574],
        [-0.623529, 0.0552475, 0.345025, 0.6595],
        [-0.498039, 0.128047, 0.492588, 0.720288],
        [-0.372549, 0.188955, 0.641309, 0.792092],
        [-0.247059, 0.327673, 0.784935, 0.873434],
        [-0.121569, 0.60824, 0.892164, 0.935547],
        [0.00392157, 0.881371, 0.912178, 0.818099],
        [0.129412, 0.951407, 0.835621, 0.449279],
        [0.254902, 0.904481, 0.690489, 0],
        [0.380392, 0.85407, 0.510864, 0],
        [0.505882, 0.777093, 0.33018, 0.00088199],
        [0.631373, 0.672862, 0.139087, 0.00269398],
        [0.756863, 0.508815, 0, 0],
        [0.882353, 0.299417, 0.000366289, 0.000547829],
        [1, 0.0157519, 0.00332021, 4.55569e-08]], float)

    colorCount = 255
    colorList = np.ones((colorCount, 4), float)
    x = colorArray[:, 0]
    for cIndex in range(3):
        colorList[:, cIndex] = np.interp(
            np.linspace(-1., 1., colorCount),
            x, colorArray[:, cIndex+1])

    colorMap = cols.LinearSegmentedColormap.from_list(
        name, colorList, N=255)

    plt.register_cmap(name, colorMap)

    name = 'erdc_iceFire_L'

    colorArray = np.array([
        [-1, 0.870485, 0.913768, 0.832905],
        [-0.87451, 0.586919, 0.887865, 0.934003],
        [-0.74902, 0.31583, 0.776442, 0.867858],
        [-0.623529, 0.18302, 0.632034, 0.787722],
        [-0.498039, 0.117909, 0.484134, 0.713825],
        [-0.372549, 0.0507239, 0.335979, 0.654741],
        [-0.247059, 0, 0.209874, 0.511832],
        [-0.121569, 0, 0.114689, 0.28935],
        [0.00392157, 0.0157519, 0.00332021, 4.55569e-08],
        [0.129412, 0.312914, 0, 0],
        [0.254902, 0.520865, 0, 0],
        [0.380392, 0.680105, 0.15255, 0.0025996],
        [0.505882, 0.785109, 0.339479, 0.000797922],
        [0.631373, 0.857354, 0.522494, 0],
        [0.756863, 0.910974, 0.699774, 0],
        [0.882353, 0.951921, 0.842817, 0.478545],
        [1, 0.881371, 0.912178, 0.818099]], float)

    colorCount = 255
    colorList = np.ones((colorCount, 4), float)
    x = colorArray[:, 0]
    for cIndex in range(3):
        colorList[:, cIndex] = np.interp(
            np.linspace(-1., 1., colorCount),
            x, colorArray[:, cIndex+1])

    colorMap = cols.LinearSegmentedColormap.from_list(
        name, colorList, N=255)

    plt.register_cmap(name, colorMap)

    name = 'BuOr'
    colors1 = plt.cm.PuOr(np.linspace(0., 1, 256))
    colors2 = plt.cm.RdBu(np.linspace(0, 1, 256))

    # combine them and build a new colormap, just the orange from the first
    # and the blue from the second
    colorList = np.vstack((colors1[0:128, :], colors2[128:256, :]))
    # reverse the order
    colorList = colorList[::-1, :]
    colorMap = cols.LinearSegmentedColormap.from_list(name, colorList)

    plt.register_cmap(name, colorMap)

    name = 'Maximenko'
    colorArray = np.array([
        [-1, 0., 0.45882352941, 0.76470588235],
        [-0.666667, 0., 0.70196078431, 0.90588235294],
        [-0.333333, 0.3294117647, 0.87058823529, 1.],
        [0., 0.76470588235, 0.94509803921, 0.98039215686],
        [0.333333, 1., 1., 0.],
        [0.666667, 1., 0.29411764705, 0.],
        [1, 1., 0., 0.]], float)

    colorCount = 255
    colorList = np.ones((colorCount, 4), float)
    x = colorArray[:, 0]
    for cIndex in range(3):
        colorList[:, cIndex] = np.interp(
            np.linspace(-1., 1., colorCount),
            x, colorArray[:, cIndex+1])

    colorMap = cols.LinearSegmentedColormap.from_list(
        name, colorList, N=255)

    plt.register_cmap(name, colorMap)

    # add the cmocean color maps
    mapNames = list(cmocean.cm.cmapnames)
    # don't bother with gray (already exists, I think)
    mapNames.pop(mapNames.index('gray'))
    for mapName in mapNames:
        plt.register_cmap(mapName, getattr(cmocean.cm, mapName))

    # add Scientific Colour-Maps 3.0 from
    # http://www.fabiocrameri.ch/colourmaps.php

    for mapName in ['berlin', 'bilbao', 'broc', 'cork', 'davos', 'devon',
                    'grayC', 'lajolla', 'lapaz', 'lisbon', 'oleron', 'oslo',
                    'roma', 'tofino', 'tokyo', 'turku', 'vik']:

        xmlFile = pkg_resources.resource_filename(
            __name__, 'ColourMapSuite3/{}/{}.xml'.format(mapName, mapName))
        _read_xml_colormap(xmlFile, mapName)

    # add SciVisColor colormaps from
    # https://sciviscolor.org/home/colormaps/

    for mapName in ['3wave-yellow-grey-blue', '3Wbgy5',
                    '4wave-grey-red-green-mgreen',  '5wave-yellow-brown-blue',
                    'blue-1', 'blue-3', 'blue-6', 'blue-8', 'blue-orange-div',
                    'brown-2', 'brown-5', 'brown-8', 'green-1', 'green-4',
                    'green-7', 'green-8', 'orange-5', 'orange-6',
                    'orange-green-blue-gray', 'purple-7', 'purple-8', 'red-1',
                    'red-3', 'red-4', 'yellow-1', 'yellow-7']:

        xmlFile = pkg_resources.resource_filename(
            __name__, 'SciVisColorColormaps/{}.xml'.format(mapName))
        _read_xml_colormap(xmlFile, mapName)


def _read_xml_colormap(xmlFile, mapName):
    '''Read in an XML colormap'''

    xml = ET.parse(xmlFile)

    root = xml.getroot()
    colormap = root.findall('ColorMap')
    if len(colormap) > 0:
        colormap = colormap[0]
        colorDict = {'red': [], 'green': [], 'blue': []}
        for point in colormap.findall('Point'):
            x = float(point.get('x'))
            color = [float(point.get('r')), float(point.get('g')),
                     float(point.get('b'))]
            colorDict['red'].append((x, color[0], color[0]))
            colorDict['green'].append((x, color[1], color[1]))
            colorDict['blue'].append((x, color[2], color[2]))
        cmap = LinearSegmentedColormap(mapName,  colorDict, 256)

        plt.register_cmap(mapName, cmap)


def _plot_color_gradients():
    '''from https://matplotlib.org/tutorials/colors/colormaps.html'''

    cmap_list = [m for m in plt.cm.cmap_d if not m.endswith("_r")]

    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))

    nrows = len(cmap_list)

    fig, axes = plt.subplots(figsize=(7.2, 0.25*nrows), nrows=nrows)
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.35, right=0.99)

    for ax, name in zip(axes, cmap_list):
        ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))
        pos = list(ax.get_position().bounds)
        x_text = pos[0] - 0.01
        y_text = pos[1] + pos[3]/2.
        fig.text(x_text, y_text, name, va='center', ha='right', fontsize=10)

    # Turn off *all* ticks & spines, not just the ones with colormaps.
    for ax in axes:
        ax.set_axis_off()

    plt.savefig('colormaps.png', dpi=100)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
