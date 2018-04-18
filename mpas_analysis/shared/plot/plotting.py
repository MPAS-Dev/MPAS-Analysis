# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#
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

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import xarray as xr
import pandas as pd
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import FuncFormatter, FixedLocator
import numpy as np
from functools import partial
from mpl_toolkits.axes_grid1 import make_axes_locatable

from mpas_analysis.shared.timekeeping.utility import days_to_datetime, \
    date_to_days

from mpas_analysis.shared.constants import constants

from six.moves import configparser


def timeseries_analysis_plot(config, dsvalues, N, title, xlabel, ylabel,
                             fileout, lineStyles, lineWidths, legendText,
                             calendar, maxPoints=None, titleFontSize=None,
                             figsize=(15, 6), dpi=None, maxXTicks=20,
                             obsMean=None, obsUncertainty=None,
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

    lineStyles, lineWidths, legendText : list of str
        control line style/width and corresponding legend text

    calendar : str
        the calendar to use for formatting the time axis

    maxPoints : list of {None, int}
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

        plt.plot(mean['Time'].values, mean.values,
                 lineStyles[dsIndex],
                 linewidth=lineWidths[dsIndex],
                 label=legendText[dsIndex])

    if obsMean is not None:
        obsCount = len(obsMean)
        assert(len(obsUncertainty) == obsCount)

        # space the observations along the time line, leaving gaps at either
        # end
        start = np.amin(minDays)
        end = np.amax(maxDays)
        obsTimes = np.linspace(start, end, obsCount+2)[1:-1]
        obsSymbols = ['o', '^', 's', 'D', '*']
        for iObs in range(obsCount):
            if obsMean[iObs] is not None:
                plt.errorbar(obsTimes[iObs], obsMean[iObs],
                             yerr=obsUncertainty[iObs],
                             fmt=obsSymbols[np.mod(iObs, len(obsSymbols))],
                             ecolor='k',
                             capthick=2, label=obsLegend[iObs])

    plt.legend(loc=legendLocation)

    ax = plt.gca()

    if titleFontSize is None:
        titleFontSize = config.get('plot', 'titleFontSize')
    axis_font = {'size': config.get('plot', 'axisFontSize')}
    title_font = {'size': titleFontSize,
                  'color': config.get('plot', 'titleFontColor'),
                  'weight': config.get('plot', 'titleFontWeight')}

    plot_xtick_format(plt, calendar, minDays, maxDays, maxXTicks)

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

    if not config.getboolean('plot', 'displayToScreen'):
        plt.close()


def timeseries_analysis_plot_polar(config, dsvalues, N, title,
                                   fileout, lineStyles, lineWidths,
                                   legendText, titleFontSize=None,
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

    lineStyles, lineWidths, legendText : list of str
        control line style/width and corresponding legend text

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
    for dsIndex in range(len(dsvalues)):
        dsvalue = dsvalues[dsIndex]
        if dsvalue is None:
            continue
        mean = pd.Series.rolling(dsvalue.to_pandas(), N, center=True).mean()
        mean = xr.DataArray.from_series(mean)
        minDays.append(mean.Time.min())
        maxDays.append(mean.Time.max())
        plt.polar((mean['Time']/365.0)*np.pi*2.0, mean,
                  lineStyles[dsIndex],
                  linewidth=lineWidths[dsIndex],
                  label=legendText[dsIndex])
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

    if not config.getboolean('plot', 'displayToScreen'):
        plt.close()


def plot_polar_comparison(
        config,
        Lons,
        Lats,
        modelArray,
        refArray,
        diffArray,
        colorMapSectionName,
        fileout,
        title=None,
        plotProjection='npstere',
        latmin=50.0,
        lon0=0,
        modelTitle='Model',
        refTitle='Observations',
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

    modelArray, refArray : float arrays
        model and observational or reference run data sets

    diffArray : float array
        difference between modelArray and refArray

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

    refTitle : str, optional
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
        x, y = m(Lons, Lats)  # compute map proj coordinates

        ax.set_title(title, y=1.06, **axis_font)
        m.drawcoastlines()
        m.fillcontinents(color='grey', lake_color='white')
        m.drawparallels(np.arange(-80., 81., 10.))
        m.drawmeridians(np.arange(-180., 181., 20.),
                        labels=[True, False, True, True])

        if levels is None:
            plotHandle = m.pcolormesh(x, y, field, cmap=colormap, norm=norm)
        else:
            plotHandle = m.contourf(x, y, field, cmap=colormap, norm=norm,
                                    levels=levels)

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

    if refArray is None:
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
    axis_font = {'size': config.get('plot', 'axisFontSize')}

    ax = plt.subplot(subplots[0])
    do_subplot(ax=ax, field=modelArray, title=modelTitle, **dictModelRef)

    if refArray is not None:
        ax = plt.subplot(subplots[1])
        do_subplot(ax=ax, field=refArray, title=refTitle, **dictModelRef)

        ax = plt.subplot(subplots[2])
        do_subplot(ax=ax, field=diffArray, title=diffTitle, **dictDiff)

    plt.tight_layout(pad=4.)
    if vertical:
        plt.subplots_adjust(top=0.9)

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    if not config.getboolean('plot', 'displayToScreen'):
        plt.close()


def plot_global_comparison(
        config,
        Lons,
        Lats,
        modelArray,
        refArray,
        diffArray,
        colorMapSectionName,
        fileout,
        title=None,
        modelTitle='Model',
        refTitle='Observations',
        diffTitle='Model-Observations',
        cbarlabel='units',
        titleFontSize=None,
        figsize=(8, 13),
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

    modelArray, refArray : float arrays
        model and observational or reference run data sets

    diffArray : float array
        difference between modelArray and refArray

    colorMapSectionName : str
        section name in ``config`` where color map info can be found.

    fileout : str
        the file name to be written

    title : str, optional
        the subtitle of the plot

    modelTitle : str, optional
        title of the model panel

    refTitle : str, optional
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
        the color contour lines (if specified)
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani

    def plot_panel(title, array, colormap, norm, levels, ticks, contours,
                   lineWidth, lineColor):
        plt.title(title, y=1.06, **axis_font)
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
    fig = plt.figure(figsize=figsize, dpi=dpi)
    if (title is not None):
        if titleFontSize is None:
            titleFontSize = config.get('plot', 'titleFontSize')
        title_font = {'size': titleFontSize,
                      'color': config.get('plot', 'titleFontColor'),
                      'weight': config.get('plot', 'titleFontWeight')}
        fig.suptitle(title, y=0.95, **title_font)
    axis_font = {'size': config.get('plot', 'axisFontSize')}

    m = Basemap(projection='cyl', llcrnrlat=-85, urcrnrlat=86, llcrnrlon=-180,
                urcrnrlon=181, resolution='l')
    x, y = m(Lons, Lats)  # compute map proj coordinates

    dictModelRef = setup_colormap(config, colorMapSectionName, suffix='Result')
    dictDiff = setup_colormap(config, colorMapSectionName, suffix='Difference')

    if refArray is not None:
        plt.subplot(3, 1, 1)

    plot_panel(modelTitle, modelArray, **dictModelRef)

    if refArray is not None:
        plt.subplot(3, 1, 2)
        plot_panel(refTitle, refArray, **dictModelRef)

        plt.subplot(3, 1, 3)
        plot_panel(diffTitle, diffArray, **dictDiff)

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    if not config.getboolean('plot', 'displayToScreen'):
        plt.close()


def plot_polar_projection_comparison(
        config,
        x,
        y,
        landMask,
        modelArray,
        refArray,
        diffArray,
        fileout,
        colorMapSectionName,
        title=None,
        modelTitle='Model',
        refTitle='Observations',
        diffTitle='Model-Observations',
        cbarlabel='units',
        titleFontSize=None,
        figsize=None,
        dpi=None,
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

    landMask : numpy ndarrays
        model and observational or reference run data sets

    modelArray, refArray : numpy ndarrays
        model and observational or reference run data sets

    diffArray : float array
        difference between modelArray and refArray

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

    refTitle : str, optional
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
    # Xylar Asay-Davis

    def plot_panel(ax, title, array, colormap, norm, levels, ticks, contours,
                   lineWidth, lineColor):
        plt.title(title, y=1.06, **axis_font)

        if levels is None:
            plotHandle = plt.pcolormesh(x, y, array, cmap=colormap, norm=norm)
        else:
            plotHandle = plt.contourf(x, y, array, cmap=colormap, norm=norm,
                                      levels=levels, extend='both')

        plt.pcolormesh(x, y, landMask, cmap=landColorMap)
        plt.contour(xCenter, yCenter, landMask.mask, (0.5,), colors='k',
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

    if refArray is None:
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

    # set up figure
    fig = plt.figure(figsize=figsize, dpi=dpi)
    if (title is not None):
        if titleFontSize is None:
            titleFontSize = config.get('plot', 'titleFontSize')
        title_font = {'size': titleFontSize,
                      'color': config.get('plot', 'titleFontColor'),
                      'weight': config.get('plot', 'titleFontWeight')}
        fig.suptitle(title, y=0.95, **title_font)
    axis_font = {'size': config.get('plot', 'axisFontSize')}

    # set up land colormap
    colorList = [(0.8, 0.8, 0.8), (0.8, 0.8, 0.8)]
    landColorMap = cols.LinearSegmentedColormap.from_list('land', colorList)

    # locations of centers for land contour
    xCenter = 0.5*(x[1:] + x[0:-1])
    yCenter = 0.5*(y[1:] + y[0:-1])

    ax = plt.subplot(subplots[0])
    plot_panel(ax, modelTitle, modelArray, **dictModelRef)

    if refArray is not None:
        ax = plt.subplot(subplots[1])
        plot_panel(ax, refTitle, refArray, **dictModelRef)

        ax = plt.subplot(subplots[2])
        plot_panel(ax, diffTitle, diffArray, **dictDiff)

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    if not config.getboolean('plot', 'displayToScreen'):
        plt.close()


def plot_1D(config, xArrays, fieldArrays, errArrays,
            lineColors, lineWidths, legendText,
            title=None, xlabel=None, ylabel=None,
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

    lineColors, legendText : list of str
        control line color and legend

    lineWidths : list of int
        control line width

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

    for dsIndex in range(len(xArrays)):
        xArray = xArrays[dsIndex]
        fieldArray = fieldArrays[dsIndex]
        errArray = errArrays[dsIndex]
        if xArray is None:
            continue
        if errArray is None:
            plt.plot(xArray, fieldArray,
                     color=lineColors[dsIndex],
                     linewidth=lineWidths[dsIndex],
                     label=legendText[dsIndex])
        else:
            plt.plot(xArray, fieldArray,
                     color=lineColors[dsIndex],
                     linewidth=lineWidths[dsIndex],
                     label=legendText[dsIndex])
            plt.fill_between(xArray, fieldArray, fieldArray+errArray,
                             facecolor=lineColors[dsIndex], alpha=0.2)
            plt.fill_between(xArray, fieldArray, fieldArray-errArray,
                             facecolor=lineColors[dsIndex], alpha=0.2)
    plt.grid()
    plt.axhline(0.0, linestyle='-', color='k')  # horizontal lines
    if dsIndex > 0:
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

    if not config.getboolean('plot', 'displayToScreen'):
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
        fileout='moc.png',
        figsize=(10, 4),
        dpi=None,
        xLim=None,
        yLim=None,
        linewidths=2,
        invertYAxis=True,
        xArrayIsTime=False,
        N=None,
        maxXTicks=20,
        calendar='gregorian'):  # {{{

    """
    Plots a data set as a x distance (latitude, longitude,
    or spherical distance) vs depth map (vertical section).

    Or, if xArrayIsTime is True, plots data set on a vertical
    Hovmoller plot (depth vs. time).

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

    linewidths : int, optional
        linewidths for contours

    invertYAxis : logical, optional
        if True, invert Y axis

    xArrayIsTime : logical, optional
        if True, format X axis for time

    N : int, optional
        the number of points over which to perform a moving average
        NOTE: this option is mostly intended for use when xArrayIsTime is True,
        although it will work with other data as well.  Also, the moving
        average calculation is based on number of points, not actual x axis
        values, so for best results, the values in the xArray should be equally
        spaced.

    maxXTicks : int, optional
        the maximum number of tick marks that will be allowed along the x axis.
        This may need to be adjusted depending on the figure size and aspect
        ratio.  NOTE:  maxXTicks is only used if xArrayIsTime is True

    calendar : str, optional
        the calendar to use for formatting the time axis
        NOTE:  calendar is only used if xArrayIsTime is True
    """
    # Authors
    # -------
    # Milena Veneziani, Mark Petersen, Xylar Asay-Davis, Greg Streletz

    # verify that the dimensions of fieldArray are consistent with those of
    # xArray and depthArray
    if len(xArray) != fieldArray.shape[1]:
        raise ValueError('size mismatch between xArray and fieldArray')
    elif len(depthArray) != fieldArray.shape[0]:
        raise ValueError('size mismatch between depthArray and fieldArray')

    # set up figure
    if dpi is None:
        dpi = config.getint('plot', 'dpi')
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

    x, y = np.meshgrid(xArray, depthArray)  # change to zMid

    colormapDict = setup_colormap(config, colorMapSectionName, suffix=suffix)

    cs = plt.contourf(x, y, fieldArray, cmap=colormapDict['colormap'],
                      norm=colormapDict['norm'],
                      levels=colormapDict['levels'], extend='both')

    contourLevels = colormapDict['contours']
    if contourLevels is not None:
        if len(contourLevels) == 0:
            # automatic calculation of contour levels
            contourLevels = None
        plt.contour(x, y, fieldArray, levels=contourLevels, colors='k',
                    linewidths=linewidths)

    cbar = plt.colorbar(cs, orientation='vertical', spacing='uniform',
                        ticks=colormapDict['ticks'],
                        boundaries=colormapDict['ticks'])

    if colorbarLabel is not None:
        cbar.set_label(colorbarLabel)

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

    if xArrayIsTime:
        minDays = [xArray[0]]
        maxDays = [xArray[-1]]
        plot_xtick_format(plt, calendar, minDays, maxDays, maxXTicks)

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    if not config.getboolean('plot', 'displayToScreen'):
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

    colorMapType

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

    if config.has_option(configSectionName,
                         'colormapIndices{}'.format(suffix)):
        (colormap, norm, levels, ticks) = _setup_indexed_colormap(
            config, configSectionName, suffix=suffix)
    elif config.has_option(configSectionName, 'normType{}'.format(suffix)):
        (colormap, norm, ticks) = _setup_colormap_and_norm(
            config, configSectionName, suffix=suffix)
        levels = None
    else:
        raise ValueError('config section {} contains neither the info '
                         'for an indexed color map nor for computing a '
                         'norm'.format(configSectionName))

    option = 'contourLevels{}'.format(suffix)
    if config.has_option(configSectionName, option):
        contours = config.getExpression(configSectionName,
                                        option,
                                        usenumpyfunc=True)
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


def plot_xtick_format(plt, calendar, minDays, maxDays, maxXTicks):
    '''
    Formats tick labels and positions along the x-axis for time series
    / index plots

    Parameters
    ----------
    plt : plt handle on which to change ticks

    calendar : specified calendar for the plot

    minDays : start time for labels

    maxDays : end time for labels
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    ax = plt.gca()

    start = days_to_datetime(np.amin(minDays), calendar=calendar)
    end = days_to_datetime(np.amax(maxDays), calendar=calendar)

    if (end.year - start.year > maxXTicks/2):
        major = [date_to_days(year=year, calendar=calendar)
                 for year in np.arange(start.year, end.year+1)]
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

    norm : ``SymLogNorm`` object
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

    norm : ``SymLogNorm`` object
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

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
