"""
Plotting utilities, including routines for plotting:
    * time series (and comparing with reference data sets)
    * regridded horizontal fields (and comparing with reference data sets)
    * vertical sections on native grid
    * NINO34 time series and spectra

Authors
-------
Xylar Asay-Davis, Milena Veneziani, Luke Van Roekel

Last Modified
-------------
04/07/2017
"""

import matplotlib.pyplot as plt
import matplotlib.colors as cols
import xarray as xr
import pandas as pd
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as cols
from matplotlib.ticker import FuncFormatter, FixedLocator
import numpy as np
from functools import partial

from ..timekeeping.utility import days_to_datetime, date_to_days

from ..constants import constants


def nino34_spectra_plot(config, f, ninoSpectra, confidence95, confidence99,
                        redNoiseSpectra, fObs, f30, ninoObs, conf95Obs, conf99Obs,
                        redNoiseObs, nino30yr, conf9530, conf9930, redNoise30,
                        title, modelTitle, obsTitle,
                        fileout, linewidths, xlabel='Period (years)',
                        ylabel=r'Power ($^o$C / cycles mo$^{-1}$)', titleFontSize=None,
                        figsize=(9, 21), dpi=300):
    """
    Plots the nino34 time series and power spectra in an image file
    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    f : numpy.array
        periods to plot on x-axis

    ninoSpectra : xarray.dataArray object
        nino34 power spectra

    confidence95 : numpy.array
        95% confidence level based on chi squared test

    confidence99 : numpy.array
        99% confidence level based on chi squared test

    redNoiseSpectra : numpy.array
        red noise fit to the ninoSpectra

    fObs : numpy.array
           periods to plot on x-axis for observations

    ninoObs : xarray.dataArray object
        nino34 power spectra from the full observational record

    conf95Obs : numpy.array
        95% confidence level based on chi squared for observations

    conf99Obs : numpy.array
        99% confidence level based on chi squared for observations

    redNoiseObs : numpy.array
        red noise fit to ninoObs

    nino30yr : xarray.dataArray object
        power spectra of the last 30 years of the observational record

    title : str
        the title of the plot

    modelTitle : str
        the title of model panel

    obsTitle : str
        the title of the obs panel

    xLabel, yLabel : str
        axis labels

    fileout : str
        the file name to be written

    linewidths : control line width

    titleFontSize : int, optional
        the size of the title font

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure

    Author
    ------
    Luke Van Roekel

    Last Modified
    -------------
    04/07/2017
    """

    fig = plt.figure(figsize=figsize, dpi=dpi)

    if titleFontSize is None:
        titleFontSize = config.get('plot', 'titleFontSize')

    axis_font = {'size': config.get('plot', 'axisFontSize')}
    title_font = {'size': titleFontSize,
                  'color': config.get('plot', 'titleFontColor'),
                  'weight': config.get('plot', 'titleFontWeight')}
    if title is not None:
        fig.suptitle(title, y=0.92, **title_font)

    ax1 = plt.subplot(3, 1, 1)

    plt.plot(fObs[2:-3], ninoObs[2:-3], 'k', linewidth=linewidths)
    plt.plot(fObs[2:-3], redNoiseObs[2:-3], 'r', linewidth=linewidths)
    plt.plot(fObs[2:-3], conf95Obs[2:-3], 'b', linewidth=linewidths)
    plt.plot(fObs[2:-3], conf99Obs[2:-3], 'g', linewidth=linewidths)
    plt.xlim(10, 1)

    plt.legend(['Nino34 spectra (Full Record)', 'Red noise fit',
               '95% confidence threshold', '99% confidence threshold'],
               loc='upper right')
    maxObs = _plot_size_y_axis(plt, fObs, c1=conf99Obs, c2=redNoiseObs)
    max30 = _plot_size_y_axis(plt, f30, c1=conf9930, c2=redNoise30)
    maxModel = _plot_size_y_axis(plt, f, c1=ninoSpectra.values, c2=confidence99,
                                 c3=redNoiseSpectra)

    maxYval = max(maxObs, max30, maxModel)
    plt.ylim(0, 0.9*maxYval)

    if obsTitle is not None:
        plt.title(obsTitle+' (Full Record)', **title_font)
    if xlabel is not None:
        plt.xlabel(xlabel, **axis_font)
    if ylabel is not None:
        plt.ylabel(ylabel, **axis_font)

    ax2 = plt.subplot(3, 1, 2)

    plt.plot(f30[2:-3], nino30yr[2:-3], 'k', linewidth=linewidths)
    plt.plot(f30[2:-3], redNoise30[2:-3], 'r', linewidth=linewidths)
    plt.plot(f30[2:-3], conf9530[2:-3], 'b', linewidth=linewidths)
    plt.plot(f30[2:-3], conf9930[2:-3], 'g', linewidth=linewidths)
    plt.xlim(10, 1)
    plt.ylim(0, 0.9*maxYval)

    plt.legend(['Nino34 spectra (1976 - 2016)', 'Red noise fit',
               '95% confidence threshold', '99% confidence threshold'],
               loc='upper right')

    if obsTitle is not None:
        plt.title(obsTitle+' (1976-2016)', **title_font)
    if xlabel is not None:
        plt.xlabel(xlabel, **axis_font)
    if ylabel is not None:
        plt.ylabel(ylabel, **axis_font)

    ax3 = plt.subplot(3, 1, 3)
    plt.plot(f[2:-3], ninoSpectra[2:-3], 'k', linewidth=linewidths)
    plt.plot(f[2:-3], redNoiseSpectra[2:-3], 'r', linewidth=linewidths)
    plt.plot(f[2:-3], confidence95[2:-3], 'b', linewidth=linewidths)
    plt.plot(f[2:-3], confidence99[2:-3], 'g', linewidth=linewidths)
    plt.xlim(10, 1)
    plt.ylim(0, 0.9*maxYval)

    # add legend
    plt.legend(['Nino34 index spectra', 'Red noise fit',
               '95% confidence threshold', '99% confidence threshold'],
               loc='upper right')

    if modelTitle is not None:
        plt.title(modelTitle, **title_font)
    if xlabel is not None:
        plt.xlabel(xlabel, **axis_font)
    if ylabel is not None:
        plt.ylabel(ylabel, **axis_font)
    if fileout is not None:
        fig.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    if not config.getboolean('plot', 'displayToScreen'):
        plt.close()


def nino34_timeseries_plot(config, nino34Index, nino34Obs, nino3430, title,
                           modelTitle, obsTitle, fileout, linewidths,
                           calendar, xlabel='Time [years]', ylabel='[$^\circ$C]',
                           titleFontSize=None, figsize=(12, 28), dpi=300,
                           maxXTicks=20):
    """
    Plots the nino34 time series and power spectra in an image file

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    nino34Index : xarray.dataArray
        nino34 timeseries to plot

    nino34Obs : xarray.dataArray
        nino34 observation

    nino3430 : xarray.dataArray
        subset of nino34 observations

    title : str
        the title of the plot

    obsTitle : str
        title of observational plot

    modelTitle : str
        title of model plot

    xLabel, yLabel : str
        axis labels

    fileout : str
        the file name to be written

    lineWidths : list of str
        control line width

    titleFontSize : int, optional
        the size of the title font

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure

    maxXTicks : int, optional
        the maximum number of tick marks that will be allowed along the x axis.
        This may need to be adjusted depending on the figure size and aspect
        ratio.

    Author
    ------
    Luke Van Roekel

    Last Modified
    -------------
    04/07/2017
    """
    fig = plt.figure(figsize=figsize, dpi=dpi)

    if titleFontSize is None:
        titleFontSize = config.get('plot', 'titleFontSize')

    axis_font = {'size': config.get('plot', 'axisFontSize')}
    title_font = {'size': titleFontSize,
                  'color': config.get('plot', 'titleFontColor'),
                  'weight': config.get('plot', 'titleFontWeight')}
    if title is not None:
        fig.suptitle(title, y=0.92, **title_font)

    # Plot Nino34 Observation Time series
    plt.subplot(3, 1, 1)
    _plot_nino_timeseries(plt, nino34Obs[2:-3].values, nino34Obs.Time[2:-3].values,
                          xlabel, ylabel, obsTitle+' (Full Record)', calendar,
                          axis_font, linewidths, maxXTicks)

    # Plot subset of the observational data set
    plt.subplot(3, 1, 2)
    _plot_nino_timeseries(plt, nino3430.values, nino3430.Time.values,
                          xlabel, ylabel, obsTitle+' (1976 - 2016)', calendar,
                          axis_font, linewidths, maxXTicks)

    # Plot Nino34 model time series
    plt.subplot(3, 1, 3)
    _plot_nino_timeseries(plt, nino34Index[2:-3].values, nino34Index.Time[2:-3].values,
                          xlabel, ylabel, modelTitle, calendar, axis_font, linewidths,
                          maxXTicks)
    minDays = nino34Index.Time[2:-3].values.min()
    maxDays = nino34Index.Time[2:-3].values.max()

    _plot_xtick_format(plt, calendar, minDays, maxDays, maxXTicks)

    if fileout is not None:
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    if not config.getboolean('plot', 'displayToScreen'):
        plt.close()


def _plot_nino_timeseries(plt, ninoIndex, time, xlabel, ylabel,
                          panelTitle, calendar, axis_font, linewidths,
                          maxXTicks):
    '''
    Plot the nino time series on a subplot

    Parameters
    ----------
    ninoIndex : numpy.array
      nino34 Index values (can be obs or model)

    time : numpy.array
      time values for the nino index

    calendar : specified calendar for the plot

    maxXTicks : int, optional
        the maximum number of tick marks that will be allowed along the x axis.
        This may need to be adjusted depending on the figure size and aspect
        ratio.

    panelTitle : string
        string to label the subplot with

    xlabel : string
        string for x-axis label

    ylabel : string
        string for y-axis label

    Author
    ------
    Luke Van Roekel

    Last Modified
    -------------
    04/07/2017
    '''
    plt.title(panelTitle, y=1.06, **axis_font)
    y1 = ninoIndex
    nt = np.size(ninoIndex)

    y2 = np.zeros(nt)

    plt.plot(time, 0.4*np.ones(nt), '--k',
             linewidth=linewidths)
    plt.plot(time, -0.4*np.ones(nt), '--k',
             linewidth=linewidths)
    plt.fill_between(time, y1, y2, where=y1 > y2,
                     facecolor='red', interpolate=True, linewidth=0)
    plt.fill_between(time, y1, y2, where=y1 < y2,
                     facecolor='blue', interpolate=True, linewidth=0)

    if xlabel is not None:
        plt.xlabel(xlabel, **axis_font)
    if ylabel is not None:
        plt.ylabel(ylabel, **axis_font)


def timeseries_analysis_plot(config, dsvalues, N, title, xlabel, ylabel,
                             fileout, lineStyles, lineWidths, calendar,
                             titleFontSize=None, figsize=(15, 6), dpi=300,
                             maxXTicks=20):

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

    lineStyles, lineWidths : list of str
        control line style/width

    titleFontSize : int, optional
        the size of the title font

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure

    maxXTicks : int, optional
        the maximum number of tick marks that will be allowed along the x axis.
        This may need to be adjusted depending on the figure size and aspect
        ratio.

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani

    Last Modified
    -------------
    03/14/2017
    """
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
        plt.plot(mean['Time'], mean,
                 lineStyles[dsIndex],
                 linewidth=lineWidths[dsIndex])
    ax = plt.gca()
    # Add a y=0 line if y ranges between positive and negative values
    yaxLimits = ax.get_ylim()
    if yaxLimits[0]*yaxLimits[1] < 0:
        indgood = np.where(np.logical_not(np.isnan(mean)))
        x = mean['Time'][indgood]
        plt.plot(x, np.zeros(np.size(x)), 'k-', linewidth=1.2)

    _plot_xtick_format(plt, calendar, minDays, maxDays, maxXTicks)

    if titleFontSize is None:
        titleFontSize = config.get('plot', 'titleFontSize')

    axis_font = {'size': config.get('plot', 'axisFontSize')}
    title_font = {'size': titleFontSize,
                  'color': config.get('plot', 'titleFontColor'),
                  'weight': config.get('plot', 'titleFontWeight')}
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
                                   calendar, titleFontSize=None,
                                   figsize=(15, 6), dpi=300):

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

    lineStyles, lineWidths : list of str
        control line style/width

    titleFontSize : int, optional
        the size of the title font

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure

    Authors
    -------
    Adrian K. Turner

    Last Modified
    -------------
    03/15/2017
    """
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
                   linewidth=lineWidths[dsIndex])

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
            (((constants.daysInMonth[month-1] + \
               constants.daysInMonth[month]) * np.pi) / 365.0)

    ax.set_xticks(majorTickLocs)
    ax.set_xticklabels([])

    ax.set_xticks(minorTickLocs, minor=True)
    ax.set_xticklabels(constants.abrevMonthNames, minor=True)

    if titleFontSize is None:
        titleFontSize = config.get('plot', 'titleFontSize')

    axis_font = {'size': config.get('plot', 'axisFontSize')}
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
    obsArray,
    diffArray,
    cmapModelObs,
    clevsModelObs,
    cmapDiff,
    clevsDiff,
    fileout,
    title=None,
    plotProjection='npstere',
    latmin=50.0,
    lon0=0,
    modelTitle='Model',
    obsTitle='Observations',
    diffTitle='Model-Observations',
    cbarlabel='units',
    titleFontSize=None,
    figsize=(8, 22),
    dpi=300):

    """
    Plots a data set around either the north or south pole.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    Lons, Lats : float arrays
        longitude and latitude arrays

    modelArray, obsArray : float arrays
        model and observational data sets

    diffArray : float array
        difference between modelArray and obsArray

    cmapModelObs : str
        colormap of model and observations panel

    clevsModelObs : int array
        colorbar values for model and observations panel

    cmapDiff : str
        colormap of difference (bias) panel

    clevsDiff : int array
        colorbar values for difference (bias) panel

    fileout : str
        the file name to be written

    title : str, optional
        the subtitle of the plot

    plotProjection : str, optional
        Basemap projection for the plot

    modelTitle : str, optional
        title of the model panel

    obsTitle : str, optional
        title of the observations panel

    diffTitle : str, optional
        title of the difference (bias) panel

    cbarlabel : str, optional
        label on the colorbar

    titleFontSize : int, optional
        size of the title font

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani

    Last Modified
    -------------
    03/17/2017
    """

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

    m = Basemap(projection=plotProjection, boundinglat=latmin,
                lon_0=lon0, resolution='l')
    x, y = m(Lons, Lats)  # compute map proj coordinates

    normModelObs = cols.BoundaryNorm(clevsModelObs, cmapModelObs.N)
    normDiff = cols.BoundaryNorm(clevsDiff, cmapDiff.N)

    plt.subplot(3, 1, 1)
    plt.title(modelTitle, y=1.06, **axis_font)
    m.drawcoastlines()
    m.fillcontinents(color='grey', lake_color='white')
    m.drawparallels(np.arange(-80., 81., 10.))
    m.drawmeridians(np.arange(-180., 181., 20.),
                    labels=[True, True, True, True])
    cs = m.contourf(x, y, modelArray, cmap=cmapModelObs, norm=normModelObs,
                    spacing='uniform', levels=clevsModelObs)
    cbar = m.colorbar(cs, location='right', pad="15%", spacing='uniform',
                      ticks=clevsModelObs, boundaries=clevsModelObs)
    cbar.set_label(cbarlabel)

    plt.subplot(3, 1, 2)
    plt.title(obsTitle, y=1.06, **axis_font)
    m.drawcoastlines()
    m.fillcontinents(color='grey', lake_color='white')
    m.drawparallels(np.arange(-80., 81., 10.))
    m.drawmeridians(np.arange(-180., 181., 20.),
                    labels=[True, True, True, True])
    cs = m.contourf(x, y, obsArray, cmap=cmapModelObs, norm=normModelObs,
                    spacing='uniform', levels=clevsModelObs)
    cbar = m.colorbar(cs, location='right', pad="15%", spacing='uniform',
                      ticks=clevsModelObs, boundaries=clevsModelObs)
    cbar.set_label(cbarlabel)

    plt.subplot(3, 1, 3)
    plt.title(diffTitle, y=1.06, **axis_font)
    m.drawcoastlines()
    m.fillcontinents(color='grey', lake_color='white')
    m.drawparallels(np.arange(-80., 81., 10.))
    m.drawmeridians(np.arange(-180., 181., 20.),
                    labels=[True, True, True, True])
    cs = m.contourf(x, y, diffArray, cmap=cmapDiff, norm=normDiff,
                    spacing='uniform', levels=clevsDiff)
    cbar = m.colorbar(cs, location='right', pad="15%", spacing='uniform',
                      ticks=clevsDiff, boundaries=clevsModelObs)
    cbar.set_label(cbarlabel)

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    if not config.getboolean('plot', 'displayToScreen'):
        plt.close()


def plot_global_comparison(
    config,
    Lons,
    Lats,
    modelArray,
    obsArray,
    diffArray,
    cmapModelObs,
    clevsModelObs,
    cmapDiff,
    clevsDiff,
    fileout,
    title=None,
    modelTitle='Model',
    obsTitle='Observations',
    diffTitle='Model-Observations',
    cbarlabel='units',
    titleFontSize=None,
    figsize=(8, 12),
    dpi=300):

    """
    Plots a data set as a longitude/latitude map.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    Lons, Lats : float arrays
        longitude and latitude arrays

    modelArray, obsArray : float arrays
        model and observational data sets

    diffArray : float array
        difference between modelArray and obsArray

    cmapModelObs : str
        colormap of model and observations panel

    clevsModelObs : int array
        colorbar values for model and observations panel

    cmapDiff : str
        colormap of difference (bias) panel

    clevsDiff : int array
        colorbar values for difference (bias) panel

    fileout : str
        the file name to be written

    title : str, optional
        the subtitle of the plot

    modelTitle : str, optional
        title of the model panel

    obsTitle : str, optional
        title of the observations panel

    diffTitle : str, optional
        title of the difference (bias) panel

    cbarlabel : str, optional
        label on the colorbar

    titleFontSize : int, optional
        size of the title font

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani

    Last Modified
    -------------
    03/13/2017
    """

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

    m = Basemap(projection='cyl', llcrnrlat=-85, urcrnrlat=86, llcrnrlon=-180,
                urcrnrlon=181, resolution='l')
    x, y = m(Lons, Lats)  # compute map proj coordinates

    normModelObs = cols.BoundaryNorm(clevsModelObs, cmapModelObs.N)
    normDiff = cols.BoundaryNorm(clevsDiff, cmapDiff.N)

    plt.subplot(3, 1, 1)
    plt.title(modelTitle, y=1.06, **axis_font)
    m.drawcoastlines()
    m.fillcontinents(color='grey', lake_color='white')
    m.drawparallels(np.arange(-80., 80., 20.),
                    labels=[True, False, False, False])
    m.drawmeridians(np.arange(-180., 180., 60.),
                    labels=[False, False, False, True])
    cs = m.contourf(x, y, modelArray, cmap=cmapModelObs, norm=normModelObs,
                    spacing='uniform', levels=clevsModelObs, extend='both')
    cbar = m.colorbar(cs, location='right', pad="5%", spacing='uniform',
                      ticks=clevsModelObs, boundaries=clevsModelObs)
    cbar.set_label(cbarlabel)

    plt.subplot(3, 1, 2)
    plt.title(obsTitle, y=1.06, **axis_font)
    m.drawcoastlines()
    m.fillcontinents(color='grey', lake_color='white')
    m.drawparallels(np.arange(-80., 80., 20.),
                    labels=[True, False, False, False])
    m.drawmeridians(np.arange(-180., 180., 40.),
                    labels=[False, False, False, True])
    cs = m.contourf(x, y, obsArray, cmap=cmapModelObs, norm=normModelObs,
                    spacing='uniform', levels=clevsModelObs, extend='both')
    cbar = m.colorbar(cs, location='right', pad="5%", spacing='uniform',
                      ticks=clevsModelObs, boundaries=clevsModelObs)
    cbar.set_label(cbarlabel)

    plt.subplot(3, 1, 3)
    plt.title(diffTitle, y=1.06, **axis_font)
    m.drawcoastlines()
    m.fillcontinents(color='grey', lake_color='white')
    m.drawparallels(np.arange(-80., 80., 20.),
                    labels=[True, False, False, False])
    m.drawmeridians(np.arange(-180., 180., 40.),
                    labels=[False, False, False, True])
    cs = m.contourf(x, y, diffArray, cmap=cmapDiff, norm=normDiff,
                    spacing='uniform', levels=clevsDiff, extend='both')
    cbar = m.colorbar(cs, location='right', pad="5%", spacing='uniform',
                      ticks=clevsDiff, boundaries=clevsModelObs)
    cbar.set_label(cbarlabel)

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    if not config.getboolean('plot', 'displayToScreen'):
        plt.close()


def _date_tick(days, pos, calendar='gregorian', includeMonth=True):
    days = np.maximum(days, 0.)
    date = days_to_datetime(days, calendar)
    if includeMonth:
        return '{:04d}-{:02d}'.format(date.year, date.month)
    else:
        return '{:04d}'.format(date.year)


def plot_vertical_section(
    config,
    xArray,
    depthArray,
    fieldArray,
    colormapName,
    colorbarLevels,
    contourLevels,
    colorbarLabel=None,
    title=None,
    xlabel=None,
    ylabel=None,
    fileout='moc.png',
    figsize=(10, 4),
    dpi=300):  # {{{

    """
    Plots a data set as a x distance (latitude, longitude,
    or spherical distance) vs depth map (vertical section).

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    xArray : float array
        x array (latitude, longitude, or spherical distance)

    depthArray : float array
        depth array [m]

    fieldArray : float array
        field array to plot

    colormapName : str
        colormap of plot

    colorbarLevels : int array
        colorbar levels of plot

    contourLevels : int levels
        levels of contours to be drawn

    colorbarLabel : str, optional
        label of the colorbar

    title : str, optional
        title of plot

    xlabel, ylabel : str, optional
        label of x- and y-axis

    fileout : str, optional
        the file name to be written

    figsize : tuple of float, optional
        size of the figure in inches

    dpi : int, optional
        number of dots per inch of the figure

    Authors
    -------
    Milena Veneziani, Mark Petersen

    Last Modified
    -------------
    03/13/2017
    """

    # set up figure
    fig = plt.figure(figsize=figsize, dpi=dpi)

    x, y = np.meshgrid(xArray, depthArray)  # change to zMid

    normModelObs = cols.BoundaryNorm(colorbarLevels, colormapName.N)

    cs = plt.contourf(x, y, fieldArray, cmap=colormapName, norm=normModelObs,
                      spacing='uniform', levels=colorbarLevels, extend='both')
    plt.contour(x, y, fieldArray, levels=contourLevels[::2], colors='k')

    cbar = plt.colorbar(cs, orientation='vertical', spacing='uniform',
                        ticks=colorbarLevels, boundaries=colorbarLevels)
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

    plt.gca().invert_yaxis()

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

    Returns
    -------
    colormap : srt
        new colormap

    colorbarLevels : int array
        colorbar levels

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani

    Last modified
    -------------
    03/17/2017
    '''

    colormap = plt.get_cmap(config.get(configSectionName,
                                       'colormapName{}'.format(suffix)))
    indices = config.getExpression(configSectionName,
                                   'colormapIndices{}'.format(suffix))
    colorbarLevels = config.getExpression(configSectionName,
                                          'colorbarLevels{}'.format(suffix))

    # set under/over values based on the first/last indices in the colormap
    underColor = colormap(indices[0])
    overColor = colormap(indices[-1])
    if len(colorbarLevels)+1 == len(indices):
        # we have 2 extra values for the under/over so make the colormap
        # without these values
        indices = indices[1:-1]
    colormap = cols.ListedColormap(colormap(indices),
                                   'colormapName{}'.format(suffix))
    colormap.set_under(underColor)
    colormap.set_over(overColor)
    return (colormap, colorbarLevels)


def _plot_size_y_axis(plt, xaxisValues, **data):
    '''
    Resize the y-axis limit based on the curves being plotted

    Parameters
    ----------
    plt : plot handle

    xaxisValues : numpy.array
       Values plotted along the x-axis

    data : dictionary entries must be numpy.array
       data for curves on plot

    Author
    ------
    Luke Van Roekel

    Last modified
    -------------
    04/07/2017
    '''

    ax = plt.gca()
    xmin = ax.get_xlim()[0]
    xmax = ax.get_xlim()[1]

    # find period/frequency bounds for chosen xmin/xmax
    minIndex = np.abs(xaxisValues - xmin).argmin()
    maxIndex = np.abs(xaxisValues - xmax).argmin()

    # find maximum value of three curves plotted
    maxCurveVal = -1E20
    for key in data:
        maxTemp = data[key][minIndex:maxIndex].max()
        maxCurveVal = max(maxTemp, maxCurveVal)

    return maxCurveVal


def _plot_xtick_format(plt, calendar, minDays, maxDays, maxXTicks):
    '''
    Formats tick labels and positions along the x-axis for time series
    / index plots

    Parameters
    ----------
    plt : plt handle on which to change ticks

    calendar : specified calendar for the plot

    minDays : start time for labels

    maxDays : end time for labels

    Author
    ------
    Xylar Asay-Davis

    '''
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

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
