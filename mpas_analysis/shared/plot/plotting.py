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
from matplotlib.ticker import FuncFormatter, FixedLocator
import numpy as np
from functools import partial

from ..timekeeping.utility import days_to_datetime, date_to_days

from ..constants import constants


def nino34_spectra_plot(config, f, ninoSpectra,
                        confidence95, confidence99, redNoiseSpectra,
                        fObs, f30, ninoObs,
                        conf95Obs, conf99Obs, redNoiseObs,
                        nino30yr, conf9530, conf9930, redNoise30,
                        title, modelTitle, obsTitle,
                        fileout, linewidths, xlabel='Period (years)',
                        ylabel=r'Power ($^o$C / cycles mo$^{-1}$)',
                        titleFontSize=None, figsize=(9, 21), dpi=300):
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
    maxModel = _plot_size_y_axis(plt, f, c1=ninoSpectra.values,
                                 c2=confidence99, c3=redNoiseSpectra)

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
                           modelTitle, obsTitle, fileout, linewidths, calendar,
                           xlabel='Time [years]', ylabel='[$^\circ$C]',
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
    _plot_nino_timeseries(plt, nino34Obs[2:-3].values,
                          nino34Obs.Time[2:-3].values,
                          xlabel, ylabel, obsTitle+' (Full Record)',
                          calendar, axis_font, linewidths, maxXTicks)

    # Plot subset of the observational data set
    plt.subplot(3, 1, 2)
    _plot_nino_timeseries(plt, nino3430.values, nino3430.Time.values,
                          xlabel, ylabel, obsTitle+' (1976 - 2016)', calendar,
                          axis_font, linewidths, maxXTicks)

    # Plot Nino34 model time series
    plt.subplot(3, 1, 3)
    _plot_nino_timeseries(plt, nino34Index[2:-3].values,
                          nino34Index.Time[2:-3].values,
                          xlabel, ylabel, modelTitle, calendar,
                          axis_font, linewidths, maxXTicks)
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
    04/20/2017
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


def plot_polar_projection_comparison(
        config,
        x,
        y,
        landMask,
        modelArray,
        obsArray,
        diffArray,
        sectionName,
        fileout,
        title=None,
        modelTitle='Model',
        obsTitle='Observations',
        diffTitle='Model-Observations',
        cbarlabel='units',
        titleFontSize=None,
        figsize=(8, 22),
        dpi=300):

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
        model and observational data sets

    modelArray, obsArray : numpy ndarrays
        model and observational data sets

    diffArray : float array
        difference between modelArray and obsArray

    sectionName : str
        sectino name in ``config`` where color map info can be found.  The
        following options must be defined for suffixes ``Result`` and
        ``Difference``:
            ``colormapName<suffix>``, ``symLogNorm<suffix>``,
            ``colorbarTicks<suffix>``
        The colorbar for each panel will be constructed from these options

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
    Xylar Asay-Davis

    Last Modified
    -------------
    04/16/2017
    """

    def plot_panel(ax, title, array, cmap, norm, ticks):
        plt.title(title, y=1.06, **axis_font)

        plt.pcolormesh(x, y, array, cmap=cmap, norm=norm)
        cbar = plt.colorbar()
        cbar.set_label(cbarlabel)
        if ticks is not None:
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(['{}'.format(tick) for tick in ticks])
        plt.pcolormesh(x, y, landMask, cmap=landColorMap)
        plt.contour(xCenter, yCenter, landMask.mask, (0.5,), colors='k',
                    linewidths=0.5)
        ax.axis('off')
        ax.set_aspect('equal')
        ax.autoscale(tight=True)

    (cmapModelObs, normModelObs) = _setup_colormap_and_norm(
        config, sectionName, suffix='Result')
    (cmapDiff, normDiff) = _setup_colormap_and_norm(
        config, sectionName, suffix='Difference')

    colorbarTicksResult = config.getExpression(sectionName,
                                               'colorbarTicksResult')
    colorbarTicksDifference = config.getExpression(sectionName,
                                                   'colorbarTicksDifference')

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

    ax = plt.subplot(3, 1, 1)
    plot_panel(ax, modelTitle, modelArray, cmapModelObs, normModelObs,
               colorbarTicksResult)

    ax = plt.subplot(3, 1, 2)
    plot_panel(ax, obsTitle, obsArray, cmapModelObs, normModelObs,
               colorbarTicksResult)

    ax = plt.subplot(3, 1, 3)
    plot_panel(ax, diffTitle, diffArray, cmapDiff, normDiff,
               colorbarTicksDifference)

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    if not config.getboolean('plot', 'displayToScreen'):
        plt.close()


def plot_1D(config, xArrays, fieldArrays, errArrays,
            lineColors, lineWidths, legendText,
            title=None, xlabel=None, ylabel=None,
            fileout='plot_1D.png',
            figsize=(10, 4), dpi=300,
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
        number of dots per inch of the figure

    xLim : float array, optional
        x range of plot

    yLim : float array, optional
        y range of plot

    invertYAxis : logical, optional
        if True, invert Y axis

    Authors
    -------
    Mark Petersen, Milena Veneziani

    Last Modified
    -------------
    04/20/2017
    """

    # set up figure
    fig = plt.figure(figsize=figsize, dpi=dpi)

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
    colormapName,
    colorbarLevels,
    contourLevels,
    colorbarLabel=None,
    title=None,
    xlabel=None,
    ylabel=None,
    fileout='moc.png',
    figsize=(10, 4),
    dpi=300,
    xLim=None,
    yLim=None,
    invertYAxis=True):  # {{{

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

    xLim : float array, optional
        x range of plot

    yLim : float array, optional
        y range of plot

    invertYAxis : logical, optional
        if True, invert Y axis

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

    _register_custom_colormaps()

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

    Authors
    -------
    Xylar Asay-Davis

    Last modified
    -------------
    04/18/2017
    '''

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

    return (colormap, norm)


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
