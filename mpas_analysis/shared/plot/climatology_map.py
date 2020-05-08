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
Functions for plotting remapped horizontal fields (and comparing with reference
 data sets)
"""
# Authors
# -------
# Xylar Asay-Davis, Milena Veneziani, Luke Van Roekel, Greg Streletz

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import matplotlib.ticker as mticker
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy
from cartopy.util import add_cyclic_point

from mpas_analysis.shared.plot.colormap import setup_colormap


def plot_polar_comparison(
        config,
        lon,
        lat,
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
        model and observational or control run data sets

    diffArray : float array
        difference between modelArray and refArray

    colorMapSectionName : str
        section name in ``config`` where color map info can be found.

    fileout : str
        the file name to be written

    title : str, optional
        the subtitle of the plot

    plotProjection : 'npstere' or 'spstere', optional
        projection for the plot (north or south pole)

    modelTitle : str, optional
        title of the model panel

    refTitle : str, optional
        title of the observations or control run panel

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

        data_crs = cartopy.crs.PlateCarree()
        ax.set_extent(extent, crs=data_crs)

        ax.set_title(title, y=1.06, **plottitle_font)

        gl = ax.gridlines(crs=data_crs, color='k', linestyle=':', zorder=5,
                          draw_labels=True)
        gl.xlocator = mticker.FixedLocator(np.arange(-180., 181., 20.))
        gl.ylocator = mticker.FixedLocator(np.arange(-80., 81., 10.))
        gl.n_steps = 100
        gl.right_labels = False
        gl.xformatter = cartopy.mpl.gridliner.LONGITUDE_FORMATTER
        gl.yformatter = cartopy.mpl.gridliner.LATITUDE_FORMATTER

        fieldPeriodic, lonPeriodic = add_cyclic_point(field, lon)

        LonsPeriodic, LatsPeriodic = np.meshgrid(lonPeriodic, lat)

        if levels is None:
            plotHandle = ax.pcolormesh(LonsPeriodic, LatsPeriodic,
                                       fieldPeriodic, cmap=colormap,
                                       norm=norm, transform=data_crs,
                                       zorder=1)
        else:
            plotHandle = ax.contourf(LonsPeriodic, LatsPeriodic,
                                     fieldPeriodic, cmap=colormap,
                                     norm=norm, levels=levels,
                                     transform=data_crs,
                                     zorder=1)

        _add_land_lakes_coastline(ax)

        if contours is not None:
            matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
            ax.contour(LonsPeriodic, LatsPeriodic, fieldPeriodic,
                       levels=contours, colors=lineColor,
                       linewidths=lineWidth, transform=data_crs)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1,
                                  axes_class=plt.Axes)
        cbar = plt.colorbar(plotHandle, cax=cax, boundaries=levels)
        cbar.set_label(cbarlabel)
        if ticks is not None:
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(['{}'.format(tick) for tick in ticks])

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

    plottitle_font = {'size': config.get('plot',
                                         'threePanelPlotTitleFontSize')}

    if plotProjection == 'npstere':
        projection = cartopy.crs.NorthPolarStereo()
        extent = [-180, 180, latmin, 90]
    elif plotProjection == 'spstere':
        projection = cartopy.crs.SouthPolarStereo()
        extent = [-180, 180, -90, latmin]
    else:
        raise ValueError('Unexpected plot projection {}'.format(
                plotProjection))

    ax = plt.subplot(subplots[0], projection=projection)
    do_subplot(ax=ax, field=modelArray, title=modelTitle, **dictModelRef)

    if refArray is not None:
        ax = plt.subplot(subplots[1], projection=projection)
        do_subplot(ax=ax, field=refArray, title=refTitle, **dictModelRef)

        ax = plt.subplot(subplots[2], projection=projection)
        do_subplot(ax=ax, field=diffArray, title=diffTitle, **dictDiff)

    fig.canvas.draw()
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

    modelArray, refArray : float arrays
        model and observational or control run data sets

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
        title of the observations or control run panel

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

    def plot_panel(ax, title, array, colormap, norm, levels, ticks, contours,
                   lineWidth, lineColor):

        ax.set_extent(extent, crs=projection)

        ax.set_title(title, y=1.06, **plottitle_font)

        gl = ax.gridlines(crs=projection, color='k', linestyle=':', zorder=5,
                          draw_labels=True)
        gl.right_labels = False
        gl.top_labels = False
        gl.xlocator = mticker.FixedLocator(np.arange(-180., 181., 60.))
        gl.ylocator = mticker.FixedLocator(np.arange(-80., 81., 20.))
        gl.xformatter = cartopy.mpl.gridliner.LONGITUDE_FORMATTER
        gl.yformatter = cartopy.mpl.gridliner.LATITUDE_FORMATTER

        if levels is None:
            plotHandle = ax.pcolormesh(Lons, Lats, array, cmap=colormap,
                                       norm=norm, transform=projection,
                                       zorder=1)
        else:
            plotHandle = ax.contourf(Lons, Lats, array, cmap=colormap,
                                     norm=norm, levels=levels, extend='both',
                                     transform=projection, zorder=1)

        _add_land_lakes_coastline(ax)

        if contours is not None:
            matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
            ax.contour(Lons, Lats, array, levels=contours, colors=lineColor,
                       linewidths=lineWidth, transform=projection)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1,
                                  axes_class=plt.Axes)
        cbar = plt.colorbar(plotHandle, cax=cax, ticks=ticks, boundaries=ticks)
        cbar.set_label(cbarlabel)

    # set up figure
    if dpi is None:
        dpi = config.getint('plot', 'dpi')
    if figsize is None:
        # set the defaults, depending on if we have 1 or 3 panels
        if refArray is None:
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

    if refArray is None:
        subplots = [111]
    else:
        subplots = [311, 312, 313]

    projection = cartopy.crs.PlateCarree()

    extent = [-180, 180, -85, 85]

    dictModelRef = setup_colormap(config, colorMapSectionName, suffix='Result')
    dictDiff = setup_colormap(config, colorMapSectionName, suffix='Difference')

    ax = plt.subplot(subplots[0], projection=projection)
    plot_panel(ax, modelTitle, modelArray, **dictModelRef)

    if refArray is not None:
        ax = plt.subplot(subplots[1], projection=projection)
        plot_panel(ax, refTitle, refArray, **dictModelRef)

        ax = plt.subplot(subplots[2], projection=projection)
        plot_panel(ax, diffTitle, diffArray, **dictDiff)

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

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

    landMask : numpy ndarrays
        model and observational or control run data sets

    modelArray, refArray : numpy ndarrays
        model and observational or control run data sets

    diffArray : float array
        difference between modelArray and refArray

    fileout : str
        the file name to be written

    colorMapSectionName : str
        section name in ``config`` where color map info can be found.

    title : str, optional
        the subtitle of the plot

    modelTitle : str, optional
        title of the model panel

    refTitle : str, optional
        title of the observations or control run panel

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
    xCenter = 0.5 * (x[1:] + x[0:-1])
    yCenter = 0.5 * (y[1:] + y[0:-1])

    ax = plt.subplot(subplots[0])
    plot_panel(ax, modelTitle, modelArray, **dictModelRef)

    if refArray is not None:
        ax = plt.subplot(subplots[1])
        plot_panel(ax, refTitle, refArray, **dictModelRef)

        ax = plt.subplot(subplots[2])
        plot_panel(ax, diffTitle, diffArray, **dictDiff)

    if (fileout is not None):
        plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

    plt.close()


def _add_land_lakes_coastline(ax):
    land_50m = cartopy.feature.NaturalEarthFeature(
            'physical', 'land', '50m', edgecolor='k',
            facecolor='gray', linewidth=0.5)
    ice_50m = cartopy.feature.NaturalEarthFeature(
            'physical', 'antarctic_ice_shelves_polys', '50m', edgecolor='k',
            facecolor='lightgray', linewidth=0.5)
    lakes_50m = cartopy.feature.NaturalEarthFeature(
            'physical', 'lakes', '50m', edgecolor='k',
            facecolor='white',
            linewidth=0.5)
    ax.add_feature(land_50m, zorder=2)
    ax.add_feature(ice_50m, zorder=3)
    ax.add_feature(lakes_50m, zorder=4)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
