# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
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
from mpl_toolkits.basemap import Basemap
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

from mpas_analysis.shared.plot.colormap import setup_colormap


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
        model and observational or control run data sets

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


# function copied from basemap because the latest version of this funciton
# is buggy but no new release seems imminent
def addcyclic(*arr, **kwargs):
    """
    Adds cyclic (wraparound) points in longitude to one or several arrays,
    the last array being longitudes in degrees. e.g.
   ``data1out, data2out, lonsout = addcyclic(data1,data2,lons)``
    ==============   ====================================================
    Keywords         Description
    ==============   ====================================================
    axis             the dimension representing longitude (default -1,
                     or right-most)
    cyclic           width of periodic domain (default 360)
    ==============   ====================================================
    """
    # get (default) keyword arguments
    axis = kwargs.get('axis', -1)
    cyclic = kwargs.get('cyclic', 360)
    # define functions

    def _addcyclic(a):
        """addcyclic function for a single data array"""
        npsel = np.ma if np.ma.is_masked(a) else np
        slicer = [slice(None)] * np.ndim(a)
        try:
            slicer[axis] = slice(0, 1)
        except IndexError:
            raise ValueError('The specified axis does not correspond to an '
                             'array dimension.')
        return npsel.concatenate((a, a[tuple(slicer)]), axis=axis)

    def _addcyclic_lon(a):
        """addcyclic function for a single longitude array"""
        # select the right numpy functions
        npsel = np.ma if np.ma.is_masked(a) else np
        # get cyclic longitudes
        clon = (np.take(a, [0], axis=axis)
                + cyclic * np.sign(np.diff(np.take(a, [0, -1], axis=axis),
                                           axis=axis)))
        # ensure the values do not exceed cyclic
        clonmod = npsel.where(clon <= cyclic, clon, np.mod(clon, cyclic))
        return npsel.concatenate((a, clonmod), axis=axis)
    # process array(s)
    if len(arr) == 1:
        return _addcyclic_lon(arr[-1])
    else:
        return list(map(_addcyclic, arr[:-1])) + [_addcyclic_lon(arr[-1])]


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
