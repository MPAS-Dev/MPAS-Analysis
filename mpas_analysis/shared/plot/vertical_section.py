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
Funcitons for plotting vertical sections, both alone and as comparisons between
runs or with observations
"""
# Authors
# -------
# Xylar Asay-Davis, Milena Veneziani, Luke Van Roekel, Greg Streletz

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import xarray as xr
import numpy as np

from mpas_analysis.shared.timekeeping.utility import date_to_days

from mpas_analysis.shared.plot.colormap import setup_colormap
from mpas_analysis.shared.plot.ticks import plot_xtick_format
from mpas_analysis.shared.plot.title import limit_title


def plot_vertical_section_comparison(
        config,
        modelArray,
        refArray,
        diffArray,
        colorMapSectionName,
        xCoords=None,
        zCoord=None,
        triangulation_args=None,
        xOutlineModel=None,
        zOutlineModel=None,
        xOutlineRef=None,
        zOutlineRef=None,
        xOutlineDiff=None,
        zOutlineDiff=None,
        colorbarLabel=None,
        xlabels=None,
        ylabel=None,
        title=None,
        modelTitle='Model',
        refTitle='Observations',
        diffTitle='Model-Observations',
        titleFontSize=None,
        defaultFontSize=None,
        plotTitleFontSize=None,
        axisFontSize=None,
        figsize=None,
        dpi=None,
        lineWidth=2,
        lineStyle='solid',
        lineColor='black',
        contourColormap=None,
        backgroundColor='grey',
        invalidColor='white',
        outlineValid=True,
        xLim=None,
        yLim=None,
        numUpperTicks=None,
        upperXAxisTickLabelPrecision=None,
        invertYAxis=True,
        xCoordIsTime=False,
        movingAveragePoints=None,
        firstYearXTicks=None,
        yearStrideXTicks=None,
        maxXTicks=20,
        calendar='gregorian',
        compareAsContours=False,
        comparisonContourLineWidth=None,
        comparisonContourLineStyle=None,
        comparisonContourLineColor=None,
        labelContours=False,
        contourLabelPrecision=1,
        resultSuffix='Result',
        diffSuffix='Difference',
        maxTitleLength=None):
    """
    Plots vertical section plots in a three-panel format, comparing model data
    (in modelArray) to some reference dataset (in refArray), which can be
    either observations or an alternative model, and also presenting the
    difference plot of the two.  If refArray is None, then only one panel
    is plotted, displaying the model data.

    If compareAsContours is true, the contours of modelArray and refArray are
    plotted on a single plot.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    modelArray, refArray : xarray.DataArray
        model and observational or control run data sets

    diffArray : float array
        difference between modelArray and refArray

    xCoords : xarray.DataArray or list of xarray.DataArray, optional
       The x coordinate(s) for the model, ref and diff arrays.  Optional second
       and third entries will be used for a second and third x axis above the
       plot.  The typical use for the second and third axis is for transects,
       for which the primary x axis represents distance along a transect, and
       the second and third x axes are used to display the corresponding
       latitudes and longitudes.

    zCoord : xarray.DataArray, optional
        The z coordinates for the model, ref and diff arrays

    triangulation_args : dict, optional
        A dict of arguments to create a matplotlib.tri.Triangulation of the
        transect that does not rely on it being on a logically rectangular grid.
        The arguments rather than the triangulation itself are passed because
        multiple triangulations with different masks are needed internally and
        there is not an obvious mechanism for copying an existing triangulation.
        If this option is provided, ``xCoords`` is only used for tick marks if
        more than one x axis is requested, and ``zCoord`` will be ignored.

    xOutlineModel, zOutlineModel : numpy.ndarray, optional
        pairs of points defining line segments that are used to outline the
        valid region of the mesh for the model panel if ``outlineValid = True``
        and ``triangulation_args`` is not ``None``

    xOutlineRef, zOutlineRef : numpy.ndarray, optional
        Same as ``xOutlineModel`` and ``zOutlineModel`` but for the reference
        panel

    xOutlineDiff, zOutlineDiff : numpy.ndarray, optional
        Same as ``xOutlineModel`` and ``zOutlineModel`` but for the difference
        panel

    colorMapSectionName : str
        section name in ``config`` where color map info can be found.

    colorbarLabel : str, optional
        the label for the colorbar.  If compareAsContours and labelContours are
        both True, colorbarLabel is used as follows (typically in order to
        indicate the units that are associated with the contour labels):
        if refArray is None, the colorbarLabel string is parenthetically
        appended to the plot title;  if refArray is not None, it is
        parenthetically appended to the legend entries of the contour
        comparison plot.

    xlabels : str or list of str, optional
        labels of x-axes.  Labels correspond to entries in ``xCoords``.

    ylabel : str, optional
        label of y-axis

    title : str, optional
        the subtitle of the plot

    modelTitle : str, optional
        title of the model panel

    refTitle : str, optional
        title of the observations or control run panel

    diffTitle : str, optional
        title of the difference (bias) panel

    titleFontSize : int, optional
        size of the title font

    defaultFontSize : int, optional
        the size of text other than the title

    plotTitleFontSize : int, optional
        size of the title font for the individual plots

    axisFontSize : int, optional
        size of the axis font

    figsize : tuple of float, optional
        the size of the figure in inches

    dpi : int, optional
        the number of dots per inch of the figure, taken from section ``plot``
        option ``dpi`` in the config file by default

    lineWidth : float, optional
        the line width of contour lines (if specified)

    lineStyle : str, optional
        the line style of contour lines (if specified); this applies to the
        contour lines on heatmaps and to the contour lines of the model field
        on contour comparison plots (the line style of the contour lines of
        the reference field on contour comparison plots is set using the
        contourComparisonLineStyle argument).

    lineColor : str, optional
        the color of contour lines (if specified); this applies to the
        contour lines on heatmaps and to the contour lines of the model field
        on contour comparison plots (the line color of the contour lines of
        the reference field on contour comparison plots is set using the
        contourComparisonLineColor argument).

    backgroundColor : str, optional
        the background color for the plot outside the limits of ``xCoord`` and
        ``zCoord``.

    invalidColor : str, optional
        the color for invalid values (NaNs and masked areas will be
        shown in this color)

    outlineValid : bool, optional
        whether to outline the boundary between the valid an invalid regions
        with a black contour

    xLim : float array, optional
        x range of plot

    yLim : float array, optional
        y range of plot

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

    xCoordIsTime : logical, optional
        if True, format the x axis for time (this applies only to the primary
        x axis, not to the optional second or third x axes)

    movingAveragePoints : int, optional
        the number of points over which to perform a moving average
        NOTE: this option is mostly intended for use when xCoordIsTime is True,
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
        aspect ratio.  NOTE:  maxXTicks is only used if xCoordIsTime is True

    calendar : str, optional
        the calendar to use for formatting the time axis
        NOTE:  calendar is only used if xCoordIsTime is True

    compareAsContours : bool, optional
       if compareAsContours is True, instead of creating a three panel plot
       showing modelArray, refArray, and their difference, the function will
       plot the contours of modelArray and refArray on a single plot (unless
       refArray is None, in which case only the contours of modelArray will be
       plotted on the single panel plot).

    comparisonContourLineWidth : float, optional
        the line width of contour lines of the comparisonFieldName field on
        a contour comparison plot

    comparisonContourLineStyle : str, optional
        the line style of contour lines of the reference field on a contour
        comparison plot

    comparisonContourLineColor : str, optional
        the line color of contour lines of the reference field on a contour
        comparison plot

    labelContours : bool, optional
        whether or not to label contour lines (if specified) with their values

    contourLabelPrecision : int, optional
        the precision (in terms of number of figures to the right of the
        decimal point) of contour labels

    resultSuffix : str, optional
        a suffix added to the config options related to colormap information
        for the main and control fields

    diffSuffix : str, optional
        a suffix added to the config options related to colormap information
        for the difference field

    maxTitleLength : int or None, optional
        the maximum number of characters in the title, beyond which it is
        truncated with a trailing ellipsis.  The default is from the
        ``maxTitleLength`` config option.

    Returns
    -------
    fig : ``matplotlib.figure.Figure``
        The figure that was plotted

    axes : list of ``matplotlib.axes.Axes``
        The subplots

    suptitle : ``matplotlib.text.Text``
        The super-title
    """
    # Authors
    # -------
    # Greg Streletz, Xylar Asay-Davis, Milena Veneziani

    if maxTitleLength is None:
        maxTitleLength = config.getint('plot', 'maxTitleLength')

    if defaultFontSize is None:
        defaultFontSize = config.getint('plot', 'defaultFontSize')
    matplotlib.rc('font', size=defaultFontSize)
    if not isinstance(xCoords, list):
        xCoords = [xCoords]

    if not isinstance(xlabels, list):
        xlabels = [xlabels]

    if refArray is None or compareAsContours:
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
            if compareAsContours and refArray is not None and \
                    contourColormap is None:
                # no color bar but there is a legend at the bottom
                if len(xCoords) == 3:
                    figsize = (8, 8)
                else:
                    figsize = (8, 7)
            else:
                # color bar and legend
                figsize = (8, 7)
        elif len(xCoords) == 3:
            figsize = (8, 17)
        else:
            figsize = (8, 13)

    fig = plt.figure(figsize=figsize, dpi=dpi)

    if title is not None:
        if titleFontSize is None:
            titleFontSize = config.get('plot', 'threePanelTitleFontSize')
        title_font = {'size': titleFontSize,
                      'color': config.get('plot', 'threePanelTitleFontColor'),
                      'weight': config.get('plot',
                                           'threePanelTitleFontWeight')}
        suptitle = fig.suptitle(title, y=0.99, **title_font)
    else:
        suptitle = None

    if plotTitleFontSize is None:
        plotTitleFontSize = config.get('plot', 'threePanelPlotTitleFontSize')

    if len(xCoords) == 3:
        if singlePanel:
            titleY = 1.64
        else:
            titleY = 1.34
    elif len(xCoords) >= 2:
        titleY = 1.20
    else:
        titleY = 1.06

    if axisFontSize is None:
        axisFontSize = config.get('plot', 'threePanelAxisFontSize')

    if not singlePanel:
        plt.subplot(3, 1, 1)

    if not compareAsContours or refArray is None:
        title = modelTitle
        contourComparisonField = None
        comparisonFieldName = None
        originalFieldName = None
    else:
        title = None
        contourComparisonField = refArray
        comparisonFieldName = refTitle
        originalFieldName = modelTitle

    axes = []

    _, ax = plot_vertical_section(
        config,
        modelArray,
        colorMapSectionName,
        xCoords=xCoords,
        zCoord=zCoord,
        triangulation_args=triangulation_args,
        xOutline=xOutlineModel,
        zOutline=zOutlineModel,
        suffix=resultSuffix,
        colorbarLabel=colorbarLabel,
        title=title,
        xlabels=xlabels,
        ylabel=ylabel,
        figsize=None,
        titleFontSize=plotTitleFontSize,
        defaultFontSize=defaultFontSize,
        titleY=titleY,
        axisFontSize=axisFontSize,
        xLim=xLim,
        yLim=yLim,
        lineWidth=lineWidth,
        lineStyle=lineStyle,
        lineColor=lineColor,
        contourColormap=contourColormap,
        numUpperTicks=numUpperTicks,
        upperXAxisTickLabelPrecision=upperXAxisTickLabelPrecision,
        invertYAxis=invertYAxis,
        xCoordIsTime=xCoordIsTime,
        movingAveragePoints=movingAveragePoints,
        firstYearXTicks=firstYearXTicks,
        yearStrideXTicks=yearStrideXTicks,
        maxXTicks=maxXTicks, calendar=calendar,
        backgroundColor=backgroundColor,
        invalidColor=invalidColor,
        outlineValid=outlineValid,
        plotAsContours=compareAsContours,
        contourComparisonField=contourComparisonField,
        comparisonFieldName=comparisonFieldName,
        originalFieldName=originalFieldName,
        comparisonContourLineWidth=comparisonContourLineWidth,
        comparisonContourLineStyle=comparisonContourLineStyle,
        comparisonContourLineColor=comparisonContourLineColor,
        labelContours=labelContours,
        contourLabelPrecision=contourLabelPrecision,
        maxTitleLength=maxTitleLength)

    axes.append(ax)

    if not singlePanel:
        plt.subplot(3, 1, 2)
        _, ax = plot_vertical_section(
            config,
            refArray,
            colorMapSectionName,
            xCoords=xCoords,
            zCoord=zCoord,
            triangulation_args=triangulation_args,
            xOutline=xOutlineRef,
            zOutline=zOutlineRef,
            suffix=resultSuffix,
            colorbarLabel=colorbarLabel,
            title=refTitle,
            xlabels=xlabels,
            ylabel=ylabel,
            figsize=None,
            titleFontSize=plotTitleFontSize,
            defaultFontSize=defaultFontSize,
            titleY=titleY,
            axisFontSize=axisFontSize,
            xLim=xLim,
            yLim=yLim,
            lineWidth=lineWidth,
            lineStyle=lineStyle,
            lineColor=lineColor,
            upperXAxisTickLabelPrecision=upperXAxisTickLabelPrecision,
            numUpperTicks=numUpperTicks,
            invertYAxis=invertYAxis,
            xCoordIsTime=xCoordIsTime,
            movingAveragePoints=movingAveragePoints,
            firstYearXTicks=firstYearXTicks,
            yearStrideXTicks=yearStrideXTicks,
            maxXTicks=maxXTicks,
            calendar=calendar,
            backgroundColor=backgroundColor,
            invalidColor=invalidColor,
            outlineValid=outlineValid,
            labelContours=labelContours,
            contourLabelPrecision=contourLabelPrecision,
            maxTitleLength=maxTitleLength)

        axes.append(ax)

        plt.subplot(3, 1, 3)
        _, ax = plot_vertical_section(
            config,
            diffArray,
            colorMapSectionName,
            xCoords=xCoords,
            zCoord=zCoord,
            triangulation_args=triangulation_args,
            xOutline=xOutlineDiff,
            zOutline=zOutlineDiff,
            suffix=diffSuffix,
            colorbarLabel=colorbarLabel,
            title=diffTitle,
            xlabels=xlabels,
            ylabel=ylabel,
            figsize=None,
            titleFontSize=plotTitleFontSize,
            defaultFontSize=defaultFontSize,
            titleY=titleY,
            axisFontSize=axisFontSize,
            xLim=xLim,
            yLim=yLim,
            lineWidth=lineWidth,
            lineStyle=lineStyle,
            lineColor=lineColor,
            upperXAxisTickLabelPrecision=upperXAxisTickLabelPrecision,
            numUpperTicks=numUpperTicks,
            invertYAxis=invertYAxis,
            xCoordIsTime=xCoordIsTime,
            movingAveragePoints=movingAveragePoints,
            firstYearXTicks=firstYearXTicks,
            yearStrideXTicks=yearStrideXTicks,
            maxXTicks=maxXTicks,
            calendar=calendar,
            backgroundColor=backgroundColor,
            invalidColor=invalidColor,
            outlineValid=outlineValid,
            labelContours=labelContours,
            contourLabelPrecision=contourLabelPrecision,
            maxTitleLength=maxTitleLength)

        axes.append(ax)

    if singlePanel:
        if len(xCoords) == 3 and refArray is None:
            plt.tight_layout(pad=0.0, h_pad=2.0, rect=[0.0, 0.0, 1.0, 0.98])
        else:
            plt.tight_layout(pad=0.0, h_pad=2.0, rect=[0.0, 0.0, 1.0, 0.95])
    else:
        plt.tight_layout(pad=0.0, h_pad=2.0, rect=[0.01, 0.0, 1.0, 0.97])

    return fig, axes, suptitle


def plot_vertical_section(
        config,
        field,
        colorMapSectionName,
        xCoords=None,
        zCoord=None,
        triangulation_args=None,
        xOutline=None,
        zOutline=None,
        suffix='',
        colorbarLabel=None,
        title=None,
        xlabels=None,
        ylabel=None,
        figsize=(10, 4),
        dpi=None,
        titleFontSize=None,
        defaultFontSize=None,
        titleY=None,
        axisFontSize=None,
        xLim=None,
        yLim=None,
        lineWidth=2,
        lineStyle='solid',
        lineColor='black',
        contourColormap=None,
        backgroundColor='grey',
        invalidColor='white',
        outlineValid=True,
        numUpperTicks=None,
        upperXAxisTickLabelPrecision=None,
        invertYAxis=True,
        xCoordIsTime=False,
        movingAveragePoints=None,
        firstYearXTicks=None,
        yearStrideXTicks=None,
        maxXTicks=20,
        calendar='gregorian',
        plotAsContours=False,
        contourComparisonField=None,
        comparisonFieldName=None,
        originalFieldName=None,
        comparisonContourLineWidth=None,
        comparisonContourLineStyle=None,
        comparisonContourLineColor=None,
        labelContours=False,
        contourLabelPrecision=1,
        maxTitleLength=None):
    """
    Plots a data set as a x distance (latitude, longitude,
    or spherical distance) vs depth map (vertical section).

    Or, if xCoordIsTime is True, plots data set on a vertical
    Hovmoller plot (depth vs. time).

    Typically, the ``field`` data are plotted using a heatmap, but if
    ``contourComparisonField`` is not None, then contours of both
    ``field`` and ``contourComparisonField`` are plotted instead.

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    field : xarray.DataArray
        field array to plot.  For contour plots, ``xCoords`` and ``zCoords``
        should broadcast to the same shape as ``field``.  For heatmap plots,
        ``xCoords`` and ``zCoords`` are the corners of the plot.  If they
        broadcast to the same shape as ``field``, ``field`` will be bilinearly
        interpolated to center values for each plot cell.  If the coordinates
        have one extra element in each direction than ``field``, ``field`` is
        assumed to contain cell values and no interpolation is performed.

    colorMapSectionName : str
        section name in ``config`` where color map info can be found.

    xCoords : xarray.DataArray or list of xarray.DataArray, optional
       The x coordinate(s) for the ``field``.  Optional second
       and third entries will be used for a second and third x axis above the
       plot.  The typical use for the second and third axis is for transects,
       for which the primary x axis represents distance along a transect, and
       the second and third x axes are used to display the corresponding
       latitudes and longitudes.

    zCoord : xarray.DataArray, optional
        The z coordinates for the ``field``

    triangulation_args : dict, optional
        A dict of arguments to create a matplotlib.tri.Triangulation of the
        transect that does not rely on it being on a logically rectangular grid.
        The arguments rather than the triangulation itself are passed because
        multiple triangulations with different masks are needed internally and
        there is not an obvious mechanism for copying an existing triangulation.
        If this option is provided, ``xCoords`` is only used for tick marks if
        more than one x axis is requested, and ``zCoord`` will be ignored.

    xOutline, zOutline : numpy.ndarray, optional
        pairs of points defining line segments that are used to outline the
        valid region of the mesh if ``outlineValid = True`` and
        ``triangulation_args`` is not ``None``



    suffix : str, optional
        the suffix used for colorbar config options

    colorbarLabel : str, optional
        the label for the colorbar.  If plotAsContours and labelContours are
        both True, colorbarLabel is used as follows (typically in order to
        indicate the units that are associated with the contour labels):
        if ``contourComparisonField`` is None, the ``colorbarLabel`` string is
        parenthetically appended to the plot title;  if
        ``contourComparisonField`` is not None, it is parenthetically appended
        to the legend entries of the contour comparison plot.

    title : str, optional
        title of plot

    xlabels : str or list of str, optional
        labels of x-axes.  Labels correspond to entries in ``xCoords``.

    ylabel : str, optional
        label of y-axis

    figsize : tuple of float, optional
        size of the figure in inches, or None if the current figure should
        be used (e.g. if this is a subplot)

    dpi : int, optional
        the number of dots per inch of the figure, taken from section ``plot``
        option ``dpi`` in the config file by default

    titleFontSize : int, optional
        size of the title font

    defaultFontSize : int, optional
        the size of text other than the title

    titleY : float, optional
        the y value to use for placing the plot title

    axisFontSize : int, optional
        size of the axis font

    xLim : float array, optional
        x range of plot

    yLim : float array, optional
        y range of plot

    lineWidth : float, optional
        the line width of contour lines (if specified)

    lineStyle : str, optional
        the line style of contour lines (if specified); this applies to the
        style of contour lines of fieldArray (the style of the contour lines
        of contourComparisonField is set using
        contourComparisonLineStyle).

    lineColor : str, optional
        the color of contour lines (if specified); this applies to the
        contour lines of fieldArray (the color of the contour lines of
        contourComparisonField is set using contourComparisonLineColor

    backgroundColor : str, optional
        the background color for the plot outside the limits of ``xCoord`` and
        ``zCoord``.

    invalidColor : str, optional
        the color for invalid values (NaNs and masked areas will be
        shown in this color)

    outlineValid : bool, optional
        whether to outline the boundary between the valid an invalid regions
        with a black contour

    numUpperTicks :  int, optional
        the approximate number of ticks to use on the upper x axis
        or axes (these are the second and third x axes, which are placed above
        the plot if they have been requested by specifying the secondXAxisData
        or thirdXAxisData arrays above)

    upperXAxisTickLabelPrecision : int, optional
        the number of decimal places (to the right
        of the decimal point) to use for values at upper axis ticks.  This
        value can be adjusted (in concert with numUpperTicks) to avoid problems
        with overlapping numbers along the upper axis.

    invertYAxis : logical, optional
        if True, invert Y axis

    xCoordIsTime : logical, optional
        if True, format the x axis for time (this applies only to the primary
        x axis, not to the optional second or third x axes)

    movingAveragePoints : int, optional
        the number of points over which to perform a moving average
        NOTE: this option is mostly intended for use when ``xCoordIsTime`` is
        True, although it will work with other data as well.  Also, the moving
        average calculation is based on number of points, not actual x axis
        values, so for best results, the values in the first entry in
        ``xCoords`` should be equally spaced.

    firstYearXTicks : int, optional
        The year of the first tick on the x axis.  By default, the first time
        entry is the first tick.

    yearStrideXTicks : int, optional
        The number of years between x ticks. By default, the stride is chosen
        automatically to have ``maxXTicks`` tick marks or fewer.

    maxXTicks : int, optional
        the maximum number of tick marks that will be allowed along the primary
        x axis.  This may need to be adjusted depending on the figure size and
        aspect ratio.  NOTE:  maxXTicks is only used if xCoordIsTime is True

    calendar : str, optional
        the calendar to use for formatting the time axis
        NOTE:  calendar is only used if xCoordIsTime is True

    plotAsContours : bool, optional
        if plotAsContours is True, instead of plotting ``field`` as a
        heatmap, the function will plot only the contours of ``field``.  In
        addition, if contourComparisonField is not None, the contours
        of this field will be plotted on the same plot.  The selection of
        contour levels is still determined as for the contours on the heatmap
        plots, via the 'contours' entry in ``colorMapSectionName``.

    contourComparisonField : float array, optional
        a comparison ``field`` array (typically observational data or results
        from another simulation run), assumed to be of the same shape as
        ``field``. If ``plotAsContours`` is ``True`` and
        ``countourComparisonFieldArray`` is not ``None``, then contours of both
        ``field`` and ``contourComparisonField`` will be plotted in order to
        enable a comparison of the two fields on the same plot.

    comparisonFieldName : str, optional
        the name for the comparison field.  If contourComparisonField is
        None, this parameter is ignored.

    originalFieldName : str, optional
        the name for the ``field`` field (for the purposes of labeling the
        contours on a contour comparison plot).  If contourComparisonField
        is None, this parameter is ignored.

    comparisonContourLineWidth : float, optional
        the line width of contour lines of the comparisonFieldName field on
        a contour comparison plot

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

    maxTitleLength : int or None, optional
        the maximum number of characters in the title, beyond which it is
        truncated with a trailing ellipsis.  The default is from the
        ``maxTitleLength`` config option.

    Returns
    -------
    fig : ``matplotlib.figure.Figure``
        The figure that was plotted

    ax : ``matplotlib.axes.Axes``
        The subplot
    """
    # Authors
    # -------
    # Milena Veneziani, Mark Petersen, Xylar Asay-Davis, Greg Streletz

    if maxTitleLength is None:
        maxTitleLength = config.getint('plot', 'maxTitleLength')

    if defaultFontSize is None:
        defaultFontSize = config.getint('plot', 'defaultFontSize')
    matplotlib.rc('font', size=defaultFontSize)
    if xCoords is not None:
        if not isinstance(xCoords, list):
            xCoords = [xCoords]

        if not isinstance(xlabels, list):
            xlabels = [xlabels]

        if len(xCoords) != len(xlabels):
            raise ValueError('Expected the same number of xCoords and xlabels')

    if triangulation_args is None:

        x, y = xr.broadcast(xCoords[0], zCoord)
        dims_in_field = all([dim in field.dims for dim in x.dims])

        if dims_in_field:
            x = x.transpose(*field.dims)
            y = y.transpose(*field.dims)
        else:
            xsize = list(x.sizes.values())
            fieldsize = list(field.sizes.values())
            if xsize[0] == fieldsize[0] + 1 and xsize[1] == fieldsize[1] + 1:
                pass
            elif xsize[0] == fieldsize[1] + 1 and xsize[1] == fieldsize[0] + 1:
                x = x.transpose(x.dims[1], x.dims[0])
                y = y.transpose(y.dims[1], y.dims[0])
            else:
                raise ValueError('Sizes of coords {}x{} and field {}x{} not '
                                 'compatible.'.format(xsize[0], xsize[1],
                                                      fieldsize[0],
                                                      fieldsize[1]))

        # compute moving averages with respect to the x dimension
        if movingAveragePoints is not None and movingAveragePoints != 1:
            dim = field.dims[0]
            field = field.rolling(dim={dim: movingAveragePoints},
                                  center=True).mean().dropna(dim, how='all')
            x = x.rolling(dim={dim: movingAveragePoints},
                          center=True).mean().dropna(dim, how='all')
            y = y.rolling(dim={dim: movingAveragePoints},
                          center=True).mean().dropna(dim, how='all')

        mask = field.notnull()
        maskedTriangulation, unmaskedTriangulation = _get_triangulation(
            x, y, mask)
        if contourComparisonField is not None:
            mask = field.notnull()
            maskedComparisonTriangulation, _ = _get_triangulation(x, y, mask)
        else:
            maskedComparisonTriangulation = None
    else:
        mask = field.notnull()
        triMask = np.logical_not(mask.values)
        # if any node of a triangle is masked, the triangle is masked
        triMask = np.amax(triMask, axis=1)
        unmaskedTriangulation = Triangulation(**triangulation_args)
        anythingToPlot = not np.all(triMask)
        if anythingToPlot:
            mask_args = dict(triangulation_args)
            mask_args['mask'] = triMask
            maskedTriangulation = Triangulation(**mask_args)
        else:
            maskedTriangulation = None
        if contourComparisonField is not None:
            mask = contourComparisonField.notnull()
            triMask = np.logical_not(mask.values)
            triMask = np.amax(triMask, axis=1)
            mask_args = dict(triangulation_args)
            mask_args['mask'] = triMask
            maskedComparisonTriangulation = Triangulation(**mask_args)
        else:
            maskedComparisonTriangulation = None

    # set up figure
    if dpi is None:
        dpi = config.getint('plot', 'dpi')
    if figsize is not None:
        fig = plt.figure(figsize=figsize, dpi=dpi)
    else:
        fig = plt.gcf()

    colormapDict = setup_colormap(config, colorMapSectionName,
                                  suffix=suffix)

    # fill the unmasked region with the invalid color so it will show through
    # any masked regions
    zeroArray = xr.zeros_like(field)
    plt.tricontourf(unmaskedTriangulation, zeroArray.values.ravel(),
                    colors=invalidColor)

    if maskedTriangulation is not None:
        # there's something to plot
        if not plotAsContours:
            # display a heatmap of fieldArray
            fieldMasked = field.where(mask, 0.0).values.ravel()

            if colormapDict['levels'] is None:

                plotHandle = plt.tripcolor(maskedTriangulation, fieldMasked,
                                           cmap=colormapDict['colormap'],
                                           norm=colormapDict['norm'],
                                           rasterized=True, shading='gouraud')
            else:
                plotHandle = plt.tricontourf(maskedTriangulation, fieldMasked,
                                             cmap=colormapDict['colormap'],
                                             norm=colormapDict['norm'],
                                             levels=colormapDict['levels'],
                                             extend='both')

            cbar = plt.colorbar(plotHandle,
                                orientation='vertical',
                                spacing='uniform',
                                aspect=9,
                                ticks=colormapDict['ticks'])

            if colorbarLabel is not None:
                cbar.set_label(colorbarLabel)

        else:
            # display a white heatmap to get a white background for non-land
            zeroArray = xr.zeros_like(field)
            plt.tricontourf(maskedTriangulation, zeroArray.values.ravel(),
                            colors='white')

    ax = plt.gca()
    ax.set_facecolor(backgroundColor)
    if outlineValid:
        if xOutline is not None and zOutline is not None:
            # also outline the domain if provided
            plt.plot(xOutline, zOutline, color='black', linewidth=1)
        else:
            # do a contour to outline the boundary between valid and invalid
            # values
            landMask = np.isnan(field.values).ravel()
            plt.tricontour(unmaskedTriangulation, landMask, levels=[0.0001],
                           colors='black', linewidths=1)

    # plot contours, if they were requested
    contourLevels = colormapDict['contours']
    fmt_string = None
    cs1 = None
    cs2 = None
    plotLegend = False

    if contourLevels is not None and maskedTriangulation is not None:
        if len(contourLevels) == 0:
            # automatic calculation of contour levels
            contourLevels = None
        mask = field.notnull()
        fieldMasked = field.where(mask, 0.0).values.ravel()

        cs1 = plt.tricontour(maskedTriangulation, fieldMasked,
                             levels=contourLevels,
                             colors=lineColor,
                             linestyles=lineStyle,
                             linewidths=lineWidth,
                             cmap=contourColormap)
        if labelContours:
            fmt_string = "%%1.%df" % int(contourLabelPrecision)
            plt.clabel(cs1, fmt=fmt_string)

        if plotAsContours and contourComparisonField is not None:
            if comparisonContourLineWidth is None:
                comparisonContourLineWidth = lineWidth
            mask = contourComparisonField.notnull()
            fieldMasked = contourComparisonField.where(mask, 0.0).values.ravel()
            cs2 = plt.tricontour(maskedComparisonTriangulation,
                                 fieldMasked,
                                 levels=contourLevels,
                                 colors=comparisonContourLineColor,
                                 linestyles=comparisonContourLineStyle,
                                 linewidths=comparisonContourLineWidth,
                                 cmap=contourColormap)

            if labelContours:
                plt.clabel(cs2, fmt=fmt_string)

        plotLegend = (((lineColor is not None and
                        comparisonContourLineColor is not None) or
                    (lineWidth is not None and
                        comparisonContourLineWidth is not None)) and
                    (plotAsContours and contourComparisonField is not None))

    if plotLegend:
        h1, _ = cs1.legend_elements()
        h2, _ = cs2.legend_elements()
        if labelContours:
            originalFieldName = originalFieldName + " (" + colorbarLabel + ")"
            comparisonFieldName = (comparisonFieldName + " (" +
                                   colorbarLabel + ")")
        ax.legend([h1[0], h2[0]], [originalFieldName, comparisonFieldName],
                  loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=1)

    if title is not None:
        if plotAsContours and labelContours \
                and contourComparisonField is None:
            title = limit_title(title,
                                maxTitleLength - (3 + len(colorbarLabel)))
            title = title + " (" + colorbarLabel + ")"
        else:
            title = limit_title(title, maxTitleLength)
        if titleFontSize is None:
            titleFontSize = config.get('plot', 'titleFontSize')
        title_font = {'size': titleFontSize,
                      'color': config.get('plot', 'titleFontColor'),
                      'weight': config.get('plot', 'titleFontWeight')}
        if titleY is not None:
            plt.title(title, y=titleY, **title_font)
        else:
            plt.title(title, **title_font)

    if axisFontSize is None:
        axisFontSize = config.get('plot', 'axisFontSize')
    axis_font = {'size': axisFontSize}

    if xlabels is not None:
        plt.xlabel(xlabels[0], **axis_font)
    if ylabel is not None:
        plt.ylabel(ylabel, **axis_font)

    if invertYAxis:
        ax.invert_yaxis()

    if xLim:
        ax.set_xlim(xLim)
    if yLim:
        ax.set_ylim(yLim)

    if xCoords is not None and xCoordIsTime:
        if firstYearXTicks is None:
            minDays = xCoords[0][0].values
        else:
            minDays = date_to_days(year=firstYearXTicks, calendar=calendar)
        maxDays = xCoords[0][-1].values

        plot_xtick_format(calendar, minDays, maxDays, maxXTicks,
                          yearStride=yearStrideXTicks)

    if contourLevels is not None:
        if contourColormap is not None:
            cbar1 = fig.colorbar(cs1, ax=ax, fraction=.05,
                                 orientation='vertical',
                                 spacing='proportional')

            if colorbarLabel is not None:
                cbar1.set_label(colorbarLabel)

    # add a second x-axis scale, if it was requested
    if xCoords is not None and len(xCoords) >= 2:
        ax2 = ax.twiny()
        ax2.set_facecolor(backgroundColor)
        if xlabels[1] is not None:
            ax2.set_xlabel(xlabels[1], **axis_font)
        xlimits = ax.get_xlim()
        ax2.set_xlim(xlimits)
        formatString = None
        xticks = None
        if numUpperTicks is not None:
            xticks = np.linspace(xlimits[0], xlimits[1], numUpperTicks)
            tickValues = np.interp(xticks, xCoords[0].values, xCoords[1].values)
            ax2.set_xticks(xticks)
            formatString = "{{0:.{:d}f}}{}".format(
                upperXAxisTickLabelPrecision, r'$\degree$')
            ax2.set_xticklabels([formatString.format(member)
                                 for member in tickValues])

        # add a third x-axis scale, if it was requested
        if len(xCoords) == 3:
            ax3 = ax.twiny()
            ax3.set_facecolor(backgroundColor)
            ax3.set_xlabel(xlabels[2], **axis_font)
            ax3.set_xlim(xlimits)
            ax3.set_xticks(xticks)
            if numUpperTicks is not None:
                tickValues = np.interp(xticks, xCoords[0].values,
                                       xCoords[2].values)
                ax3.set_xticklabels([formatString.format(member)
                                     for member in tickValues])
                ax3.spines['top'].set_position(('outward', 36))

    return fig, ax


def _get_triangulation(x, y, mask):
    """divide each quad in the x/y mesh into 2 triangles"""

    nx = x.sizes[x.dims[0]] - 1
    ny = x.sizes[x.dims[1]] - 1
    nTriangles = 2 * nx * ny

    mask = mask.values
    mask = np.logical_and(np.logical_and(mask[0:-1, 0:-1], mask[1:, 0:-1]),
                          np.logical_and(mask[0:-1, 1:], mask[1:, 1:]))
    triMask = np.zeros((nx, ny, 2), bool)
    triMask[:, :, 0] = np.logical_not(mask)
    triMask[:, :, 1] = triMask[:, :, 0]

    triMask = triMask.ravel()

    xIndices, yIndices = np.meshgrid(np.arange(nx), np.arange(ny),
                                     indexing='ij')

    tris = np.zeros((nx, ny, 2, 3), int)
    # upper triangles:
    tris[:, :, 0, 0] = (ny + 1) * xIndices + yIndices
    tris[:, :, 0, 1] = (ny + 1) * (xIndices + 1) + yIndices
    tris[:, :, 0, 2] = (ny + 1) * xIndices + yIndices + 1
    # lower triangle
    tris[:, :, 1, 0] = (ny + 1) * xIndices + yIndices + 1
    tris[:, :, 1, 1] = (ny + 1) * (xIndices + 1) + yIndices
    tris[:, :, 1, 2] = (ny + 1) * (xIndices + 1) + yIndices + 1

    tris = tris.reshape((nTriangles, 3))

    x = x.values.ravel()
    y = y.values.ravel()

    anythingToPlot = not np.all(triMask)
    if anythingToPlot:
        maskedTriangulation = Triangulation(x=x, y=y, triangles=tris,
                                            mask=triMask)
    else:
        maskedTriangulation = None

    unmaskedTriangulation = Triangulation(x=x, y=y, triangles=tris)

    return maskedTriangulation, unmaskedTriangulation
