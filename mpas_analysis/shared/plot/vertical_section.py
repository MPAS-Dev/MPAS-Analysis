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
Funcitons for plotting vertical sections, both alone and as comparisons between
runs or with observations
"""
# Authors
# -------
# Xylar Asay-Davis, Milena Veneziani, Luke Van Roekel, Greg Streletz

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import matplotlib
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np

from mpas_analysis.shared.timekeeping.utility import date_to_days

from mpas_analysis.shared.plot.colormap import setup_colormap
from mpas_analysis.shared.plot.ticks import plot_xtick_format
from mpas_analysis.shared.plot.title import limit_title


def plot_vertical_section_comparison(
        config,
        xArray,
        depthArray,
        modelArray,
        refArray,
        diffArray,
        colorMapSectionName,
        colorbarLabel=None,
        xlabel=None,
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
        movingAveragePoints=None,
        firstYearXTicks=None,
        yearStrideXTicks=None,
        maxXTicks=20,
        calendar='gregorian',
        compareAsContours=False,
        comparisonContourLineStyle=None,
        comparisonContourLineColor=None,
        labelContours=False,
        contourLabelPrecision=1,
        resultSuffix='Result',
        diffSuffix='Difference',
        maxTitleLength=70):
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

    xArray : float array
        x array (latitude, longitude, spherical distance, or distance along
        a transect;  or, time for Hovmoller plots)

    depthArray : float array
        depth array [m]

    modelArray, refArray : float arrays
        model and observational or control run data sets

    diffArray : float array
        difference between modelArray and refArray

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

    xlabel, ylabel : str, optional
        label of x- and y-axis

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

    lineWidth : int, optional
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

    movingAveragePoints : int, optional
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
       showing modelArray, refArray, and their difference, the function will
       plot the contours of modelArray and refArray on a single plot (unless
       refArray is None, in which case only the contours of modelArray will be
       plotted on the single panel plot).

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

    maxTitleLength : int, optional
        the maximum number of characters in the title, beyond which it is
        truncated with a trailing ellipsis

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

    if defaultFontSize is None:
        defaultFontSize = config.getint('plot', 'defaultFontSize')
    matplotlib.rc('font', size=defaultFontSize)

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
            if compareAsContours and refArray is not None:
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

    if not compareAsContours or refArray is None:
        title = modelTitle
        contourComparisonFieldArray = None
        comparisonFieldName = None
        originalFieldName = None
    else:
        title = None
        contourComparisonFieldArray = refArray
        comparisonFieldName = refTitle
        originalFieldName = modelTitle

    axes = []

    _, ax = plot_vertical_section(
        config,
        xArray,
        depthArray,
        modelArray,
        colorMapSectionName,
        suffix=resultSuffix,
        colorbarLabel=colorbarLabel,
        title=title,
        xlabel=xlabel,
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
        secondXAxisData=secondXAxisData,
        secondXAxisLabel=secondXAxisLabel,
        thirdXAxisData=thirdXAxisData,
        thirdXAxisLabel=thirdXAxisLabel,
        numUpperTicks=numUpperTicks,
        upperXAxisTickLabelPrecision=upperXAxisTickLabelPrecision,
        invertYAxis=invertYAxis,
        xArrayIsTime=xArrayIsTime,
        movingAveragePoints=movingAveragePoints,
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
        contourLabelPrecision=contourLabelPrecision,
        maxTitleLength=maxTitleLength)

    axes.append(ax)

    if not singlePanel:
        plt.subplot(3, 1, 2)
        _, ax = plot_vertical_section(
            config,
            xArray,
            depthArray,
            refArray,
            colorMapSectionName,
            suffix=resultSuffix,
            colorbarLabel=colorbarLabel,
            title=refTitle,
            xlabel=xlabel,
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
            secondXAxisData=secondXAxisData,
            secondXAxisLabel=secondXAxisLabel,
            thirdXAxisData=thirdXAxisData,
            thirdXAxisLabel=thirdXAxisLabel,
            upperXAxisTickLabelPrecision=upperXAxisTickLabelPrecision,
            numUpperTicks=numUpperTicks,
            invertYAxis=invertYAxis,
            xArrayIsTime=xArrayIsTime,
            movingAveragePoints=movingAveragePoints,
            firstYearXTicks=firstYearXTicks,
            yearStrideXTicks=yearStrideXTicks,
            maxXTicks=maxXTicks,
            calendar=calendar,
            backgroundColor=backgroundColor,
            labelContours=labelContours,
            contourLabelPrecision=contourLabelPrecision,
            maxTitleLength=maxTitleLength)

        axes.append(ax)

        plt.subplot(3, 1, 3)
        _, ax = plot_vertical_section(
            config,
            xArray,
            depthArray,
            diffArray,
            colorMapSectionName,
            suffix=diffSuffix,
            colorbarLabel=colorbarLabel,
            title=diffTitle,
            xlabel=xlabel,
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
            secondXAxisData=secondXAxisData,
            secondXAxisLabel=secondXAxisLabel,
            thirdXAxisData=thirdXAxisData,
            thirdXAxisLabel=thirdXAxisLabel,
            upperXAxisTickLabelPrecision=upperXAxisTickLabelPrecision,
            numUpperTicks=numUpperTicks,
            invertYAxis=invertYAxis,
            xArrayIsTime=xArrayIsTime,
            movingAveragePoints=movingAveragePoints,
            firstYearXTicks=firstYearXTicks,
            yearStrideXTicks=yearStrideXTicks,
            maxXTicks=maxXTicks,
            calendar=calendar,
            backgroundColor=backgroundColor,
            labelContours=labelContours,
            contourLabelPrecision=contourLabelPrecision,
            maxTitleLength=maxTitleLength)

        axes.append(ax)

    if singlePanel:
        if thirdXAxisData is not None and refArray is None:
            plt.tight_layout(pad=0.0, h_pad=2.0, rect=[0.0, 0.0, 1.0, 0.98])
        else:
            plt.tight_layout(pad=0.0, h_pad=2.0, rect=[0.0, 0.0, 1.0, 0.95])
    else:
        plt.tight_layout(pad=0.0, h_pad=2.0, rect=[0.01, 0.0, 1.0, 0.97])

    return fig, axes, suptitle


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
        backgroundColor='grey',
        secondXAxisData=None,
        secondXAxisLabel=None,
        thirdXAxisData=None,
        thirdXAxisLabel=None,
        numUpperTicks=None,
        upperXAxisTickLabelPrecision=None,
        invertYAxis=True,
        xArrayIsTime=False,
        movingAveragePoints=None,
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
        contourLabelPrecision=1,
        maxTitleLength=70):  # {{{
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

    movingAveragePoints : int, optional
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

    maxTitleLength : int, optional
        the maximum number of characters in the title, beyond which it is
        truncated with a trailing ellipsis

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

    if defaultFontSize is None:
        defaultFontSize = config.getint('plot', 'defaultFontSize')
    matplotlib.rc('font', size=defaultFontSize)

    # compute moving averages with respect to the x dimension
    if movingAveragePoints is not None and movingAveragePoints != 1:
        N = movingAveragePoints
        movingAverageDepthSlices = []
        for nVertLevel in range(len(depthArray)):
            depthSlice = fieldArray[[nVertLevel]][0]
            # in case it's not an xarray already
            depthSlice = xr.DataArray(depthSlice)
            mean = pd.Series.rolling(depthSlice.to_series(), N,
                                     center=True).mean()
            mean = xr.DataArray.from_series(mean)
            mean = mean[int(N / 2.0):-int(round(N / 2.0) - 1)]
            movingAverageDepthSlices.append(mean)
        xArray = xArray[int(N / 2.0):-int(round(N / 2.0) - 1)]
        fieldArray = xr.DataArray(movingAverageDepthSlices)

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
    if secondXAxisData is None and thirdXAxisData is not None:
        raise ValueError('secondXAxisData cannot be None if thirdXAxisData '
                         'is not None')
    if secondXAxisData is not None:
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
    if thirdXAxisData is not None:
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
    if figsize is not None:
        fig = plt.figure(figsize=figsize, dpi=dpi)
    else:
        fig = plt.gcf()

    colormapDict = setup_colormap(config, colorMapSectionName, suffix=suffix)

    if not plotAsContours:    # display a heatmap of fieldArray

        if colormapDict['levels'] is None:
            # interpFieldArray contains the values at centers of grid cells,
            # for pcolormesh plots (using bilinear interpolation)
            interpFieldArray = \
                0.5 * (0.5 * (fieldArray[1:, 1:] + fieldArray[0:-1, 1:]) +
                       0.5 * (fieldArray[1:, 0:-1] + fieldArray[0:-1, 0:-1]))

            plotHandle = plt.pcolormesh(x, y, interpFieldArray,
                                        cmap=colormapDict['colormap'],
                                        norm=colormapDict['norm'],
                                        rasterized=True)
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
        plt.contourf(x, y, zeroArray, colors='white')

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

    if title is not None:
        if plotAsContours and labelContours \
           and contourComparisonFieldArray is None:
            title = limit_title(title, maxTitleLength-(3+len(colorbarLabel)))
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

    if (xlabel is not None) or (ylabel is not None):
        if axisFontSize is None:
            axisFontSize = config.get('plot', 'axisFontSize')
        axis_font = {'size': axisFontSize}

    if xlabel is not None:
        plt.xlabel(xlabel, **axis_font)
    if ylabel is not None:
        plt.ylabel(ylabel, **axis_font)

    if invertYAxis:
        ax.invert_yaxis()

    if xLim:
        ax.set_xlim(xLim)
    if yLim:
        ax.set_ylim(yLim)

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
        xlimits = ax.get_xlim()
        ax2.set_xlim(xlimits)
        xticks = np.linspace(xlimits[0], xlimits[1], numUpperTicks)
        tickValues = np.interp(xticks, x.flatten()[:num_x], secondXAxisData)
        ax2.set_xticks(xticks)
        formatString = "{{0:.{:d}f}}{}".format(
            upperXAxisTickLabelPrecision, r'$\degree$')
        ax2.set_xticklabels([formatString.format(member)
                             for member in tickValues])

        # add a third x-axis scale, if it was requested
        if thirdXAxisData is not None:
            ax3 = ax.twiny()
            ax3.set_facecolor(backgroundColor)
            ax3.set_xlabel(thirdXAxisLabel, **axis_font)
            ax3.set_xlim(xlimits)
            ax3.set_xticks(xticks)
            tickValues = np.interp(xticks, x.flatten()[:num_x], thirdXAxisData)
            ax3.set_xticklabels([formatString.format(member)
                                 for member in tickValues])
            ax3.spines['top'].set_position(('outward', 36))

    return fig, ax  # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
