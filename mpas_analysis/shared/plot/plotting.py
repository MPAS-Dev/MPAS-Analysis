import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as cols
import matplotlib.dates as mdates
import numpy as np

"""
Plotting utilities, including routine for plotting:
    * time series (and comparing with reference data sets)
    * regridded horizontal fields (and comparing with reference data sets)

Author: Xylar Asay-Davis

Last Modified: 02/11/2017
"""


def timeseries_analysis_plot(config, dsvalues, N, title, xlabel, ylabel,
                             fileout, lineStyles, lineWidths,
                             titleFontSize=None, figsize=(15, 6), dpi=300):

    """
    Plots the list of time series data sets and stores the result in an image
    file.

    config is an instance of ConfigParser

    dsvalues is a list of xarray DataSets to be plotted

    N is the numer of time points over which to perform a moving average

    title is a string with the title of the plot

    xlabel and ylabel are strings with axes labels

    fileout is the file name to be written

    lineStyles and lineWidths are lists of strings controling line style/width

    titleFontSize is the size of the title font

    figsize is the size of the figure in inches

    dpi is the number of dots per inch of the figure

    Author: Xylar Asay-Davis

    Last Modified: 02/11/2017
    """
    fig = plt.figure(figsize=figsize, dpi=dpi)

    for dsIndex in range(len(dsvalues)):
        dsvalue = dsvalues[dsIndex]
        if dsvalue is None:
            continue
        mean = pd.Series.rolling(dsvalue.to_pandas(), N, center=True).mean()
        mean = xr.DataArray.from_series(mean)
        plt.plot_date(mean['Time'], mean,
                      fmt=lineStyles[dsIndex],
                      xdate=True,
                      ydate=False,
                      linewidth=lineWidths[dsIndex])

    ax = plt.gca()
    # TODO: it would be good to change labels to include months if the plot is
    # shorter than a couple of years
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    fig.autofmt_xdate()

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
    title = None,
    plotProjection = "npstere",
    latmin =  50.0,
    lon0 = 0,
    modelTitle = "Model",
    obsTitle = "Observations",
    diffTitle = "Model-Observations",
    cbarlabel = "units",
    titleFontSize = None,
    figsize = (8,22),
    dpi = 300):

    """
    Plots a data set around either the north or south pole.

    config is an instance of ConfigParser

    Lons and Lats are the longitude and latitude arrays

    modelArray and obsArray are the model and observational data sets

    diffArray is the difference between modelArray and obsArray

    cmapModelObs is the colormap of model and observations

    clevsModleObs are contour values for model and observations

    cmapDiff is the colormap of difference

    clevsDiff are contour values fordifference

    fileout is the file name to be written

    title is a string with the subtitle of the plot

    plotProjection is a Basemap projection for the plot

    modelTitle is the title of the model plot

    obsTitle is the title of the observations plot

    diffTitle is the title of the difference plot

    cbarlabel is the label on the colorbar

    titleFontSize is the size of the title font

    figsize is the size of the figure in inches

    dpi is the number of dots per inch of the figure

    Author: Xylar Asay-Davis

    Last Modified: 02/02/2017
    """

    # set up figure
    fig = plt.figure(figsize=figsize, dpi=dpi)
    if (title is not None):
        if titleFontSize is None:
            titleFontSize = config.get('plot', 'titleFontSize')
        title_font = {'size': titleFontSize,
                      'color':config.get('plot', 'titleFontColor'),
                      'weight':config.get('plot', 'titleFontWeight')}
        fig.suptitle(title, y=0.95, **title_font)
    axis_font = {'size':config.get('plot', 'axisFontSize')}

    m = Basemap(projection=plotProjection,boundinglat=latmin,lon_0=lon0,resolution='l')
    x, y = m(Lons, Lats) # compute map proj coordinates

    normModelObs = cols.BoundaryNorm(clevsModelObs, cmapModelObs.N)
    normDiff = cols.BoundaryNorm(clevsDiff, cmapDiff.N)

    plt.subplot(3,1,1)
    plt.title(modelTitle, y=1.06, **axis_font)
    m.drawcoastlines()
    m.fillcontinents(color='grey',lake_color='white')
    m.drawparallels(np.arange(-80.,81.,10.))
    m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,True,True,True])
    cs = m.contourf(x,y,modelArray,cmap=cmapModelObs,norm=normModelObs,spacing='uniform',levels=clevsModelObs)
    cbar = m.colorbar(cs,location='right',pad="15%",spacing='uniform',ticks=clevsModelObs,boundaries=clevsModelObs)
    #cbar = m.colorbar(cs,location='right',pad="15%",spacing='uniform',extendfrac='auto',
    #                  extendrect='True',ticks=clevsModelObs, boundaries=clevsModelObs)
    cbar.set_label(cbarlabel)

    plt.subplot(3,1,2)
    plt.title(obsTitle, y=1.06, **axis_font)
    m.drawcoastlines()
    m.fillcontinents(color='grey',lake_color='white')
    m.drawparallels(np.arange(-80.,81.,10.))
    m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,True,True,True])
    cs = m.contourf(x,y,obsArray,cmap=cmapModelObs,norm=normModelObs,spacing='uniform',levels=clevsModelObs)
    cbar = m.colorbar(cs,location='right',pad="15%",spacing='uniform',ticks=clevsModelObs,boundaries=clevsModelObs)
    #cbar = m.colorbar(cs,location='right',pad="15%",spacing='uniform',extendfrac='auto',
    #                  extendrect='True',ticks=clevsModelObs, boundaries=clevsModelObs)
    cbar.set_label(cbarlabel)

    plt.subplot(3,1,3)
    plt.title(diffTitle, y=1.06, **axis_font)
    m.drawcoastlines()
    m.fillcontinents(color='grey',lake_color='white')
    m.drawparallels(np.arange(-80.,81.,10.))
    m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,True,True,True])
    cs = m.contourf(x,y,diffArray,cmap=cmapDiff,norm=normDiff,spacing='uniform',levels=clevsDiff)
    cbar = m.colorbar(cs,location='right',pad="15%",spacing='uniform',ticks=clevsDiff,boundaries=clevsModelObs)
    #cbar = m.colorbar(cs,location='right',pad="15%",spacing='uniform',extendfrac='auto',
    #                  extendrect='True',ticks=clevsDiff, boundaries=clevsDiff)
    cbar.set_label(cbarlabel)

    if (fileout is not None):
        plt.savefig(fileout,dpi=dpi,bbox_inches='tight',pad_inches=0.1)

    if not config.getboolean('plot','displayToScreen'):
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
    title = None,
    modelTitle = "Model",
    obsTitle = "Observations",
    diffTitle = "Model-Observations",
    cbarlabel = "units",
    titleFontSize = None,
    figsize = (8,12),
    dpi = 300):

    """
    Plots a data set as a longitude/latitude map.

    config is an instance of ConfigParser

    Lons and Lats are the longitude and latitude arrays

    modelArray and obsArray are the model and observational data sets

    diffArray is the difference between modelArray and obsArray

    cmapModelObs is the colormap of model and observations

    clevsModleObs are contour values for model and observations

    cmapDiff is the colormap of difference

    clevsDiff are contour values fordifference

    fileout is the file name to be written

    title is a string with the subtitle of the plot

    modelTitle is the title of the model plot

    obsTitle is the title of the observations plot

    diffTitle is the title of the difference plot

    cbarlabel is the label on the colorbar

    titleFontSize is the size of the title font

    figsize is the size of the figure in inches

    dpi is the number of dots per inch of the figure

    Author: Xylar Asay-Davis

    Last Modified: 02/02/2017
    """
    # set up figure
    fig = plt.figure(figsize=figsize, dpi=dpi)
    if (title is not None):
        if titleFontSize is None:
            titleFontSize = config.get('plot', 'titleFontSize')
        title_font = {'size': titleFontSize,
                      'color':config.get('plot', 'titleFontColor'),
                      'weight':config.get('plot', 'titleFontWeight')}
        fig.suptitle(title, y=0.95, **title_font)
    axis_font = {'size':config.get('plot', 'axisFontSize')}

    m = Basemap(projection='cyl',llcrnrlat=-85,urcrnrlat=86,llcrnrlon=-180,urcrnrlon=181,resolution='l')
    #m = Basemap(projection='robin',lon_0=200,resolution='l') # this doesn't work because lons are -180 to 180..
    x, y = m(Lons, Lats) # compute map proj coordinates

    normModelObs = cols.BoundaryNorm(clevsModelObs, cmapModelObs.N)
    normDiff = cols.BoundaryNorm(clevsDiff, cmapDiff.N)

    plt.subplot(3,1,1)
    plt.title(modelTitle, y=1.06, **axis_font)
    m.drawcoastlines()
    m.fillcontinents(color='grey',lake_color='white')
    m.drawparallels(np.arange(-80.,80.,20.),labels=[True,False,False,False])
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[False,False,False,True])
    cs = m.contourf(x,y,modelArray,cmap=cmapModelObs,norm=normModelObs,spacing='uniform',levels=clevsModelObs,extend='both')
    cbar = m.colorbar(cs,location='right',pad="5%",spacing='uniform',ticks=clevsModelObs,boundaries=clevsModelObs)
    #cbar = m.colorbar(cs,location='right',pad="5%",spacing='uniform',extendfrac='auto',
    #                  extendrect='True',ticks=clevsModelObs, boundaries=clevsModelObs)
    cbar.set_label(cbarlabel)

    plt.subplot(3,1,2)
    plt.title(obsTitle, y=1.06, **axis_font)
    m.drawcoastlines()
    m.fillcontinents(color='grey',lake_color='white')
    m.drawparallels(np.arange(-80.,80.,20.),labels=[True,False,False,False])
    m.drawmeridians(np.arange(-180.,180.,40.),labels=[False,False,False,True])
    cs = m.contourf(x,y,obsArray,cmap=cmapModelObs,norm=normModelObs,spacing='uniform',levels=clevsModelObs,extend='both')
    cbar = m.colorbar(cs,location='right',pad="5%",spacing='uniform',ticks=clevsModelObs,boundaries=clevsModelObs)
    #cbar = m.colorbar(cs,location='right',pad="5%",spacing='uniform',extendfrac='auto',
    #                  extendrect='True',ticks=clevsModelObs, boundaries=clevsModelObs)
    cbar.set_label(cbarlabel)

    plt.subplot(3,1,3)
    plt.title(diffTitle, y=1.06, **axis_font)
    m.drawcoastlines()
    m.fillcontinents(color='grey',lake_color='white')
    m.drawparallels(np.arange(-80.,80.,20.),labels=[True,False,False,False])
    m.drawmeridians(np.arange(-180.,180.,40.),labels=[False,False,False,True])
    cs = m.contourf(x,y,diffArray,cmap=cmapDiff,norm=normDiff,spacing='uniform',levels=clevsDiff,extend='both')
    cbar = m.colorbar(cs,location='right',pad="5%",spacing='uniform',ticks=clevsDiff,boundaries=clevsModelObs)
    #cbar = m.colorbar(cs,location='right',pad="5%",spacing='uniform',extendfrac='auto',
    #                  extendrect='True',ticks=clevsDiff, boundaries=clevsDiff)
    cs.cmap.set_over((1., 1., 1.))
    cs.cmap.set_under((0., 0., 0.))
    cbar.set_label(cbarlabel)

    if (fileout is not None):
        plt.savefig(fileout,dpi=dpi,bbox_inches='tight',pad_inches=0.1)

    if not config.getboolean('plot','displayToScreen'):
        plt.close()
