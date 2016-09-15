import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as cols
import numpy as np

def timeseries_analysis_plot(config, dsvalues, N, title, xlabel, ylabel, fileout,
                             lineStyles, lineWidths, title_font_size=None,
                             figsize=(15,6), dpi=300):

    plt.figure(figsize=figsize, dpi=dpi)

    for dsIndex in range(len(dsvalues)):
        dsvalue = dsvalues[dsIndex]
        if dsvalue is None:
            continue
        mean = pd.Series.rolling(dsvalue.to_pandas(), N, center=True).mean()
        mean.plot(style=lineStyles[dsIndex], lw=lineWidths[dsIndex])

    if title_font_size is None:
        title_font_size = config.get('plot', 'title_font_size')

    axis_font = {'size':config.get('plot', 'axis_font_size')}
    title_font = {'size': title_font_size,
                  'color':config.get('plot', 'title_font_color'),
                  'weight':config.get('plot', 'title_font_weight')}
    if (title != None):
        plt.title(title, **title_font)
    if (xlabel != None):
        plt.xlabel(xlabel, **axis_font)
    if (ylabel != None):
        plt.ylabel(ylabel, **axis_font)
    if (fileout is not None):
        plt.savefig(fileout,dpi=dpi,bbox_inches='tight',pad_inches=0.1)

    if not config.getboolean('plot','displayToScreen'):
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
    title_font_size = None,
    figsize = (8,22),
    dpi = 300):

    # set up figure
    fig = plt.figure(figsize=figsize, dpi=dpi)
    if (title is not None):
        if title_font_size is None:
            title_font_size = config.get('plot', 'title_font_size')
        title_font = {'size': title_font_size,
                      'color':config.get('plot', 'title_font_color'),
                      'weight':config.get('plot', 'title_font_weight')}
        fig.suptitle(title, y=0.95, **title_font)
    axis_font = {'size':config.get('plot', 'axis_font_size')}

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
    title_font_size = None,
    figsize = (8,12),
    dpi = 300):

    # set up figure
    fig = plt.figure(figsize=figsize, dpi=dpi)
    if (title is not None):
        if title_font_size is None:
            title_font_size = config.get('plot', 'title_font_size')
        title_font = {'size': title_font_size,
                      'color':config.get('plot', 'title_font_color'),
                      'weight':config.get('plot', 'title_font_weight')}
        fig.suptitle(title, y=0.95, **title_font)
    axis_font = {'size':config.get('plot', 'axis_font_size')}
        
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
