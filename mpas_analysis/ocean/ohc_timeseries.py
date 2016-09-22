import numpy as np
from netCDF4 import Dataset as netcdf_dataset
import xarray as xr
import pandas as pd
import datetime

from ..shared.mpas_xarray.mpas_xarray import preprocess_mpas, remove_repeated_time_index

from ..shared.plot.plotting import timeseries_analysis_plot


def ohc_timeseries(config):
    # read parameters from config file
    indir = config.get('paths','archive_dir_ocn')
    meshfile = config.get('data','mpas_meshfile')
    casename = config.get('case','casename')
    ref_casename_v0 = config.get('case','ref_casename_v0')
    indir_v0data = config.get('paths','ref_archive_v0_ocndir')

    compare_with_obs = config.getboolean('ohc_timeseries','compare_with_obs')

    plots_dir = config.get('paths','plots_dir')

    yr_offset = config.getint('time','yr_offset')

    N_movavg = config.getint('ohc_timeseries','N_movavg')

    regions = config.getlist('regions','regions',listType=str)
    plot_titles = config.getlist('regions','plot_titles',listType=str)
    iregions = config.getlist('ohc_timeseries','regionIndicesToPlot',listType=int)

    # Define/read in general variables
    print "  Read in depth and compute specific depth indexes..."
    f = netcdf_dataset(meshfile,mode='r')
    depth = f.variables["refBottomDepth"][:] # reference depth [m]
    cp = f.getncattr("config_specific_heat_sea_water") # specific heat [J/(kg*degC)]
    rho = f.getncattr("config_density0") # [kg/m3]
    fac = 1e-22*rho*cp;

    k700m = np.where(depth>700.)[0][0]-1
    k2000m = np.where(depth>2000.)[0][0]-1

    kbtm = len(depth)-1

    # Load data
    print "  Load ocean data..."
    infiles = '%s/am.mpas-o.timeSeriesStats.????-??*nc'%indir

    # Load data:
    ds = xr.open_mfdataset(infiles, preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset,
                         timeSeriesStats=True,
                         timestr='time_avg_daysSinceStartOfSim',
                         onlyvars=['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
                                   'time_avg_avgValueWithinOceanLayerRegion_sumLayerMaskValue',
                                   'time_avg_avgValueWithinOceanLayerRegion_avgLayerArea',
                                   'time_avg_avgValueWithinOceanLayerRegion_avgLayerThickness']))


    ds = remove_repeated_time_index(ds)

    # Select year-1 data and average it (for later computing anomalies)
    time_start = datetime.datetime(yr_offset+1,1,1)
    time_end = datetime.datetime(yr_offset+1,12,31)
    ds_yr1 = ds.sel(Time=slice(time_start,time_end))
    mean_yr1 = ds_yr1.mean('Time')


    print "  Compute temperature anomalies..."
    avgLayerTemperature = ds.time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature
    avgLayerTemperature_yr1 = mean_yr1.time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature

    avgLayTemp_anomaly = avgLayerTemperature - avgLayerTemperature_yr1

    year_start = (pd.to_datetime(ds.Time.min().values)).year
    year_end   = (pd.to_datetime(ds.Time.max().values)).year
    time_start = datetime.datetime(year_start,1,1)
    time_end   = datetime.datetime(year_end,12,31)

    if ref_casename_v0 != "None":
        print "  Load in OHC for ACMEv0 case..."
        infiles_v0data = "".join([indir_v0data,'/OHC.',ref_casename_v0,'.year*.nc'])
        ds_v0 = xr.open_mfdataset(infiles_v0data,preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset))
        ds_v0 = remove_repeated_time_index(ds_v0)
        ds_v0_tslice = ds_v0.sel(Time=slice(time_start,time_end))

    sumLayerMaskValue = ds.time_avg_avgValueWithinOceanLayerRegion_sumLayerMaskValue
    avgLayerArea = ds.time_avg_avgValueWithinOceanLayerRegion_avgLayerArea
    avgLayerThickness = ds.time_avg_avgValueWithinOceanLayerRegion_avgLayerThickness

    print "  Compute OHC and make plots..."
    for index in range(len(iregions)):
        iregion = iregions[index]

        # Compute volume of each layer in the region:
        layerArea = sumLayerMaskValue[:,iregion,:] * avgLayerArea[:,iregion,:]
        layerVolume = layerArea * avgLayerThickness[:,iregion,:]

        # Compute OHC:
        ohc = layerVolume * avgLayTemp_anomaly[:,iregion,:]
        # OHC over 0-bottom depth range:
        ohc_tot = ohc.sum('nVertLevels')
        ohc_tot = fac*ohc_tot

        # OHC over 0-700m depth range:
        ohc_700m = fac*ohc[:,0:k700m].sum('nVertLevels')

        # OHC over 700m-2000m depth range:
        ohc_2000m = fac*ohc[:,k700m+1:k2000m].sum('nVertLevels')

        # OHC over 2000m-bottom depth range:
        ohc_btm = ohc[:,k2000m+1:kbtm].sum('nVertLevels')
        ohc_btm = fac*ohc_btm

        title = 'OHC, %s, 0-bottom (thick-), 0-700m (thin-), 700-2000m (--), 2000m-bottom (-.) \n %s'%(plot_titles[iregion], casename)

        xlabel = "Time [years]"
        ylabel = "[x$10^{22}$ J]"

        if ref_casename_v0 != "None":
            figname = "%s/ohc_%s_%s_%s.png" % (plots_dir,regions[iregion],casename,ref_casename_v0)
            ohc_v0_tot = ds_v0_tslice.ohc_tot
            ohc_v0_700m = ds_v0_tslice.ohc_700m
            ohc_v0_2000m = ds_v0_tslice.ohc_2000m
            ohc_v0_btm = ds_v0_tslice.ohc_btm
            title = "".join([title," (r), ",ref_casename_v0," (b)"])
            timeseries_analysis_plot(config, [ohc_tot, ohc_700m, ohc_2000m, ohc_btm,
                                              ohc_v0_tot, ohc_v0_700m, ohc_v0_2000m, ohc_v0_btm],
                                     N_movavg, title, xlabel, ylabel, figname,
                                     lineStyles = ['r-', 'r-', 'r--', 'r-.', 'b-', 'b-', 'b--', 'b-.'],
                                     lineWidths = [2, 1, 1.5, 1.5, 2, 1, 1.5, 1.5])

        if not compare_with_obs and ref_casename_v0 == "None":
            figname = "%s/ohc_%s_%s.png" % (plots_dir,regions[iregion],casename)
            timeseries_analysis_plot(config, [ohc_tot, ohc_700m, ohc_2000m, ohc_btm],
                                     N_movavg, title, xlabel, ylabel, figname,
                                     lineStyles = ['r-', 'r-', 'r--', 'r-.'],
                                     lineWidths = [2, 1, 1.5, 1.5])
