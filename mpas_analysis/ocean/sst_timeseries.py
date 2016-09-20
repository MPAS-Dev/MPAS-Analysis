import xarray as xr
import pandas as pd
import datetime

from ..shared.mpas_xarray.mpas_xarray import preprocess_mpas, remove_repeated_time_index

from ..shared.plot.plotting import timeseries_analysis_plot

def sst_timeseries(config):

    # Define/read in general variables
    print "  Load SST data..."
    indir = config.get('paths','archive_dir_ocn')
    infiles = '%s/am.mpas-o.timeSeriesStats.????-??*nc'%indir

    casename = config.get('case','casename')
    ref_casename_v0 = config.get('case','ref_casename_v0')
    indir_v0data = config.get('paths','ref_archive_v0_ocndir')

    compare_with_obs = config.getboolean('sst_timeseries','compare_with_obs')

    plots_dir = config.get('paths','plots_dir')

    yr_offset = config.getint('time','yr_offset')

    N_movavg = config.getint('sst_timeseries','N_movavg')

    regions = config.getlist('regions','regions',listType=str)
    plot_titles = config.getlist('regions','plot_titles',listType=str)
    iregions = config.getlist('sst_timeseries','regionIndicesToPlot',listType=int)

    # Load data:
    ds = xr.open_mfdataset(infiles, preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset,
                           timeSeriesStats=True, timestr='time_avg_daysSinceStartOfSim',
                           onlyvars=['time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature']))
    ds = remove_repeated_time_index(ds)

    SSTregions = ds.time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature

    year_start = (pd.to_datetime(ds.Time.min().values)).year
    year_end   = (pd.to_datetime(ds.Time.max().values)).year
    time_start = datetime.datetime(year_start,1,1)
    time_end   = datetime.datetime(year_end,12,31)


    if ref_casename_v0 != "None":
        print "  Load in SST for ACMEv0 case..."
        infiles_v0data = "".join([indir_v0data,'/SST.',ref_casename_v0,'.year*.nc'])
        ds_v0 = xr.open_mfdataset(infiles_v0data,preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset))
        ds_v0 = remove_repeated_time_index(ds_v0)
        ds_v0_tslice = ds_v0.sel(Time=slice(time_start,time_end))

    print "  Make plots..."
    for index in range(len(iregions)):
        iregion = iregions[index]

        title = plot_titles[iregion]
        title = "SST, %s, %s (r-)" % (title, casename)
        xlabel = "Time [years]"
        ylabel = "[$^\circ$ C]"

        SST = SSTregions[:,iregion]

        if ref_casename_v0 != "None":
            figname = "%s/sst_%s_%s_%s.png" % (plots_dir,regions[iregion],casename,ref_casename_v0)
            SST_v0 = ds_v0_tslice.SST

            title = "%s\n %s (b-)" % (title, ref_casename_v0)
            timeseries_analysis_plot(config, [SST,SST_v0], N_movavg,
                                     title, xlabel, ylabel, figname,
                                     lineStyles = ['r-','b-'],
                                     lineWidths = [1.2,1.2])
        else:
            figname = "%s/sst_%s_%s.png" % (plots_dir,regions[iregion],casename)
            timeseries_analysis_plot(config, [SST], N_movavg, title, xlabel, ylabel, figname,
                                     lineStyles = ['r-'], lineWidths = [1.2])
