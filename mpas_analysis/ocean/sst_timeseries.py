import xarray as xr
import pandas as pd
import datetime
from netCDF4 import Dataset as netcdf_dataset

from ..shared.mpas_xarray.mpas_xarray import preprocess_mpas, \
    remove_repeated_time_index

from ..shared.plot.plotting import timeseries_analysis_plot

from ..shared.io import StreamsFile

from ..shared.timekeeping.Date import Date


def sst_timeseries(config, streamMap=None, variableMap=None):
    """
    Performs analysis of the time-series output of sea-surface temperature
    (SST).

    config is an instance of `MpasAnalysisConfigParser` containing
    configuration options.

    If present, `streamMap` is a dictionary of MPAS-O stream names that map to
    their mpas_analysis counterparts.

    If present, `variableMap` is a dictionary of MPAS-O variable names that map
    to their mpas_analysis counterparts.

    Author: Xylar Asay-Davis, Milena Veneziani

    Last Modified: 01/26/2017
    """

    # Define/read in general variables
    print '  Load SST data...'
    # read parameters from config file
    indir = config.get('paths', 'archive_dir_ocn')

    streams_filename = config.get('input', 'ocean_streams_filename')
    streams = StreamsFile(streams_filename, streamsdir=indir)

    try:
        restartFile = streams.readpath('restart')[0]
    except ValueError:
        raise IOError('No MPAS-O restart file found: need at least one '
                      'restart file for ocn_modelvsobs calculation')

    f = netcdf_dataset(restartFile, mode='r')
    inRefDate = ''.join(f.variables["simulationStartTime"][:]).strip()
    # convert the MPAs date to a CF-compliant date
    inRefDate = inRefDate.replace('_', ' ')
    simStartYear = int(inRefDate[0:4])

    yr_offset = config.getint('time', 'yr_offset')
    outRefDate = '{:04d}-01-01 00:00:00'.format(yr_offset+simStartYear)

    print inRefDate, outRefDate

    # get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('time', 'timeseries_start_date')
    endDate = config.get('time', 'timeseries_end_date')
    streamName = streams.find_stream(streamMap['timeSeriesStats'])
    infiles = streams.readpath(streamName, startDate=startDate,
                               endDate=endDate)
    print 'Reading files {} through {}'.format(infiles[0], infiles[-1])

    casename = config.get('case', 'casename')
    ref_casename_v0 = config.get('case', 'ref_casename_v0')
    indir_v0data = config.get('paths', 'ref_archive_v0_ocndir')

    plots_dir = config.get('paths', 'plots_dir')

    N_movavg = config.getint('sst_timeseries', 'N_movavg')

    regions = config.getExpression('regions', 'regions')
    plot_titles = config.getExpression('regions', 'plot_titles')
    iregions = config.getExpression('sst_timeseries', 'regionIndicesToPlot')

    # Load data:
    varList = ['avgSurfaceTemperature']
    ds = xr.open_mfdataset(
        infiles,
        preprocess=lambda x: preprocess_mpas(x, inrefdate=inRefDate,
                                             outrefdate=outRefDate,
                                             timestr='Time',
                                             onlyvars=varList,
                                             varmap=variableMap),
        decode_times=False)
    ds = remove_repeated_time_index(ds)

    # convert the start and end dates to netcdftime.datetime objects using
    # the Date class and adding the yr_offset
    offset = Date(dateString='{:04d}-00-00'.format(yr_offset), isInterval=True)
    offset = offset.to_timedelta()
    time_start = Date(startDate).to_netcdftime() + offset
    time_end = Date(endDate).to_netcdftime() + offset

    # select only the data in the specified range of years
    ds = ds.sel(Time=slice(time_start, time_end))

    SSTregions = ds.avgSurfaceTemperature

    year_start = (pd.to_datetime(ds.Time.min().values)).year
    year_end = (pd.to_datetime(ds.Time.max().values)).year
    time_start = datetime.datetime(year_start, 1, 1)
    time_end = datetime.datetime(year_end, 12, 31)

    if ref_casename_v0 != 'None':
        print '  Load in SST for ACMEv0 case...'
        inRefDate = '0001-01-01 00:00:00'

        infiles_v0data = '{}/SST.{}.year*.nc'.format(indir_v0data,
                                                     ref_casename_v0)
        ds_v0 = xr.open_mfdataset(
            infiles_v0data,
            preprocess=lambda x: preprocess_mpas(x, inrefdate=inRefDate,
                                                 outrefdate=outRefDate),
            decode_times=False)
        ds_v0 = remove_repeated_time_index(ds_v0)
        year_end_v0 = (pd.to_datetime(ds_v0.Time.max().values)).year
        if year_start <= year_end_v0:
            ds_v0_tslice = ds_v0.sel(Time=slice(time_start, time_end))
        else:
            print '   Warning: v0 time series lies outside current bounds ' \
                'of v1 time series. Skipping it.'
            ref_casename_v0 = 'None'

    print '  Make plots...'
    for index in range(len(iregions)):
        iregion = iregions[index]

        title = plot_titles[iregion]
        title = 'SST, %s, %s (r-)' % (title, casename)
        xlabel = 'Time [years]'
        ylabel = '[$^\circ$ C]'

        SST = SSTregions[:, iregion]

        if ref_casename_v0 != 'None':
            figname = '{}/sst_{}_{}_{}.png'.format(plots_dir, regions[iregion],
                                                   casename, ref_casename_v0)
            SST_v0 = ds_v0_tslice.SST

            title = '{}\n {} (b-)'.format(title, ref_casename_v0)
            timeseries_analysis_plot(config, [SST, SST_v0], N_movavg,
                                     title, xlabel, ylabel, figname,
                                     lineStyles=['r-', 'b-'],
                                     lineWidths=[1.2, 1.2])
        else:
            figname = '{}/sst_{}_{}.png'.format(plots_dir, regions[iregion],
                                                casename)
            timeseries_analysis_plot(config, [SST], N_movavg, title, xlabel,
                                     ylabel, figname, lineStyles=['r-'],
                                     lineWidths=[1.2])
