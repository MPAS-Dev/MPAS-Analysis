#!/usr/bin/env python
"""
General comparison of 2-d model fields against data.  Currently only supports
mixed layer depths (mld) and sea surface temperature (sst)

Author: Luke Van Roekel, Milena Veneziani, Xylar Asay-Davis
Last Modified: 12/06/2016
"""

import matplotlib.pyplot as plt
import matplotlib.colors as cols

import numpy as np
import xarray as xr
import datetime
from netCDF4 import Dataset as netcdf_dataset

from ..shared.mpas_xarray.mpas_xarray import preprocess_mpas, \
    remove_repeated_time_index
from ..shared.plot.plotting import plot_global_comparison
from ..shared.interpolation.interpolate import interp_fields, init_tree
from ..shared.constants import constants

from ..shared.io import StreamsFile


def ocn_modelvsobs(config, field, streamMap=None, variableMap=None):

    """
    Plots a comparison of ACME/MPAS output to SST or MLD observations

    config is an instance of MpasAnalysisConfigParser containing configuration
    options.

    field is the name of a field to be analyize (currently one of 'mld' or
    'sst')

    If present, streamMap is a dictionary of MPAS-O stream names that map to
    their mpas_analysis counterparts.

    If present, variableMap is a dictionary of MPAS-O variable names that map
    to their mpas_analysis counterparts.

    Authors: Luke Van Roekel, Milena Veneziani, Xylar Asay-Davis
    Modified: 12/08/2016
    """

    # read parameters from config file
    indir = config.get('paths', 'archive_dir_ocn')

    streams_filename = config.get('input', 'ocean_streams_filename')
    streams = StreamsFile(streams_filename, streamsdir=indir)

    # get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('time', 'climo_start_date')
    endDate = config.get('time', 'climo_end_date')
    streamName = streams.find_stream(streamMap['timeSeriesStats'])
    infiles = streams.readpath(streamName, startDate=startDate,
                               endDate=endDate)
    print 'Reading files {} through {}'.format(infiles[0], infiles[-1])

    plots_dir = config.get('paths', 'plots_dir')
    obsdir = config.get('paths', 'obs_' + field + 'dir')
    casename = config.get('case', 'casename')
    meshfile = config.get('data', 'mpas_meshfile')
    climo_yr1 = config.getint('time', 'climo_yr1')
    climo_yr2 = config.getint('time', 'climo_yr2')
    yr_offset = config.getint('time', 'yr_offset')

    outputTimes = config.getExpression(field + '_modelvsobs',
                                       'comparisonTimes')

    f = netcdf_dataset(meshfile, mode='r')
    lonCell = f.variables["lonCell"][:]
    latCell = f.variables["latCell"][:]

    varList = [field]

    if field == 'mld':

        selvals = None

        # Load MLD observational data
        obs_filename = "{}/holtetalley_mld_climatology.nc".format(obsdir)
        dsData = xr.open_mfdataset(obs_filename)

        # Increment month value to be consistent with the model output
        dsData.iMONTH.values += 1

        # Rename the time dimension to be consistent with the SST dataset
        dsData.rename({'month': 'calmonth'}, inplace=True)
        dsData.rename({'iMONTH': 'month'}, inplace=True)

        obsFieldName = 'mld_dt_mean'

        # Reorder dataset for consistence
        dsData = dsData.transpose('month', 'iLON', 'iLAT')

        # Set appropriate MLD figure labels
        obsTitleLabel = "Observations (HolteTalley density threshold MLD)"
        fileOutLabel = "mldHolteTalleyARGO"
        unitsLabel = 'm'

    elif field == 'sst':

        selvals = {'nVertLevels': 0}

        obs_filename = \
            "{}/MODEL.SST.HAD187001-198110.OI198111-201203.nc".format(obsdir)
        dsData = xr.open_mfdataset(obs_filename)
        # Select years for averaging (pre-industrial or present-day)
        # This seems fragile as definitions can change
        if yr_offset < 1900:
            time_start = datetime.datetime(1870, 1, 1)
            time_end = datetime.datetime(1900, 12, 31)
            preIndustrial_txt = "pre-industrial 1870-1900"
        else:
            time_start = datetime.datetime(1990, 1, 1)
            time_end = datetime.datetime(2011, 12, 31)
            preIndustrial_txt = "present-day 1990-2011"

        ds_tslice = dsData.sel(time=slice(time_start, time_end))
        monthly_clim_data = ds_tslice.groupby('time.month').mean('time')

        # Rename the observation data for code compactness
        dsData = monthly_clim_data.transpose('month', 'lon', 'lat')
        obsFieldName = 'SST'

        # Set appropriate figure labels for SST
        obsTitleLabel = \
            "Observations (Hadley/OI, {})".format(preIndustrial_txt)
        fileOutLabel = "sstHADOI"
        unitsLabel = r'$^o$C'

    ds = xr.open_mfdataset(
        infiles,
        preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset,
                                             timestr='Time',
                                             onlyvars=varList,
                                             selvals=selvals,
                                             variable_map=variableMap))
    ds = remove_repeated_time_index(ds)

    time_start = datetime.datetime(yr_offset+climo_yr1, 1, 1)
    time_end = datetime.datetime(yr_offset+climo_yr2, 12, 31)
    ds_tslice = ds.sel(Time=slice(time_start, time_end))
    monthly_clim = ds_tslice.groupby('Time.month').mean('Time')

    latData, lonData = np.meshgrid(dsData.lat.values, dsData.lon.values)
    latData = latData.flatten()
    lonData = lonData.flatten()

    daysarray = np.ones((12, dsData[obsFieldName].values.shape[1],
                         dsData[obsFieldName].values.shape[2]))

    for i, dval in enumerate(constants.dinmonth):
        daysarray[i, :, :] = dval
        inds = np.where(np.isnan(dsData[obsFieldName][i, :, :].values))
        daysarray[i, inds[0], inds[1]] = np.NaN

    # initialize interpolation variables
    d2, inds2, lonTarg, latTarg = init_tree(np.rad2deg(lonCell),
                                            np.rad2deg(latCell),
                                            constants.lonmin,
                                            constants.lonmax,
                                            constants.latmin,
                                            constants.latmax,
                                            constants.dLongitude,
                                            constants.dLatitude)
    d, inds, lonTargD, latTargD = init_tree(lonData, latData,
                                            constants.lonmin,
                                            constants.lonmax,
                                            constants.latmin,
                                            constants.latmax,
                                            constants.dLongitude,
                                            constants.dLatitude)
    nLon = lonTarg.shape[0]
    nLat = lonTarg.shape[1]

    modelOutput = np.zeros((len(outputTimes), nLon, nLat))
    observations = np.zeros((len(outputTimes), nLon, nLat))
    bias = np.zeros((len(outputTimes), nLon, nLat))

    # Interpolate and compute biases
    for i, timestring in enumerate(outputTimes):
        monthsvalue = constants.monthdictionary[timestring]

        if isinstance(monthsvalue, (int, long)):
            modeldata = monthly_clim.sel(month=monthsvalue)[field].values
            obsdata = dsData.sel(month=monthsvalue)[obsFieldName].values
        else:

            modeldata = (np.sum(
                constants.dinmonth[monthsvalue-1] *
                monthly_clim.sel(month=monthsvalue)[field].values.T, axis=1) /
                np.sum(constants.dinmonth[monthsvalue-1]))
            obsdata = (np.nansum(
                daysarray[monthsvalue-1, :, :] *
                dsData.sel(month=monthsvalue)[obsFieldName].values, axis=0) /
                np.nansum(daysarray[monthsvalue-1, :, :], axis=0))

        modelOutput[i, :, :] = interp_fields(modeldata, d2, inds2, lonTarg)
        observations[i, :, :] = interp_fields(obsdata.flatten(), d, inds,
                                              lonTargD)

    for i in range(len(outputTimes)):
        bias[i, :, :] = modelOutput[i, :, :] - observations[i, :, :]

    clevsModelObs = config.getExpression(field + '_modelvsobs',
                                         'clevsModelObs')
    cmap = plt.get_cmap(config.get(field + '_modelvsobs',
                                   'cmapModelObs'))
    cmapIndices = config.getExpression(field + '_modelvsobs',
                                       'cmapIndicesModelObs')
    cmapModelObs = cols.ListedColormap(cmap(cmapIndices), "cmapModelObs")
    clevsDiff = config.getExpression(field + '_modelvsobs',
                                     'clevsDiff')
    cmap = plt.get_cmap(config.get(field + '_modelvsobs', 'cmapDiff'))
    cmapIndices = config.getExpression(field + '_modelvsobs',
                                       'cmapIndicesDiff')
    cmapDiff = cols.ListedColormap(cmap(cmapIndices), "cmapDiff")

    for i in range(len(outputTimes)):
        fileout = "{}/{}_{}_{}_years{:04d}-{:04d}.png".format(
            plots_dir, fileOutLabel, casename, outputTimes[i], climo_yr1,
            climo_yr2)
        title = "{} ({}, years {:04d}-{:04d})".format(
            field.upper(), outputTimes[i], climo_yr1, climo_yr2)
        plot_global_comparison(config,
                               lonTarg,
                               latTarg,
                               modelOutput[i, :, :],
                               observations[i, :, :],
                               bias[i, :, :],
                               cmapModelObs,
                               clevsModelObs,
                               cmapDiff,
                               clevsDiff,
                               fileout=fileout,
                               title=title,
                               modelTitle="{}".format(casename),
                               obsTitle=obsTitleLabel,
                               diffTitle="Model-Observations",
                               cbarlabel=unitsLabel)
