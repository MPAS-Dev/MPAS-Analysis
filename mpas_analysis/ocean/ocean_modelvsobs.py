#!/usr/bin/env python
"""
General comparison of 2-d model fields against data.  Currently only supports
mixed layer depths (mld) and sea surface temperature (sst)

Author: Luke Van Roekel, Milena Veneziani
Last Modified: 10/24/2016
"""

import matplotlib.pyplot as plt
import matplotlib.colors as cols

import numpy as np
import xarray as xr
import datetime
from netCDF4 import Dataset as netcdf_dataset

from ..shared.mpas_xarray.mpas_xarray import preprocess_mpas, remove_repeated_time_index
from ..shared.plot.plotting import plot_global_comparison
from ..shared.interpolation.interpolate import interp_fields, init_tree
from ..shared.constants import constants

from ..shared.io import StreamsFile
from ..shared.io.utility import paths

def ocn_modelvsobs(config, field):

    """
    Plots a comparison of ACME/MPAS output to SST or MLD observations

    Authors: Luke Van Roekel, Milena Veneziani, Xylar Asay-Davis
    Modified: 10/27/2016
    """

    # read parameters from config file
    indir = config.get('paths', 'archive_dir_ocn')

    streams_filename = config.get('input', 'ocean_streams_filename')
    streams = StreamsFile('{}/{}'.format(indir, streams_filename))

    # read the file template for timeSeriesStatsOutput, convert it to fnmatch
    # expression and make it an absolute path
    infiles = streams.readpath('timeSeriesStatsOutput', 'filename_template')
    # find files matching the fnmatch experession
    infiles = paths(infiles)

    plots_dir = config.get('paths', 'plots_dir')
    obsdir = config.get('paths', 'obs_' + field.lower() + 'dir')
    casename = config.get('case', 'casename')
    meshfile = config.get('data', 'mpas_meshfile')
    climo_yr1 = config.getint('time', 'climo_yr1')
    climo_yr2 = config.getint('time', 'climo_yr2')
    yr_offset = config.getint('time', 'yr_offset')

    #Seems like the following line should be a config.get option and not
    #read every time series file when only a subset is taken
    #infiles = "".join([indir,"/am.mpas-o.timeSeriesStats.????-*.nc"])

    outputTimes = config.getlist(field.lower() + '_modelvsobs', 'comparisonTimes')

    f = netcdf_dataset(meshfile,mode='r')
    lonCell = f.variables["lonCell"][:]
    latCell = f.variables["latCell"][:]

    if field.lower() == 'mld':
        obs_filename = "%s/holtetalley_mld_climatology.nc" % obsdir

        ds = xr.open_mfdataset(infiles, preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset,
					   timeSeriesStats=True,
                               timestr='time_avg_daysSinceStartOfSim',
                               onlyvars=['time_avg_dThreshMLD']))
        ds = remove_repeated_time_index(ds)
        ds.rename({'time_avg_dThreshMLD':'mpasData'}, inplace = True)

        #Load MLD observational data
        dsData = xr.open_mfdataset(obs_filename)

        #Increment month value to be consistent with the model output
        dsData.iMONTH.values += 1

        #Rename the time dimension to be consistent with the SST dataset
        dsData.rename({'month':'calmonth'}, inplace = True)
        dsData.rename({'iMONTH':'month'}, inplace = True)

        #rename appropriate observational data for compactness
        dsData.rename({'mld_dt_mean':'observationData'}, inplace = True)

        #Reorder dataset for consistence
        dsData = dsData.transpose('month', 'iLON', 'iLAT')

        #Set appropriate MLD figure labels
        obsTitleLabel = "Observations (HolteTalley density threshold MLD)"
        fileOutLabel = "mldHolteTalleyARGO"
        unitsLabel = 'm'

    elif field.lower() == 'sst':

        ds = xr.open_mfdataset(infiles, preprocess=lambda x: preprocess_mpas(x, yearoffset=yr_offset,
		                    timeSeriesStats=True,
		                    timestr='time_avg_daysSinceStartOfSim',
		                    onlyvars=['time_avg_activeTracers_temperature'],
		                    selvals={'nVertLevels':1}))
        ds = remove_repeated_time_index(ds)
        ds.rename({'time_avg_activeTracers_temperature':'mpasData'}, inplace = True)

        obs_filename = "%s/MODEL.SST.HAD187001-198110.OI198111-201203.nc" % obsdir
        dsData = xr.open_mfdataset(obs_filename)
        #Select years for averaging (pre-industrial or present-day)
        #This seems fragile as definitions can change
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

        #Rename the observation data for code compactness
        dsData = monthly_clim_data.transpose('month', 'lon', 'lat')
        dsData.rename({'SST':'observationData'}, inplace = True)

        #Set appropriate figure labels for SST
        obsTitleLabel = "Observations (Hadley/OI, %s)" % preIndustrial_txt
        fileOutLabel = "sstHADOI"
        unitsLabel = r'$^o$C'

    time_start = datetime.datetime(yr_offset+climo_yr1, 1, 1)
    time_end = datetime.datetime(yr_offset+climo_yr2, 12, 31)
    ds_tslice = ds.sel(Time=slice(time_start, time_end))
    monthly_clim = ds_tslice.groupby('Time.month').mean('Time')

    latData, lonData = np.meshgrid(dsData.lat.values, dsData.lon.values)
    latData = latData.flatten()
    lonData = lonData.flatten()

    daysarray = np.ones((12, dsData.observationData.values.shape[1], dsData.observationData.values.shape[2]))

    for i, dval in enumerate(constants.dinmonth):
        daysarray[i, :, :] = dval
        inds = np.where(np.isnan(dsData.observationData[i, :, :].values))
        daysarray[i, inds[0], inds[1]] = np.NaN


    # initialize interpolation variables
    d2, inds2, lonTarg, latTarg = init_tree(np.rad2deg(lonCell), np.rad2deg(latCell), constants.lonmin,
                                            constants.lonmax, constants.latmin, constants.latmax,
                                            constants.dLongitude, constants.dLatitude)
    d, inds, lonTargD, latTargD = init_tree(lonData, latData, constants.lonmin, constants.lonmax,
								constants.latmin, constants.latmax, constants.dLongitude,
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
            modeldata = monthly_clim.sel(month=monthsvalue).mpasData.values
            obsdata = dsData.sel(month = monthsvalue).observationData.values
        else:
            modeldata = np.sum(constants.dinmonth[monthsvalue-1]*monthly_clim.sel(month=monthsvalue).
        		               mpasData.values.T, axis=1) / np.sum(constants.dinmonth[monthsvalue-1])
            obsdata = np.nansum(daysarray[monthsvalue-1, :, :]*dsData.sel(month = monthsvalue).
                                observationData.values, axis=0) / np.nansum(daysarray[monthsvalue-1, :, :],
                                axis=0)

        modelOutput[i, :, :] = interp_fields(modeldata, d2, inds2, lonTarg)
        observations[i, :, :] = interp_fields(obsdata.flatten(), d, inds, lonTargD)

    for i in range(len(outputTimes)):
        bias[i, :, :] = modelOutput[i, :, :] - observations[i, :, :]

    clevsModelObs = config.getlist(field.lower() + '_modelvsobs', 'clevsModelObs', listType=float)
    cmap = plt.get_cmap(config.get(field.lower() + '_modelvsobs', 'cmapModelObs'))
    cmapIndices = config.getlist(field.lower() + '_modelvsobs', 'cmapIndicesModelObs', listType=int)
    cmapModelObs = cols.ListedColormap(cmap(cmapIndices), "cmapModelObs")
    clevsDiff = config.getlist(field.lower() + '_modelvsobs', 'clevsDiff', listType=float)
    cmap = plt.get_cmap(config.get(field.lower() + '_modelvsobs', 'cmapDiff'))
    cmapIndices = config.getlist(field.lower() + '_modelvsobs', 'cmapIndicesDiff', listType=int)
    cmapDiff = cols.ListedColormap(cmap(cmapIndices), "cmapDiff")

    for i in range(len(outputTimes)):
        plot_global_comparison(config,
                               lonTarg,
                               latTarg,
                               modelOutput[i, :, :],
                               observations[i, :, :],
                               bias[i, :,:],
                               cmapModelObs,
                               clevsModelObs,
                               cmapDiff,
                               clevsDiff,
                               fileout = "%s/%s_%s_%s_years%04d-%04d.png" % (plots_dir, fileOutLabel,
				                                casename, outputTimes[i], climo_yr1, climo_yr2),
                               title = "%s (%s, years %04d-%04d)" % (field.upper(), outputTimes[i], climo_yr1, climo_yr2),
                               modelTitle = "%s" % casename,
                               obsTitle = obsTitleLabel,
                               diffTitle = "Model-Observations",
                               cbarlabel = unitsLabel)
