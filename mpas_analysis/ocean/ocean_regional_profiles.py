# -*- coding: utf-8 -*-
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2018 Los Alamos National Security, LLC. All rights reserved.
# Copyright (c) 2018 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2018 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
#

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr
import numpy as np
import os

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.io.utility import build_config_full_path, \
    get_files_year_month, make_directories
from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf
from mpas_analysis.shared.timekeeping.utility import days_to_datetime
from mpas_analysis.shared.regions import ComputeRegionMasksSubtask, \
    get_feature_list
from mpas_analysis.shared.climatology import compute_climatology
from mpas_analysis.shared.constants import constants


class OceanRegionalProfiles(AnalysisTask):  # {{{
    '''
    Compute and plot vertical profiles of regionally analyzed data.  The
    mean and standard deviation of the data are computed over each region.
    The mean isdisplayed as a Hovmoller plot.  The mean and std. dev. are
    further computed in time (within the requested seasons) and this result
    is plotted as a vertical profile with shading showing variability.
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, controlConfig=None):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(OceanRegionalProfiles, self).__init__(
            config=config,
            taskName='oceanRegionalProfiles',
            componentName='ocean',
            tags=['profiles', 'climatology'])

        startYear = config.getint('climatology', 'startYear')
        endYear = config.getint('climatology', 'endYear')

        self.fields = config.getExpression('oceanRegionalProfiles', 'fields')

        self.seasons = config.getExpression('oceanRegionalProfiles', 'seasons')

        self.regionMaskSuffix = config.get('oceanRegionalProfiles',
                                           'regionMaskSuffix')

        self.regionNames = config.getExpression('oceanRegionalProfiles',
                                                'regionNames')

        regionMaskDirectory = build_config_full_path(config,
                                                     'diagnostics',
                                                     'regionMaskSubdirectory')
        masksFile = '{}/{}.geojson'.format(regionMaskDirectory,
                                           self.regionMaskSuffix)

        masksSubtask = ComputeRegionMasksSubtask(
                self, masksFile,
                outFileSuffix=self.regionMaskSuffix,
                featureList=self.regionNames)

        if 'all' in self.regionNames:
            self.regionNames = get_feature_list(config, masksFile)

        self.masksSubtask = masksSubtask

        years = range(startYear, endYear+1)

        # in the end, we'll combine all the time series into one, but we create
        # this task first so it's easier to tell it to run after all the
        # compute tasks
        combineSubtask = CombineRegionalProfileTimeSeriesSubtask(
                self, startYears=years, endYears=years)

        # run one subtask per year
        for year in years:
            computeSubtask = ComputeRegionalProfileTimeSeriesSubtask(
                    self, startYear=year, endYear=year)
            computeSubtask.run_after(masksSubtask)
            combineSubtask.run_after(computeSubtask)

        # plotSubtask = PlotRegionalProfileTimeSeriesSubtask(
        #         self, controlConfig)
        # plotSubtask.run_after(combineSubtask)

        # }}}
    # }}}


class ComputeRegionalProfileTimeSeriesSubtask(AnalysisTask):  # {{{
    '''
    Compute regional statistics on each layer and time point of a set of
    MPAS fields

    Attributes
    ----------
    parentTask : ``OceanRegionalProfiles``
        The main task of which this is a subtask

    startyear, endYear : int
        The beginning and end of the time series to compute
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, startYear, endYear):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``OceanRegionalProfiles``
            The main task of which this is a subtask

        startyear, endYear : int
            The beginning and end of the time series to compute
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(ComputeRegionalProfileTimeSeriesSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='computeRegionalProfileTimeSeries_{:04d}-{:04d}'
                        ''.format(startYear, endYear))

        parentTask.add_subtask(self)
        self.parentTask = parentTask
        self.startYear = startYear
        self.endYear = endYear
        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError
            if timeSeriesStatsMonthly is not enabled in the MPAS run
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ComputeRegionalProfileTimeSeriesSubtask, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        # }}}

    def run_task(self):  # {{{
        '''
        Process MOC analysis member data if available, or compute MOC at
        post-processing if not.
        '''
        # Authors
        # -------
        # Milena Veneziani, Mark Petersen, Phillip J. Wolfram, Xylar Asay-Davis

        self.logger.info("\nCompute time series of regional profiles...")

        startDate = '{:04d}-01-01_00:00:00'.format(self.startYear)
        endDate = '{:04d}-12-31_23:59:59'.format(self.endYear)

        timeSeriesName = self.parentTask.regionMaskSuffix

        outputDirectory = '{}/{}/'.format(
                build_config_full_path(self.config, 'output',
                                       'timeseriesSubdirectory'),
                timeSeriesName)
        try:
            os.makedirs(outputDirectory)
        except OSError:
            pass

        outputFileName = '{}/{}_{:04d}-{:04d}.nc'.format(
                outputDirectory, timeSeriesName, self.startYear, self.endYear)

        inputFiles = sorted(self.historyStreams.readpath(
                'timeSeriesStatsMonthlyOutput', startDate=startDate,
                endDate=endDate, calendar=self.calendar))

        years, months = get_files_year_month(inputFiles,
                                             self.historyStreams,
                                             'timeSeriesStatsMonthlyOutput')

        variableList = [field['mpas'] for field in self.parentTask.fields]

        outputExists = os.path.exists(outputFileName)
        outputValid = outputExists
        if outputExists:
            with open_mpas_dataset(fileName=outputFileName,
                                   calendar=self.calendar,
                                   timeVariableNames=None,
                                   variableList=None,
                                   startDate=startDate,
                                   endDate=endDate) as dsIn:

                for inIndex in range(dsIn.dims['Time']):

                    mask = np.logical_and(
                            dsIn.year[inIndex].values == years,
                            dsIn.month[inIndex].values == months)
                    if np.count_nonzero(mask) == 0:
                        outputValid = False
                        break

        if outputValid:
            self.logger.info('  Time series exists -- Done.')
            return

        # get areaCell
        restartFileName = \
            self.runStreams.readpath('restart')[0]

        dsRestart = xr.open_dataset(restartFileName)
        dsRestart = dsRestart.isel(Time=0)
        areaCell = dsRestart.areaCell

        nVertLevels = dsRestart.sizes['nVertLevels']

        vertIndex = \
            xr.DataArray.from_dict({'dims': ('nVertLevels',),
                                    'data': np.arange(nVertLevels)})

        vertMask = vertIndex < dsRestart.maxLevelCell

        # get region masks
        regionMaskFileName = self.parentTask.masksSubtask.maskFileName
        dsRegionMask = xr.open_dataset(regionMaskFileName)

        # figure out the indices of the regions to plot
        regionNames = [bytes.decode(name) for name in
                       dsRegionMask.regionNames.values]
        regionIndices = []
        for regionToPlot in self.parentTask.regionNames:
            for index, regionName in enumerate(regionNames):
                if regionToPlot == regionName:
                    regionIndices.append(index)
                    break

        # select only those regions we want to plot
        dsRegionMask = dsRegionMask.isel(nRegions=regionIndices)
        cellMasks = dsRegionMask.regionCellMasks

        totalArea = (cellMasks*areaCell*vertMask).sum('nCells')

        datasets = []
        for timeIndex, fileName in enumerate(inputFiles):

            dsLocal = open_mpas_dataset(
                fileName=fileName,
                calendar=self.calendar,
                variableList=variableList,
                startDate=startDate,
                endDate=endDate)
            dsLocal = dsLocal.isel(Time=0)
            time = dsLocal.Time.values
            date = days_to_datetime(time, calendar=self.calendar)

            self.logger.info('    date: {:04d}-{:02d}'.format(date.year,
                                                              date.month))

            # for each region and variable, compute area-weighted sum and
            # squared sum
            for field in self.parentTask.fields:
                variableName = field['mpas']
                prefix = field['prefix']
                self.logger.info('      {}'.format(field['titleName']))

                var = dsLocal[variableName].where(vertMask)

                meanName = '{}_mean'.format(prefix)
                dsLocal[meanName] = \
                    (cellMasks*areaCell*var).sum('nCells')/totalArea

                meanSquaredName = '{}_meanSquared'.format(prefix)
                dsLocal[meanSquaredName] = \
                    (cellMasks*areaCell*var**2).sum('nCells')/totalArea

            # drop the original variables
            dsLocal = dsLocal.drop(variableList)

            datasets.append(dsLocal)

        # combine data sets into a single data set
        dsOut = xr.concat(datasets, 'Time')

        dsOut['totalArea'] = totalArea
        dsOut['year'] = (('Time'), years)
        dsOut['year'].attrs['units'] = 'years'
        dsOut['month'] = (('Time'), months)
        dsOut['month'].attrs['units'] = 'months'
        write_netcdf(dsOut, outputFileName)
        # }}}
    # }}}


class CombineRegionalProfileTimeSeriesSubtask(AnalysisTask):  # {{{
    '''
    Combine individual time series into a single data set
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, startYears, endYears):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``StreamfunctionMOC``
            The main task of which this is a subtask

        startyear, endYear : list of int
            The beginning and end of each time series to combine
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(CombineRegionalProfileTimeSeriesSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='combineRegionalProfileTimeSeries')

        parentTask.add_subtask(self)
        self.parentTask = parentTask
        self.startYears = startYears
        self.endYears = endYears
        # }}}

    def run_task(self):  # {{{
        '''
        Combine the time series
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        timeSeriesName = self.parentTask.regionMaskSuffix

        outputDirectory = '{}/{}/'.format(
                build_config_full_path(self.config, 'output',
                                       'timeseriesSubdirectory'),
                timeSeriesName)

        outputFileName = '{}/{}_{:04d}-{:04d}.nc'.format(
                outputDirectory, timeSeriesName, self.startYears[0],
                self.endYears[-1])

        if os.path.exists(outputFileName):
            ds = xr.open_dataset(outputFileName, decode_times=False)
        else:

            inFileNames = []
            for startYear, endYear in zip(self.startYears, self.endYears):
                inFileName = '{}/{}_{:04d}-{:04d}.nc'.format(
                        outputDirectory, timeSeriesName, startYear, endYear)
                inFileNames.append(inFileName)

            ds = xr.open_mfdataset(inFileNames, concat_dim='Time',
                                   decode_times=False)

            ds['totalArea'] = ds['totalArea'].isel(Time=0)

            write_netcdf(ds, outputFileName)

        outputDirectory = build_config_full_path(self.config, 'output',
                                                 'profilesSubdirectory')

        make_directories(outputDirectory)

        for season in self.parentTask.seasons:
            outputFileName = '{}/{}_{}_{:04d}-{:04d}.nc'.format(
                    outputDirectory, timeSeriesName, season,
                    self.startYears[0], self.endYears[-1])
            if not os.path.exists(outputFileName):
                monthValues = constants.monthDictionary[season]
                dsSeason = compute_climatology(ds, monthValues,
                                               calendar=self.calendar,
                                               maskVaries=False)

                for field in self.parentTask.fields:
                    prefix = field['prefix']

                    mean = dsSeason['{}_mean'.format(prefix)]
                    meanSquared = dsSeason['{}_meanSquared'.format(prefix)]
                    stdName = '{}_std'.format(prefix)

                    dsSeason[stdName] = np.sqrt(meanSquared - mean**2)

                write_netcdf(dsSeason, outputFileName)

        # }}}
    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
