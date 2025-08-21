# -*- coding: utf-8 -*-
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
#

import xarray as xr
import numpy as np
import os
import matplotlib.pyplot as plt

from geometric_features import FeatureCollection, read_feature_collection
from geometric_features.aggregation import get_aggregator_by_name

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.io.utility import build_config_full_path, \
    get_files_year_month, make_directories, decode_strings
from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf_with_fill
from mpas_analysis.shared.io.utility import get_region_mask
from mpas_analysis.shared.timekeeping.utility import days_to_datetime
from mpas_analysis.shared.climatology import compute_climatology
from mpas_analysis.shared.constants import constants
from mpas_analysis.shared.html import write_image_xml
from mpas_analysis.shared.plot import savefig, add_inset
from mpas_analysis.shared.regions.compute_region_masks_subtask import get_feature_list


class OceanRegionalProfiles(AnalysisTask):
    """
    Compute and plot vertical profiles of regionally analyzed data.  The
    mean and standard deviation of the data are computed over each region.
    The mean and std. dev. are computed in time (within the requested seasons)
    and this result is plotted as a vertical profile with shading showing
    variability.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, regionMasksTask, controlConfig=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  mpas_tools.config.MpasConfigParser
            Contains configuration options

        regionMasksTask : ``ComputeRegionMasks``
            A task for computing region masks

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(OceanRegionalProfiles, self).__init__(
            config=config,
            taskName='oceanRegionalProfiles',
            componentName='ocean',
            tags=['profiles', 'climatology'])

        self.combineSubtasks = dict()
        self.computeSubtasks = dict()
        self.masksSubtasks = dict()

        startYear = config.getint('climatology', 'startYear')
        endYear = config.getint('climatology', 'endYear')

        regionGroups = config.getexpression('oceanRegionalProfiles',
                                            'regionGroups')

        for regionGroup in regionGroups:
            regionGroupSection = 'profiles{}'.format(
                regionGroup.replace(' ', ''))

            fields = config.getexpression(regionGroupSection, 'fields')

            max_bottom_depth = config.getexpression(regionGroupSection, 'maxDepth')
            seasons = config.getexpression(regionGroupSection, 'seasons')

            regionNames = config.getexpression(regionGroupSection,
                                               'regionNames')

            if len(regionNames) == 0:
                return
            if 'all' in regionNames:
                aggregationFunction, prefix, date = get_aggregator_by_name(
                     regionGroup)
                date = date
                outFileSuffix = '{}{}'.format(prefix, date)
                geojsonFileName = \
                    get_region_mask(self.config,
                                    '{}.geojson'.format(outFileSuffix))
                regionNames = get_feature_list(geojsonFileName)

            self.add_region_group(regionMasksTask, regionGroup, regionNames,
                                  fields, startYear, endYear, max_bottom_depth,
                                  seasons)

            combineSubtask = \
                self.combineSubtasks[regionGroup][(startYear, endYear)]

            masksSubtask = self.masksSubtasks[regionGroup]

            timeSeriesName = regionGroup.replace(' ', '')

            for field in fields:
                for regionName in regionNames:
                    for season in seasons:
                        plotSubtask = PlotRegionalProfileTimeSeriesSubtask(
                            self, masksSubtask, season, regionName, field,
                            timeSeriesName, startYear, endYear, controlConfig)
                        plotSubtask.run_after(combineSubtask)
                        self.add_subtask(plotSubtask)

    def add_region_group(self, regionMasksTask, regionGroup, regionNames,
                         fields, startYear, endYear, max_bottom_depth=None,
                         seasons=None):
        """
        Add years to the profiles to compute

        Parameters
        ----------
        startYear : int
            The start year of the time series

        endYear : int
            The end year

        """
        if regionGroup in self.masksSubtasks:
            masksSubtask = self.masksSubtasks[regionGroup]
        else:
            masksSubtask = regionMasksTask.add_mask_subtask(regionGroup)
            self.masksSubtasks[regionGroup] = masksSubtask

        if regionGroup not in self.computeSubtasks:
            self.computeSubtasks[regionGroup] = dict()
        if regionGroup not in self.combineSubtasks:
            self.combineSubtasks[regionGroup] = dict()

        timeSeriesName = regionGroup.replace(' ', '')

        if seasons is None:
            seasons = []

        key = (startYear, endYear)
        years = range(startYear, endYear + 1)
        if key in self.combineSubtasks[regionGroup]:
            combineSubtask = self.combineSubtasks[regionGroup][key]
            # add any missing fields and seasons
            _update_fields(combineSubtask.fields, fields)
            combineSubtask.seasons = list(set(seasons + combineSubtask.seasons))
        else:
            combineSubtask = CombineRegionalProfileTimeSeriesSubtask(
                self, regionGroup, timeSeriesName, seasons, fields,
                startYears=years, endYears=years)
            self.combineSubtasks[regionGroup][key] = combineSubtask

        # run one subtask per year
        for year in years:
            key = (year, year)
            if key in self.computeSubtasks[regionGroup]:
                computeSubtask = self.computeSubtasks[regionGroup][key]
                _update_fields(computeSubtask.fields, fields)
                combineSubtask.run_after(computeSubtask)
            else:
                computeSubtask = ComputeRegionalProfileTimeSeriesSubtask(
                    self, masksSubtask, regionGroup, regionNames, fields,
                    startYear=year, endYear=year,
                    max_bottom_depth=max_bottom_depth)
                computeSubtask.run_after(masksSubtask)
                combineSubtask.run_after(computeSubtask)
                self.computeSubtasks[regionGroup][key] = computeSubtask


class ComputeRegionalProfileTimeSeriesSubtask(AnalysisTask):
    """
    Compute regional statistics on each layer and time point of a set of
    MPAS fields

    Attributes
    ----------
    parentTask : ``OceanRegionalProfiles``
        The main task of which this is a subtask

    startYear, endYear : int
        The beginning and end of the time series to compute

    max_bottom_depth : float
        The maximum bottom depth of cells to include in the profile statistics
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, masksSubtask, regionGroup, regionNames,
                 fields, startYear, endYear, max_bottom_depth):
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : mpas_analysis.ocean.OceanRegionalProfiles
            The main task of which this is a subtask

        masksSubtask : mpas_analysis.shared.regions.ComputeRegionMasksSubtask
            A task for computing region masks

        regionGroup : str
            The name of the region group for which the region masks are defined

        regionNames : list
            The list of region names to compute and plot

        fields : list
            A list of dictionaries defining the fields to compute profile
            time series for

        startYear, endYear : int
            The beginning and end of the time series to compute

        max_bottom_depth : float
            The maximum bottom depth of cells to include in the profile
            statistics
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        subtaskName = 'compute{}Profiles_{:04d}-{:04d}'.format(
            regionGroup.replace(' ', ''), startYear, endYear)
        # first, call the constructor from the base class (AnalysisTask)
        super(ComputeRegionalProfileTimeSeriesSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName=subtaskName)

        parentTask.add_subtask(self)
        self.masksSubtask = masksSubtask
        if 'all' in regionNames:
            regionNames = get_feature_list(self.masksSubtask.geojsonFileName)
        self.regionNames = regionNames
        self.fields = fields
        self.startYear = startYear
        self.endYear = endYear
        self.max_bottom_depth = max_bottom_depth

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError
            if timeSeriesStatsMonthly is not enabled in the MPAS run
        """
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

    def run_task(self):
        """
        Compute time series of regional profiles
        """
        # Authors
        # -------
        # Milena Veneziani, Mark Petersen, Phillip J. Wolfram, Xylar Asay-Davis

        self.logger.info("\nCompute time series of regional profiles...")

        startDate = '{:04d}-01-01_00:00:00'.format(self.startYear)
        endDate = '{:04d}-12-31_23:59:59'.format(self.endYear)

        timeSeriesName = self.masksSubtask.regionGroup.replace(' ', '')

        outputDirectory = '{}/{}/'.format(
            build_config_full_path(self.config, 'output',
                                   'timeseriesSubdirectory'),
            timeSeriesName)
        try:
            os.makedirs(outputDirectory)
        except OSError:
            pass

        outputFileName = '{}/regionalProfiles_{}_{:04d}-{:04d}.nc'.format(
            outputDirectory, timeSeriesName, self.startYear, self.endYear)

        inputFiles = sorted(self.historyStreams.readpath(
            'timeSeriesStatsMonthlyOutput', startDate=startDate,
            endDate=endDate, calendar=self.calendar))

        years, months = get_files_year_month(inputFiles,
                                             self.historyStreams,
                                             'timeSeriesStatsMonthlyOutput')

        variableList = [field['mpas'] for field in self.fields]

        outputExists = os.path.exists(outputFileName)
        outputValid = outputExists
        if outputExists:
            with open_mpas_dataset(fileName=outputFileName,
                                   calendar=self.calendar,
                                   timeVariableNames=None,
                                   variableList=None,
                                   startDate=startDate,
                                   endDate=endDate) as dsIn:

                for inIndex in range(dsIn.sizes['Time']):

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
        meshFilename = self.get_mesh_filename()

        dsMesh = xr.open_dataset(meshFilename)
        dsMesh = dsMesh.isel(Time=0)
        areaCell = dsMesh.areaCell
        landIceFloatingMask = dsMesh.landIceFloatingMask.isel(Time=0)
        landIceFloatingMask = -1*(landIceFloatingMask-1)

        nVertLevels = dsMesh.sizes['nVertLevels']

        vertIndex = \
            xr.DataArray.from_dict({'dims': ('nVertLevels',),
                                    'data': np.arange(nVertLevels)})

        vertMask = vertIndex < dsMesh.maxLevelCell
        if self.max_bottom_depth is not None:
            depthMask = dsMesh.bottomDepth < self.max_bottom_depth
            vertDepthMask = np.logical_and(vertMask, depthMask)
        else:
            vertDepthMask = vertMask

        # get region masks
        regionMaskFileName = self.masksSubtask.maskFileName
        dsRegionMask = xr.open_dataset(regionMaskFileName)

        # figure out the indices of the regions to plot
        regionNames = decode_strings(dsRegionMask.regionNames)

        regionIndices = []
        for regionToPlot in self.regionNames:
            for index, regionName in enumerate(regionNames):
                if regionToPlot == regionName:
                    regionIndices.append(index)
                    break

        # select only those regions we want to plot
        dsRegionMask = dsRegionMask.isel(nRegions=regionIndices)
        cellMasks = dsRegionMask.regionCellMasks
        regionNamesVar = dsRegionMask.regionNames

        totalArea = self._masked_area_sum(cellMasks, areaCell, vertDepthMask)

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
            for field in self.fields:
                variableName = field['mpas']
                prefix = field['prefix']
                self.logger.info('      {}'.format(field['titleName']))

                var_mpas = dsLocal[variableName]
                print('DIMS-------------')
                print(var_mpas.dims)
                print(landIceFloatingMask.dims)
                var_mpas_masked = var_mpas*landIceFloatingMask
                print(var_mpas_masked.dims)
                var = var_mpas_masked.where(vertDepthMask)
                print(var.dims)

                meanName = '{}_mean'.format(prefix)
                dsLocal[meanName] = \
                    self._masked_area_sum(cellMasks, areaCell, var) / totalArea

                meanSquaredName = '{}_meanSquared'.format(prefix)
                dsLocal[meanSquaredName] = \
                    self._masked_area_sum(cellMasks, areaCell, var**2) / \
                    totalArea

            # drop the original variables
            dsLocal = dsLocal.drop_vars(variableList)

            datasets.append(dsLocal)

        # combine data sets into a single data set
        dsOut = xr.concat(datasets, 'Time')

        dsOut.coords['regionNames'] = regionNamesVar
        dsOut['totalArea'] = totalArea
        dsOut.coords['year'] = (('Time',), years)
        dsOut['year'].attrs['units'] = 'years'
        dsOut.coords['month'] = (('Time',), months)
        dsOut['month'].attrs['units'] = 'months'

        # Note: restart file, not a mesh file because we need refBottomDepth,
        # not in a mesh file
        meshFilename = self.get_mesh_filename()

        with xr.open_dataset(meshFilename) as dsMesh:
            depths = dsMesh.refBottomDepth.values
            z = np.zeros(depths.shape)
            z[0] = -0.5 * depths[0]
            z[1:] = -0.5 * (depths[0:-1] + depths[1:])

        dsOut.coords['z'] = (('nVertLevels',), z)
        dsOut['z'].attrs['units'] = 'meters'

        write_netcdf_with_fill(dsOut, outputFileName)

    @staticmethod
    def _masked_area_sum(cellMasks, areaCell, var):
        """sum a variable over the masked areas"""
        nRegions = cellMasks.sizes['nRegions']
        totals = []
        for index in range(nRegions):
            mask = cellMasks.isel(nRegions=slice(index, index+1))
            totals.append((mask * areaCell * var).sum('nCells'))

        total = xr.concat(totals, 'nRegions')
        return total


class CombineRegionalProfileTimeSeriesSubtask(AnalysisTask):
    """
    Combine individual time series into a single data set
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, regionGroup, timeSeriesName, seasons, fields,
                 startYears, endYears):
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``OceanRegionalProfiles``
            The main task of which this is a subtask

        regionGroup : str
            The name of the region group for which the region masks are defined

        seasons : list
            A list of seasons to compute statistic on

        fields : list
            A list of dictionaries defining the fields to compute profile
            time series for

        startYears, endYears : list
            The beginning and end of each time series to combine
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        subtaskName = 'combine{}Profiles_{:04d}-{:04d}'.format(
            regionGroup.replace(' ', ''), startYears[0], endYears[-1])
        # first, call the constructor from the base class (AnalysisTask)
        super(CombineRegionalProfileTimeSeriesSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName=subtaskName)

        parentTask.add_subtask(self)
        self.startYears = startYears
        self.endYears = endYears
        self.timeSeriesName = timeSeriesName
        self.seasons = seasons
        self.fields = fields

    def run_task(self):
        """
        Combine the time series
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        timeSeriesName = self.timeSeriesName

        outputDirectory = '{}/{}/'.format(
            build_config_full_path(self.config, 'output',
                                   'timeseriesSubdirectory'),
            timeSeriesName)

        outputFileName = '{}/regionalProfiles_{}_{:04d}-{:04d}.nc'.format(
            outputDirectory, timeSeriesName, self.startYears[0],
            self.endYears[-1])

        useExisting = False
        ds = None
        if os.path.exists(outputFileName):
            ds = xr.open_dataset(outputFileName, decode_times=False)
            if ds.sizes['Time'] > 0:
                useExisting = True
            else:
                ds.close()

        if not useExisting:

            inFileNames = []
            for startYear, endYear in zip(self.startYears, self.endYears):
                inFileName = '{}/regionalProfiles_{}_{:04d}-{:04d}.nc'.format(
                    outputDirectory, timeSeriesName, startYear, endYear)
                inFileNames.append(inFileName)

            ds = xr.open_mfdataset(inFileNames, combine='nested',
                                   concat_dim='Time', decode_times=False)

            ds.load()

            ds['totalArea'] = ds['totalArea'].isel(Time=0)

            write_netcdf_with_fill(ds, outputFileName)

        regionNames = ds['regionNames']
        ds = ds.drop_vars('regionNames')

        profileMask = ds['totalArea'] > 0

        outputDirectory = build_config_full_path(self.config, 'output',
                                                 'profilesSubdirectory')

        make_directories(outputDirectory)

        for season in self.seasons:
            outputFileName = \
                '{}/{}_{}_{:04d}-{:04d}.nc'.format(
                    outputDirectory, timeSeriesName, season,
                    self.startYears[0], self.endYears[-1])
            if not os.path.exists(outputFileName):
                monthValues = constants.monthDictionary[season]
                dsSeason = compute_climatology(ds, monthValues,
                                               calendar=self.calendar,
                                               maskVaries=False)

                for field in self.fields:
                    prefix = field['prefix']

                    mean = dsSeason['{}_mean'.format(prefix)].where(
                        profileMask)
                    meanSquared = \
                        dsSeason['{}_meanSquared'.format(prefix)].where(
                            profileMask)
                    stdName = '{}_std'.format(prefix)

                    dsSeason[stdName] = np.sqrt(meanSquared - mean**2).where(
                        profileMask)
                    dsSeason['{}_mean'.format(prefix)] = mean

                dsSeason.coords['regionNames'] = regionNames
                write_netcdf_with_fill(dsSeason, outputFileName)


class PlotRegionalProfileTimeSeriesSubtask(AnalysisTask):
    """
    Plot a profile averaged over an ocean region and in time, along with
    variability in both space and time.

    Attributes
    ----------
    parentTask : ``AnalysisTask``
        The parent task of which this is a subtask

    season : str
        The season being plotted

    regionName : str
        The region being plotted

    field : dict
        Information about the field (e.g. temperature) being plotted

    controlconfig : mpas_tools.config.MpasConfigParser
        Configuration options for a control run (if any)
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, masksSubtask, season, regionName, field,
                 timeSeriesName, startYear, endYear, controlConfig):

        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : OceanRegionalProfiles
            The parent task of which this is a subtask

        masksSubtask : mpas_analysis.shared.regions.ComputeRegionMasksSubtask
            A task for computing region masks

        season : str
            The season being plotted

        regionName : str
            The region being plotted

        field : dict
            Information about the field (e.g. temperature) being plotted

        timeSeriesName : str
            The name of the time series, related to the name of the region
            group but appropriate for a file prefix or suffix

        startYear, endYear : int
            The beginning and end of the time series to compute

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        subtaskName = '{}_{}_{}'.format(field['prefix'],
                                        regionName.replace(' ', '_'),
                                        season)
        # first, call the constructor from the base class (AnalysisTask)
        super(PlotRegionalProfileTimeSeriesSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName=subtaskName)

        self.controlConfig = controlConfig

        self.masksSubtask = masksSubtask
        self.timeSeriesName = timeSeriesName
        self.startYear = startYear
        self.endYear = endYear
        self.season = season
        self.regionName = regionName
        self.field = field
        self.filePrefix = None
        self.xmlFileNames = []

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(PlotRegionalProfileTimeSeriesSubtask, self).setup_and_check()

        self.filePrefix = 'regionalProfile_{}_{}_{}_years{:04d}-{:04d}'.format(
            self.field['prefix'], self.regionName.replace(' ', '_'),
            self.season, self.startYear,
            self.endYear)
        self.xmlFileNames = ['{}/{}.xml'.format(self.plotsDirectory,
                                                self.filePrefix)]

    def run_task(self):
        """
        Plot a depth profile with variability
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config
        startYear = self.startYear
        endYear = self.endYear

        regionMaskFile = self.masksSubtask.geojsonFileName

        fcAll = read_feature_collection(regionMaskFile)

        fc = FeatureCollection()
        for feature in fcAll.features:
            if feature['properties']['name'] == self.regionName:
                fc.add_feature(feature)
                break

        inDirectory = build_config_full_path(config, 'output',
                                             'profilesSubdirectory')
        timeSeriesName = self.timeSeriesName
        inFileName = '{}/{}_{}_{:04d}-{:04d}.nc'.format(
            inDirectory, timeSeriesName, self.season,
            self.startYear, self.endYear)

        regionGroup = self.masksSubtask.regionGroup
        regionGroupSection = 'profiles{}'.format(
            regionGroup.replace(' ', ''))

        ds = xr.open_dataset(inFileName)
        allRegionNames = decode_strings(ds.regionNames)

        regionIndex = allRegionNames.index(self.regionName)
        ds = ds.isel(nRegions=regionIndex)
        meanFieldName = '{}_mean'.format(self.field['prefix'])
        stdFieldName = '{}_std'.format(self.field['prefix'])

        mainRunName = config.get('runs', 'mainRunName')
        profileGalleryGroup = config.get(regionGroupSection,
                                         'profileGalleryGroup')

        titleFieldName = self.field['titleName']
        regionName = self.regionName.replace('_', ' ')

        xLabel = '{} ({})'.format(titleFieldName, self.field['units'])
        yLabel = 'depth (m)'
        outFileName = '{}/{}.png'.format(self.plotsDirectory, self.filePrefix)
        lineColors = ['k']
        lineWidths = [1.6]
        zArrays = [ds.z.values]
        fieldArrays = [ds[meanFieldName].values]
        errArrays = [ds[stdFieldName].values]
        if self.controlConfig is None:
            title = '{} {}, years {:04d}-{:04d}\n{}'.format(
                regionName, self.season, startYear, endYear, mainRunName)
            legendText = [None]
        else:
            controlStartYear = self.controlConfig.getint('climatology',
                                                         'startYear')
            controlEndYear = self.controlConfig.getint('climatology',
                                                       'endYear')
            controlRunName = self.controlConfig.get('runs', 'mainRunName')

            if controlStartYear == startYear and controlEndYear == endYear:
                title = '{} {}, years {:04d}-{:04d}'.format(
                    regionName, self.season, startYear, endYear)
                legendText = [mainRunName, controlRunName]
            elif mainRunName == controlRunName:
                title = '{} {}\n{}'.format(
                    regionName, self.season, mainRunName)
                legendText = ['{:04d}-{:04d}'.format(startYear, endYear),
                              '{:04d}-{:04d}'.format(controlStartYear,
                                                     controlEndYear)]
            else:
                title = '{} {}   '.format(regionName, self.season)
                legendText = ['{} {:04d}-{:04d}'.format(mainRunName, startYear,
                                                        endYear),
                              '{} {:04d}-{:04d}'.format(controlRunName,
                                                        controlStartYear,
                                                        controlEndYear)]

            controlDirectory = build_config_full_path(
                self.controlConfig, 'output',
                'profilesSubdirectory')

            controlFileName = \
                '{}/{}_{}_{:04d}-{:04d}.nc'.format(
                    controlDirectory, timeSeriesName, self.season,
                    controlStartYear, controlEndYear)

            dsControl = xr.open_dataset(controlFileName)
            allRegionNames = decode_strings(dsControl.regionNames)
            regionIndex = allRegionNames.index(self.regionName)
            dsControl = dsControl.isel(nRegions=regionIndex)

            lineColors.append('r')
            lineWidths.append(1.2)
            zArrays.append(dsControl.z.values)
            fieldArrays.append(dsControl[meanFieldName].values)
            errArrays.append(dsControl[stdFieldName].values)

        depthRange = config.getexpression(regionGroupSection, 'depthRange')
        if len(depthRange) == 0:
            depthRange = None

        fig = self.plot(zArrays, fieldArrays, errArrays,
                        lineColors=lineColors, lineWidths=lineWidths,
                        legendText=legendText, title=title, xLabel=xLabel,
                        yLabel=yLabel, yLim=depthRange)

        # do this before the inset because otherwise it moves the inset
        # and cartopy doesn't play too well with tight_layout anyway
        plt.tight_layout()

        add_inset(fig, fc, width=1.0, height=1.0)

        savefig(outFileName, config, tight=False)

        caption = '{} {} vs depth'.format(regionName, titleFieldName)
        write_image_xml(
            config=config,
            filePrefix=self.filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup=profileGalleryGroup,
            groupLink='ocnregprofs',
            imageDescription=caption,
            imageCaption=caption,
            gallery=titleFieldName,
            thumbnailDescription='{} {}'.format(regionName, self.season))

    def plot(self, zArrays, fieldArrays, errArrays, lineColors, lineWidths,
             legendText, title, xLabel, yLabel, xLim=None, yLim=None,
             figureSize=(10, 4), dpi=None):
        """
        Plots a 1D line plot with error bars if available.

        Parameters
        ----------
        zArrays : list of float arrays
            x array (latitude, or any other x axis except time)

        fieldArrays : list of float arrays
            y array (any field as function of x)

        errArrays : list of float arrays
            error array (y errors)

        lineColors, legendText : list of str
            control line color and corresponding legend text.

        lineWidths : list of float
            control line width

        title : str
            title of plot

        xLabel, yLabel : str
            label of x- and y-axis

        xLim : float array, optional
            x range of plot

        yLim : float array, optional
            y range of plot

        figureSize : tuple of float, optional
            size of the figure in inches

        dpi : int, optional
            the number of dots per inch of the figure, taken from section `
            `plot`` option ``dpi`` in the config file by default
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config

        # set up figure
        if dpi is None:
            dpi = config.getint('plot', 'dpi')
        fig = plt.figure(figsize=figureSize, dpi=dpi)

        plotLegend = False
        for dsIndex in range(len(zArrays)):
            zArray = zArrays[dsIndex]
            fieldArray = fieldArrays[dsIndex]
            errArray = errArrays[dsIndex]
            if zArray is None:
                continue

            if legendText is None:
                label = None
            else:
                label = legendText[dsIndex]
                plotLegend = True
            if lineColors is None:
                color = 'k'
            else:
                color = lineColors[dsIndex]
            if lineWidths is None:
                linewidth = 1.
            else:
                linewidth = lineWidths[dsIndex]

            plt.plot(fieldArray, zArray, color=color, linewidth=linewidth,
                     label=label)
            if errArray is not None:
                plt.fill_betweenx(zArray, fieldArray, fieldArray + errArray,
                                  facecolor=color, alpha=0.2)
                plt.fill_betweenx(zArray, fieldArray, fieldArray - errArray,
                                  facecolor=color, alpha=0.2)

        if plotLegend and len(zArrays) > 1:
            plt.legend(loc='lower left')

        axis_font = {'size': config.get('plot', 'axisFontSize')}
        title_font = {'size': config.get('plot', 'titleFontSize'),
                      'color': config.get('plot', 'titleFontColor'),
                      'weight': config.get('plot', 'titleFontWeight')}
        if title is not None:
            plt.title(title, **title_font)
        if xLabel is not None:
            plt.xlabel(xLabel, **axis_font)
        if yLabel is not None:
            plt.ylabel(yLabel, **axis_font)

        if xLim:
            plt.xlim(xLim)
        if yLim:
            plt.ylim(yLim)
        return fig


def _update_fields(fields, newFields):
    for outer in range(len(newFields)):
        found = False
        for inner in range(len(fields)):
            if fields[inner]['prefix'] == newFields[outer]['prefix']:
                for item in ['mpas', 'units', 'titleName']:
                    if fields[inner][item] != newFields[outer][item]:
                        raise ValueError(
                            'item {} in fields is not consistent between '
                            'profiles and Hovmoller tasks, which will have '
                            'unexpected consequences')
                found = True
                break
        if not found:
            fields.append(newFields[outer])
