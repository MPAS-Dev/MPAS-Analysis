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
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os
import xarray
import numpy
import matplotlib.pyplot as plt

from geometric_features import FeatureCollection, read_feature_collection

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.plot import timeseries_analysis_plot, savefig, \
    add_inset

from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf

from mpas_analysis.shared.io.utility import build_config_full_path, \
    build_obs_path, get_files_year_month, decode_strings, get_region_mask

from mpas_analysis.shared.html import write_image_xml

from mpas_analysis.ocean.utility import compute_zmid

from mpas_analysis.shared.constants import constants


class TimeSeriesOceanRegions(AnalysisTask):  # {{{
    """
    Performs analysis of the time-series output of regionoal mean temperature,
    salinity, etc.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, regionMasksTask, controlConfig=None):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        regionMasksTask : ``ComputeRegionMasks``
            A task for computing region masks

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(TimeSeriesOceanRegions, self).__init__(
            config=config,
            taskName='timeSeriesOceanRegions',
            componentName='ocean',
            tags=['timeSeries', 'regions', 'antarctic'])

        startYear = config.getint('timeSeries', 'startYear')
        endYear = config.get('timeSeries', 'endYear')
        if endYear == 'end':
            # a valid end year wasn't found, so likely the run was not found,
            # perhaps because we're just listing analysis tasks
            endYear = startYear
        else:
            endYear = int(endYear)

        regionGroups = config.getExpression(self.taskName, 'regionGroups')

        obsDicts = {
            'SOSE': {
                'suffix': 'SOSE',
                'gridName': 'SouthernOcean_0.167x0.167degree',
                'gridFileName': 'SOSE/SOSE_2005-2010_monthly_pot_temp_'
                                'SouthernOcean_0.167x0.167degree_20180710.nc',
                'TFileName': 'SOSE/SOSE_2005-2010_monthly_pot_temp_'
                             'SouthernOcean_0.167x0.167degree_20180710.nc',
                'SFileName': 'SOSE/SOSE_2005-2010_monthly_salinity_'
                             'SouthernOcean_0.167x0.167degree_20180710.nc',
                'volFileName': 'SOSE/SOSE_volume_'
                               'SouthernOcean_0.167x0.167degree_20190815.nc',
                'lonVar': 'lon',
                'latVar': 'lat',
                'TVar': 'theta',
                'SVar': 'salinity',
                'volVar': 'volume',
                'zVar': 'z',
                'tDim': 'Time',
                'legend': 'SOSE 2005-2010 ANN mean'},
            'WOA18': {
                'suffix': 'WOA18',
                'gridName': 'Global_0.25x0.25degree',
                'gridFileName': 'WOA18/woa18_decav_04_TS_mon_20190829.nc',
                'TFileName': 'WOA18/woa18_decav_04_TS_mon_20190829.nc',
                'SFileName': 'WOA18/woa18_decav_04_TS_mon_20190829.nc',
                'volFileName': None,
                'lonVar': 'lon',
                'latVar': 'lat',
                'TVar': 't_an',
                'SVar': 's_an',
                'volVar': 'volume',
                'zVar': 'depth',
                'tDim': 'month',
                'legend': 'WOA18 1955-2017 ANN mean'}}

        for regionGroup in regionGroups:
            sectionSuffix = regionGroup[0].upper() + \
                regionGroup[1:].replace(' ', '')
            sectionName = 'timeSeries{}'.format(sectionSuffix)

            regionNames = config.getExpression(sectionName, 'regionNames')
            if len(regionNames) == 0:
                # no regions in this group were requested
                continue

            masksSubtask = regionMasksTask.add_mask_subtask(
                regionGroup=regionGroup)

            regionNames = masksSubtask.expand_region_names(regionNames)

            years = list(range(startYear, endYear + 1))

            obsList = config.getExpression(sectionName, 'obs')
            groupObsDicts = {}

            for obsName in obsList:
                localObsDict = dict(obsDicts[obsName])
                obsFileName = build_obs_path(
                    config, component=self.componentName,
                    relativePath=localObsDict['gridFileName'])
                obsMasksSubtask = regionMasksTask.add_mask_subtask(
                    regionGroup=regionGroup, obsFileName=obsFileName,
                    lonVar=localObsDict['lonVar'],
                    latVar=localObsDict['latVar'],
                    meshName=localObsDict['gridName'])

                obsDicts[obsName]['maskTask'] = obsMasksSubtask

                localObsDict['maskTask'] = obsMasksSubtask
                groupObsDicts[obsName] = localObsDict

            # in the end, we'll combine all the time series into one, but we
            # create this task first so it's easier to tell it to run after all
            # the compute tasks
            combineSubtask = CombineRegionalProfileTimeSeriesSubtask(
                self, startYears=years, endYears=years,
                regionGroup=regionGroup)

            depthMasksSubtask = ComputeRegionDepthMasksSubtask(
                self,  masksSubtask=masksSubtask, regionGroup=regionGroup,
                regionNames=regionNames)
            depthMasksSubtask.run_after(masksSubtask)

            # run one subtask per year
            for year in years:
                computeSubtask = ComputeRegionTimeSeriesSubtask(
                    self, startYear=year, endYear=year,
                    masksSubtask=masksSubtask, regionGroup=regionGroup,
                    regionNames=regionNames)
                self.add_subtask(computeSubtask)
                computeSubtask.run_after(depthMasksSubtask)
                computeSubtask.run_after(masksSubtask)
                combineSubtask.run_after(computeSubtask)

            self.add_subtask(combineSubtask)

            for index, regionName in enumerate(regionNames):

                fullSuffix = sectionSuffix + '_' + regionName.replace(' ', '')

                obsSubtasks = {}
                for obsName in obsList:
                    localObsDict = dict(groupObsDicts[obsName])

                    obsSubtask = ComputeObsRegionalTimeSeriesSubtask(
                        self, regionGroup, regionName, fullSuffix, localObsDict)
                    obsSubtasks[obsName] = obsSubtask

                plotRegionSubtask = PlotRegionTimeSeriesSubtask(
                    self, regionGroup, regionName, index, controlConfig,
                    sectionName, fullSuffix, obsSubtasks,
                    masksSubtask.geojsonFileName)
                plotRegionSubtask.run_after(combineSubtask)
                self.add_subtask(plotRegionSubtask)

        # }}}

    # }}}


class ComputeRegionDepthMasksSubtask(AnalysisTask):  # {{{
    """
    Compute masks for regional and depth mean

    Attributes
    ----------
    masksSubtask : ``ComputeRegionMasksSubtask``
        A task for creating mask files for each region to plot

    regionGroup : str
        The name of the region group being computed (e.g. "Antarctic Basins")

    regionNames : list of str
        The names of the regions to compute
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, masksSubtask, regionGroup, regionNames):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``TimeSeriesOceanRegions``
            The main task of which this is a subtask

        masksSubtask : ``ComputeRegionMasksSubtask``
            A task for creating mask files for each region to plot

        regionGroup : str
            The name of the region group being computed (e.g. "Antarctic
            Basins")

        regionNames : list of str
            The names of the regions to compute
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        suffix = regionGroup[0].upper() + regionGroup[1:].replace(' ', '')

        # first, call the constructor from the base class (AnalysisTask)
        super(ComputeRegionDepthMasksSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='computeDepthMask{}'.format(suffix))

        parentTask.add_subtask(self)
        self.masksSubtask = masksSubtask
        self.regionGroup = regionGroup
        self.regionNames = regionNames
        # }}}

    def run_task(self):  # {{{
        """
        Compute the regional-mean time series
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config

        self.logger.info("\nCompute depth mask for regional means...")

        regionGroup = self.regionGroup
        sectionSuffix = regionGroup[0].upper() + \
            regionGroup[1:].replace(' ', '')
        timeSeriesName = sectionSuffix
        sectionName = 'timeSeries{}'.format(sectionSuffix)

        outputDirectory = '{}/{}/'.format(
            build_config_full_path(config, 'output', 'timeseriesSubdirectory'),
            timeSeriesName)
        try:
            os.makedirs(outputDirectory)
        except OSError:
            pass

        outFileName = '{}/depthMasks_{}.nc'.format(outputDirectory,
                                                   timeSeriesName)

        if os.path.exists(outFileName):
            self.logger.info('  Mask file exists -- Done.')
            return

        # Load mesh related variables
        try:
            restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for ocean region time series')

        if config.has_option(sectionName, 'zmin'):
            config_zmin = config.getfloat(sectionName, 'zmin')
        else:
            config_zmin = None

        if config.has_option(sectionName, 'zmax'):
            config_zmax = config.getfloat(sectionName, 'zmax')
        else:
            config_zmax = None

        dsRestart = xarray.open_dataset(restartFileName).isel(Time=0)
        zMid = compute_zmid(dsRestart.bottomDepth, dsRestart.maxLevelCell,
                            dsRestart.layerThickness)
        areaCell = dsRestart.areaCell
        if 'landIceMask' in dsRestart:
            # only the region outside of ice-shelf cavities
            openOceanMask = dsRestart.landIceMask == 0
        else:
            openOceanMask = None

        regionMaskFileName = self.masksSubtask.maskFileName
        dsRegionMask = xarray.open_dataset(regionMaskFileName)
        maskRegionNames = decode_strings(dsRegionMask.regionNames)

        regionIndices = []
        for regionName in self.regionNames:
            for index, otherName in enumerate(maskRegionNames):
                if regionName == otherName:
                    regionIndices.append(index)
                    break

        # select only those regions we want to plot
        dsRegionMask = dsRegionMask.isel(nRegions=regionIndices)

        nRegions = dsRegionMask.sizes['nRegions']

        datasets = []
        for regionIndex in range(nRegions):
            self.logger.info('    region: {}'.format(
                self.regionNames[regionIndex]))
            dsRegion = dsRegionMask.isel(nRegions=regionIndex)
            cellMask = dsRegion.regionCellMasks == 1

            if openOceanMask is not None:
                cellMask = numpy.logical_and(cellMask, openOceanMask)

            totalArea = areaCell.where(cellMask).sum()
            self.logger.info('      totalArea: {} mil. km^2'.format(
                1e-12 * totalArea.values))

            if config_zmin is None:
                if 'zminRegions' in dsRegion:
                    zmin = dsRegion.zminRegions.values
                else:
                    # the old naming convention, used in some pre-generated
                    # mask files
                    zmin = dsRegion.zmin.values
            else:
                zmin = config_zmin

            if config_zmax is None:
                if 'zmaxRegions' in dsRegion:
                    zmax = dsRegion.zmaxRegions.values
                else:
                    # the old naming convention, used in some pre-generated
                    # mask files
                    zmax = dsRegion.zmax.values
            else:
                zmax = config_zmax
            depthMask = numpy.logical_and(zMid >= zmin, zMid <= zmax)
            dsOut = xarray.Dataset()
            dsOut['zmin'] = ('nRegions', [zmin])
            dsOut['zmax'] = ('nRegions', [zmax])
            dsOut['totalArea'] = totalArea
            dsOut['cellMask'] = cellMask
            dsOut['depthMask'] = depthMask
            datasets.append(dsOut)

        dsOut = xarray.concat(objs=datasets, dim='nRegions')
        zbounds = numpy.zeros((nRegions, 2))
        zbounds[:, 0] = dsOut.zmin.values
        zbounds[:, 1] = dsOut.zmax.values
        dsOut['zbounds'] = (('nRegions', 'nbounds'), zbounds)
        dsOut['areaCell'] = areaCell
        dsOut['regionNames'] = dsRegionMask.regionNames
        write_netcdf(dsOut, outFileName)
        # }}}
    # }}}


class ComputeRegionTimeSeriesSubtask(AnalysisTask):  # {{{
    """
    Compute regional and depth mean at a function of time for a set of MPAS
    fields

    Attributes
    ----------
    startYear, endYear : int
        The beginning and end of the time series to compute

    masksSubtask : ``ComputeRegionMasksSubtask``
        A task for creating mask files for each region to plot

    regionGroup : str
        The name of the region group being computed (e.g. "Antarctic Basins")

    regionNames : list of str
        The names of the regions to compute
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, startYear, endYear, masksSubtask,
                 regionGroup, regionNames):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : TimeSeriesOceanRegions
            The main task of which this is a subtask

        startYear, endYear : int
            The beginning and end of the time series to compute

        masksSubtask : ``ComputeRegionMasksSubtask``
            A task for creating mask files for each region to plot

        regionGroup : str
            The name of the region group being computed (e.g. "Antarctic
            Basins")

        regionNames : list of str
            The names of the regions to compute
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        suffix = regionGroup[0].upper() + regionGroup[1:].replace(' ', '')

        # first, call the constructor from the base class (AnalysisTask)
        super(ComputeRegionTimeSeriesSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='compute{}_{:04d}-{:04d}'.format(suffix, startYear,
                                                         endYear))

        parentTask.add_subtask(self)
        self.startYear = startYear
        self.endYear = endYear
        self.masksSubtask = masksSubtask
        self.regionGroup = regionGroup
        self.regionNames = regionNames
        # }}}

    def setup_and_check(self):  # {{{
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
        super(ComputeRegionTimeSeriesSubtask, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        # }}}

    def run_task(self):  # {{{
        """
        Compute the regional-mean time series
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config

        self.logger.info("\nCompute time series of regional means...")

        startDate = '{:04d}-01-01_00:00:00'.format(self.startYear)
        endDate = '{:04d}-12-31_23:59:59'.format(self.endYear)

        regionGroup = self.regionGroup
        sectionSuffix = regionGroup[0].upper() + \
            regionGroup[1:].replace(' ', '')
        timeSeriesName = sectionSuffix
        sectionName = 'timeSeries{}'.format(sectionSuffix)

        outputDirectory = '{}/{}/'.format(
            build_config_full_path(config, 'output', 'timeseriesSubdirectory'),
            timeSeriesName)
        try:
            os.makedirs(outputDirectory)
        except OSError:
            pass

        outFileName = '{}/{}_{:04d}-{:04d}.nc'.format(
            outputDirectory, timeSeriesName, self.startYear, self.endYear)

        inputFiles = sorted(self.historyStreams.readpath(
            'timeSeriesStatsMonthlyOutput', startDate=startDate,
            endDate=endDate, calendar=self.calendar))

        years, months = get_files_year_month(inputFiles,
                                             self.historyStreams,
                                             'timeSeriesStatsMonthlyOutput')

        variables = config.getExpression(sectionName, 'variables')

        variableList = [var['mpas'] for var in variables] + \
            ['timeMonthly_avg_layerThickness']

        outputExists = os.path.exists(outFileName)
        outputValid = outputExists
        if outputExists:
            with open_mpas_dataset(fileName=outFileName,
                                   calendar=self.calendar,
                                   timeVariableNames=None,
                                   variableList=None,
                                   startDate=startDate,
                                   endDate=endDate) as dsOut:

                for inIndex in range(dsOut.dims['Time']):

                    mask = numpy.logical_and(
                        dsOut.year[inIndex].values == years,
                        dsOut.month[inIndex].values == months)
                    if numpy.count_nonzero(mask) == 0:
                        outputValid = False
                        break

        if outputValid:
            self.logger.info('  Time series exists -- Done.')
            return

        regionMaskFileName = '{}/depthMasks_{}.nc'.format(outputDirectory,
                                                          timeSeriesName)
        dsRegionMask = xarray.open_dataset(regionMaskFileName)
        nRegions = dsRegionMask.sizes['nRegions']
        areaCell = dsRegionMask.areaCell

        datasets = []
        nTime = len(inputFiles)
        for tIndex in range(nTime):
            self.logger.info('  {}/{}'.format(tIndex + 1, nTime))

            dsIn = open_mpas_dataset(
                fileName=inputFiles[tIndex],
                calendar=self.calendar,
                variableList=variableList,
                startDate=startDate,
                endDate=endDate).isel(Time=0)
            dsIn.load()

            layerThickness = dsIn.timeMonthly_avg_layerThickness

            innerDatasets = []
            for regionIndex in range(nRegions):
                self.logger.info('    region: {}'.format(
                    self.regionNames[regionIndex]))
                dsRegion = dsRegionMask.isel(nRegions=regionIndex)
                dsRegion.load()
                cellMask = dsRegion.cellMask
                totalArea = dsRegion.totalArea
                depthMask = dsRegion.depthMask.where(cellMask, drop=True)
                localArea = areaCell.where(cellMask, drop=True)
                localThickness = layerThickness.where(cellMask, drop=True)

                volCell = (localArea*localThickness).where(depthMask)
                volCell = volCell.transpose('nCells', 'nVertLevels')
                totalVol = volCell.sum(dim='nVertLevels').sum(dim='nCells')
                self.logger.info('      totalVol (mil. km^3): {}'.format(
                    1e-15*totalVol.values))

                dsOut = xarray.Dataset()
                dsOut['totalVol'] = totalVol
                dsOut.totalVol.attrs['units'] = 'm^3'

                for var in variables:
                    outName = var['name']
                    self.logger.info('      {}'.format(outName))
                    mpasVarName = var['mpas']
                    timeSeries = dsIn[mpasVarName].where(cellMask, drop=True)
                    units = timeSeries.units
                    description = timeSeries.long_name

                    is3d = 'nVertLevels' in timeSeries.dims
                    if is3d:
                        timeSeries = \
                            (volCell*timeSeries.where(depthMask)).sum(
                                dim='nVertLevels').sum(dim='nCells') / totalVol
                    else:
                        timeSeries = \
                            (localArea*timeSeries).sum(
                                dim='nCells') / totalArea

                    dsOut[outName] = timeSeries
                    dsOut[outName].attrs['units'] = units
                    dsOut[outName].attrs['description'] = description
                    dsOut[outName].attrs['is3d'] = str(is3d)

                innerDatasets.append(dsOut)

            datasets.append(innerDatasets)

        # combine data sets into a single data set
        dsOut = xarray.combine_nested(datasets, ['Time', 'nRegions'],
                                      combine_attrs='identical')

        dsOut['totalArea'] = dsRegionMask.totalArea
        dsOut.totalArea.attrs['units'] = 'm^2'
        dsOut['zbounds'] = dsRegionMask.zbounds
        dsOut.zbounds.attrs['units'] = 'm'
        dsOut.coords['regionNames'] = dsRegionMask.regionNames
        dsOut.coords['year'] = (('Time',), years)
        dsOut['year'].attrs['units'] = 'years'
        dsOut.coords['month'] = (('Time',), months)
        dsOut['month'].attrs['units'] = 'months'

        write_netcdf(dsOut, outFileName)  # }}}
    # }}}


class CombineRegionalProfileTimeSeriesSubtask(AnalysisTask):  # {{{
    """
    Combine individual time series into a single data set
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, startYears, endYears, regionGroup):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : TimeSeriesOceanRegions
            The main task of which this is a subtask

        startYears, endYears : list of int
            The beginning and end of each time series to combine

        regionGroup : str
            The name of the region group being computed (e.g. "Antarctic
            Basins")

        """
        # Authors
        # -------
        # Xylar Asay-Davis

        taskSuffix = regionGroup[0].upper() + regionGroup[1:].replace(' ', '')
        subtaskName = 'combine{}TimeSeries'.format(taskSuffix)

        # first, call the constructor from the base class (AnalysisTask)
        super(CombineRegionalProfileTimeSeriesSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName=subtaskName)

        self.startYears = startYears
        self.endYears = endYears
        self.regionGroup = regionGroup
        # }}}

    def run_task(self):  # {{{
        """
        Combine the time series
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        regionGroup = self.regionGroup
        timeSeriesName = regionGroup.replace(' ', '')

        outputDirectory = '{}/{}/'.format(
            build_config_full_path(self.config, 'output',
                                   'timeseriesSubdirectory'),
            timeSeriesName)

        outFileName = '{}/{}_{:04d}-{:04d}.nc'.format(
            outputDirectory, timeSeriesName, self.startYears[0],
            self.endYears[-1])

        if not os.path.exists(outFileName):
            inFileNames = []
            for startYear, endYear in zip(self.startYears, self.endYears):
                inFileName = '{}/{}_{:04d}-{:04d}.nc'.format(
                    outputDirectory, timeSeriesName, startYear, endYear)
                inFileNames.append(inFileName)

            ds = xarray.open_mfdataset(inFileNames, combine='nested',
                                       concat_dim='Time', decode_times=False)

            ds.load()

            # a few variables have become time dependent and shouldn't be
            for var in ['totalArea', 'zbounds']:
                ds[var] = ds[var].isel(Time=0, drop=True)

            write_netcdf(ds, outFileName)
        # }}}
    # }}}


class ComputeObsRegionalTimeSeriesSubtask(AnalysisTask):
    """
    Compute the regional mean of the obs climatology

    Attributes
    ----------
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, regionGroup, regionName, fullSuffix,
                 obsDict):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        regionGroup : str
            Name of the collection of region to plot

        regionName : str
            Name of the region to plot

        fullSuffix : str
            The regionGroup and regionName combined and modified to be
            appropriate as a task or file suffix

        obsDict : dict
            Information on the observations to compare against
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(ComputeObsRegionalTimeSeriesSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='compute{}_{}'.format(fullSuffix, obsDict['suffix']))

        self.regionGroup = regionGroup
        self.regionName = regionName
        self.obsDict = obsDict
        self.prefix = fullSuffix

        timeSeriesName = regionGroup.replace(' ', '')
        outputDirectory = '{}/{}/'.format(
            build_config_full_path(self.config, 'output',
                                   'timeseriesSubdirectory'),
            timeSeriesName)

        self.outFileName = '{}/TS_{}_{}.nc'.format(
            outputDirectory, obsDict['suffix'], self.prefix)

        self.run_after(obsDict['maskTask'])
        # }}}

    def run_task(self):  # {{{
        """
        Compute time-series output of properties in an ocean region.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.logger.info("\nAveraging T and S for {}...".format(
            self.regionName))

        obsDict = self.obsDict
        config = self.config

        regionGroup = self.regionGroup
        timeSeriesName = regionGroup.replace(' ', '')

        sectionSuffix = regionGroup[0].upper() + \
            regionGroup[1:].replace(' ', '')
        sectionName = 'timeSeries{}'.format(sectionSuffix)

        outputDirectory = '{}/{}/'.format(
            build_config_full_path(self.config, 'output',
                                   'timeseriesSubdirectory'),
            timeSeriesName)

        try:
            os.makedirs(outputDirectory)
        except OSError:
            pass

        outFileName = '{}/TS_{}_{}.nc'.format(
            outputDirectory, obsDict['suffix'], self.prefix)

        if os.path.exists(outFileName):
            return

        regionMaskFileName = obsDict['maskTask'].maskFileName

        print(regionMaskFileName)
        print(xarray.open_dataset(regionMaskFileName))

        dsRegionMask = \
            xarray.open_dataset(regionMaskFileName).stack(
                    nCells=(obsDict['latVar'], obsDict['lonVar']))
        dsRegionMask = dsRegionMask.reset_index('nCells').drop_vars(
            [obsDict['latVar'], obsDict['lonVar']])

        maskRegionNames = decode_strings(dsRegionMask.regionNames)
        regionIndex = maskRegionNames.index(self.regionName)

        dsMask = dsRegionMask.isel(nRegions=regionIndex)

        cellMask = dsMask.regionCellMasks == 1

        if config.has_option(sectionName, 'zmin'):
            zmin = config.getfloat(sectionName, 'zmin')
        else:
            zmin = dsMask.zminRegions.values

        if config.has_option(sectionName, 'zmax'):
            zmax = config.getfloat(sectionName, 'zmax')
        else:
            zmax = dsMask.zmaxRegions.values

        TVarName = obsDict['TVar']
        SVarName = obsDict['SVar']
        zVarName = obsDict['zVar']
        lonVarName = obsDict['lonVar']
        latVarName = obsDict['latVar']
        volVarName = obsDict['volVar']
        tDim = obsDict['tDim']

        obsFileName = build_obs_path(
            config, component=self.componentName,
            relativePath=obsDict['TFileName'])
        self.logger.info('  Reading from {}...'.format(obsFileName))

        ds = xarray.open_dataset(obsFileName)
        if obsDict['SFileName'] != obsDict['TFileName']:
            obsFileName = build_obs_path(
                config, component=self.componentName,
                relativePath=obsDict['SFileName'])
            self.logger.info('  Reading from {}...'.format(obsFileName))
            dsS = xarray.open_dataset(obsFileName)
            ds[SVarName] = dsS[SVarName]

        if obsDict['volFileName'] is None:
            # compute volume from lat, lon, depth bounds
            self.logger.info('  Computing volume...'.format(obsFileName))
            latBndsName = ds[latVarName].attrs['bounds']
            lonBndsName = ds[lonVarName].attrs['bounds']
            zBndsName = ds[zVarName].attrs['bounds']
            latBnds = ds[latBndsName]
            lonBnds = ds[lonBndsName]
            zBnds = ds[zBndsName]
            dLat = numpy.deg2rad(latBnds[:, 1] - latBnds[:, 0])
            dLon = numpy.deg2rad(lonBnds[:, 1] - lonBnds[:, 0])
            lat = numpy.deg2rad(ds[latVarName])
            dz = zBnds[:, 1] - zBnds[:, 0]
            radius = 6378137.0
            area = radius**2*numpy.cos(lat)*dLat*dLon
            volume = dz*area
            ds[volVarName] = volume

        elif obsDict['volFileName'] != obsDict['TFileName']:
            obsFileName = build_obs_path(
                config, component=self.componentName,
                relativePath=obsDict['volFileName'])
            self.logger.info('  Reading from {}...'.format(obsFileName))
            dsVol = xarray.open_dataset(obsFileName)
            ds[volVarName] = dsVol[volVarName]

        if 'positive' in ds[zVarName].attrs and \
                ds[zVarName].attrs['positive'] == 'down':
            attrs = ds[zVarName].attrs
            ds[zVarName] = -ds[zVarName]
            ds[zVarName].attrs = attrs
            ds[zVarName].attrs['positive'] = 'up'

        TMean = numpy.zeros(ds.sizes[tDim])
        SMean = numpy.zeros(ds.sizes[tDim])

        depthMask = numpy.logical_and(ds[zVarName] >= zmin,
                                      ds[zVarName] <= zmax)

        for tIndex in range(ds.sizes[tDim]):
            dsMonth = ds.isel({tDim: tIndex})
            dsMonth = dsMonth.stack(nCells=(obsDict['latVar'],
                                            obsDict['lonVar']))
            dsMonth = dsMonth.reset_index('nCells').drop_vars(
                [obsDict['latVar'], obsDict['lonVar']])

            dsMonth = dsMonth.where(cellMask, drop=True)

            dsMonth = dsMonth.where(depthMask)

            mask = dsMonth[TVarName].notnull()
            TSum = (dsMonth[TVarName]*dsMonth[volVarName]).sum(dim=('nCells',
                                                                    zVarName))
            volSum = (mask*dsMonth[volVarName]).sum(dim=('nCells', zVarName))
            TMean[tIndex] = TSum/volSum

            mask = dsMonth[SVarName].notnull()
            SSum = (dsMonth[SVarName]*dsMonth[volVarName]).sum(dim=('nCells',
                                                                    zVarName))
            volSum = (mask*dsMonth[volVarName]).sum(dim=('nCells', zVarName))
            SMean[tIndex] = SSum/volSum

        dsOut = xarray.Dataset()
        dsOut['temperature'] = ('Time', TMean)
        dsOut['salinity'] = ('Time', SMean)
        dsOut['zbounds'] = ('nBounds', [zmin, zmax])
        dsOut['month'] = ('Time', numpy.array(ds.month.values, dtype=float))
        dsOut['year'] = ('Time', numpy.ones(ds.sizes[tDim]))
        write_netcdf(dsOut, outFileName)

        # }}}

    # }}}


class PlotRegionTimeSeriesSubtask(AnalysisTask):
    """
    Plots time-series output of properties in an ocean region.

    Attributes
    ----------
    regionGroup : str
        Name of the collection of region to plot

    regionName : str
        Name of the region to plot

    regionIndex : int
        The index into the dimension ``nRegions`` of the region to plot

    sectionName : str
        The section of the config file to get options from

    controlConfig : ``MpasAnalysisConfigParser``
        The configuration options for the control run (if any)

    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, regionGroup, regionName, regionIndex,
                 controlConfig, sectionName, fullSuffix, obsSubtasks,
                 geojsonFileName):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask :  TimeSeriesOceanRegions
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        regionGroup : str
            Name of the collection of region to plot

        regionName : str
            Name of the region to plot

        regionIndex : int
            The index into the dimension ``nRegions`` of the region to plot

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)

        sectionName : str
            The config section with options for this regionGroup

        fullSuffix : str
            The regionGroup and regionName combined and modified to be
            appropriate as a task or file suffix

        obsSubtasks : dict of ``AnalysisTasks``
            Subtasks for computing the mean observed T and S in the region

        geojsonFileName : str
            The geojson file including the feature to plot
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(PlotRegionTimeSeriesSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='plot{}'.format(fullSuffix))

        self.regionGroup = regionGroup
        self.regionName = regionName
        self.regionIndex = regionIndex
        self.sectionName = sectionName
        self.controlConfig = controlConfig
        self.prefix = fullSuffix
        self.obsSubtasks = obsSubtasks
        self.geojsonFileName = geojsonFileName

        for obsName in obsSubtasks:
            self.run_after(obsSubtasks[obsName])

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        IOError
            If files are not present
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #   self.inDirectory, self.plotsDirectory, self.namelist, self.streams
        #   self.calendar
        super(PlotRegionTimeSeriesSubtask, self).setup_and_check()

        self.variables = self.config.getExpression(self.sectionName,
                                                   'variables')

        self.xmlFileNames = []
        for var in self.variables:
            self.xmlFileNames.append('{}/{}_{}.xml'.format(
                self.plotsDirectory, self.prefix, var['name']))
        return  # }}}

    def run_task(self):  # {{{
        """
        Plots time-series output of properties in an ocean region.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.logger.info("\nPlotting time series of ocean properties of {}"
                         "...".format(self.regionName))

        self.logger.info('  Load time series...')

        config = self.config
        calendar = self.calendar

        fcAll = read_feature_collection(self.geojsonFileName)

        fc = FeatureCollection()
        for feature in fcAll.features:
            if feature['properties']['name'] == self.regionName:
                fc.add_feature(feature)
                break

        baseDirectory = build_config_full_path(
            config, 'output', 'timeSeriesSubdirectory')

        startYear = config.getint('timeSeries', 'startYear')
        endYear = config.getint('timeSeries', 'endYear')
        regionGroup = self.regionGroup
        timeSeriesName = regionGroup.replace(' ', '')

        inFileName = '{}/{}/{}_{:04d}-{:04d}.nc'.format(
            baseDirectory, timeSeriesName, timeSeriesName, startYear, endYear)

        dsIn = xarray.open_dataset(inFileName).isel(nRegions=self.regionIndex)

        zbounds = dsIn.zbounds.values

        controlConfig = self.controlConfig
        plotControl = controlConfig is not None
        if plotControl:
            controlRunName = controlConfig.get('runs', 'mainRunName')
            baseDirectory = build_config_full_path(
                controlConfig, 'output', 'timeSeriesSubdirectory')

            startYear = controlConfig.getint('timeSeries', 'startYear')
            endYear = controlConfig.getint('timeSeries', 'endYear')

            inFileName = '{}/{}/{}_{:04d}-{:04d}.nc'.format(
                baseDirectory, timeSeriesName, timeSeriesName, startYear,
                endYear)
            dsRef = xarray.open_dataset(inFileName).isel(
                nRegions=self.regionIndex)

            zboundsRef = dsRef.zbounds.values

        mainRunName = config.get('runs', 'mainRunName')
        movingAverageMonths = 1

        self.logger.info('  Make plots...')

        groupLink = self.regionGroup.replace(' ', '')

        for var in self.variables:
            varName = var['name']
            mainArray = dsIn[varName]
            is3d = mainArray.attrs['is3d'] == 'True'
            if is3d:
                title = 'Volume-Mean {} in {}'.format(
                    var['title'],  self.regionName)
            else:
                title = 'Area-Mean {} in {}'.format(var['title'],
                                                    self.regionName)

            if plotControl:
                refArray = dsRef[varName]
            xLabel = 'Time (yr)'
            yLabel = '{} ({})'.format(var['title'], var['units'])

            filePrefix = '{}_{}'.format(self.prefix, varName)
            outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

            fields = [mainArray]
            lineColors = ['k']
            lineWidths = [2.5]
            legendText = [mainRunName]
            if plotControl:
                fields.append(refArray)
                lineColors.append('r')
                lineWidths.append(1.2)
                legendText.append(controlRunName)

            if varName in ['temperature', 'salinity']:
                obsColors = ['b', 'g', 'm']
                daysInMonth = constants.daysInMonth
                for obsName in self.obsSubtasks:
                    obsFileName = self.obsSubtasks[obsName].outFileName
                    obsDict = self.obsSubtasks[obsName].obsDict
                    dsObs = xarray.open_dataset(obsFileName)
                    endMonthDays = numpy.cumsum(daysInMonth)
                    midMonthDays = endMonthDays - 0.5*daysInMonth

                    obsTime = []
                    fieldMean = \
                        numpy.sum(dsObs[varName].values*daysInMonth)/365.
                    for year in range(startYear, endYear+1):
                        obsTime.append(midMonthDays + 365.*(year-1.))
                    obsTime = numpy.array(obsTime).ravel()
                    obsField = fieldMean*numpy.ones(obsTime.shape)
                    da = xarray.DataArray(data=obsField, dims='Time',
                                          coords=[('Time', obsTime)])
                    fields.append(da)
                    lineColors.append(obsColors.pop(0))
                    lineWidths.append(1.2)
                    legendText.append(obsDict['legend'])

            if is3d:
                if not plotControl or numpy.all(zbounds == zboundsRef):
                    title = '{} ({} < z < {} m)'.format(title, zbounds[0],
                                                        zbounds[1])
                else:
                    legendText[0] = '{} ({} < z < {} m)'.format(
                        legendText[0], zbounds[0], zbounds[1])
                    legendText[1] = '{} ({} < z < {} m)'.format(
                        legendText[1], zboundsRef[0], zboundsRef[1])

            fig = timeseries_analysis_plot(
                config, fields, calendar=calendar, title=title, xlabel=xLabel,
                ylabel=yLabel, movingAveragePoints=movingAverageMonths,
                lineColors=lineColors, lineWidths=lineWidths,
                legendText=legendText)

            # do this before the inset because otherwise it moves the inset
            # and cartopy doesn't play too well with tight_layout anyway
            plt.tight_layout()

            add_inset(fig, fc, width=2.0, height=2.0)

            savefig(outFileName, tight=False)

            caption = 'Regional mean of {}'.format(title)
            write_image_xml(
                config=config,
                filePrefix=filePrefix,
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup='{} Time Series'.format(self.regionGroup),
                groupLink=groupLink,
                gallery=var['title'],
                thumbnailDescription=self.regionName,
                imageDescription=caption,
                imageCaption=caption)

        # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
