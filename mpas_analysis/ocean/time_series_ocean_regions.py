# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os
import xarray
import numpy
import dask
import multiprocessing
from multiprocessing.pool import ThreadPool

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.plot.plotting import timeseries_analysis_plot

from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf

from mpas_analysis.shared.io.utility import build_config_full_path, \
    get_files_year_month

from mpas_analysis.shared.html import write_image_xml

from mpas_analysis.shared.regions import ComputeRegionMasksSubtask, \
    get_feature_list

from mpas_analysis.ocean.utility import compute_zmid


class TimeSeriesOceanRegions(AnalysisTask):  # {{{
    """
    Performs analysis of the time-series output of regionoal mean temperature,
    salinity, etc.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, controlConfig=None):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

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
            tags=['timeSeries', 'regions'])

        startYear = config.getint('timeSeries', 'startYear')
        endYear = config.getint('timeSeries', 'endYear')

        regionMaskDirectory = build_config_full_path(config,
                                                     'diagnostics',
                                                     'regionMaskSubdirectory')

        regionGroups = config.getExpression(self.taskName, 'regionGroups')

        parallelTaskCount = config.getWithDefault('execute',
                                                  'parallelTaskCount',
                                                  default=1)

        for regionGroup in regionGroups:
            sectionSuffix = regionGroup[0].upper() + \
                regionGroup[1:].replace(' ', '')
            fileSuffix = sectionSuffix[0].lower() + sectionSuffix[1:]
            sectionName = 'timeSeries{}'.format(sectionSuffix)

            regionMaskFile = config.getExpression(sectionName, 'regionMask')

            regionMaskFile = '{}/{}'.format(regionMaskDirectory,
                                            regionMaskFile)
            regionNames = config.getExpression(sectionName, 'regionNames')

            if 'all' in regionNames and os.path.exists(regionMaskFile):
                regionNames = get_feature_list(regionMaskFile)

            masksSubtask = ComputeRegionMasksSubtask(
                self, regionMaskFile, outFileSuffix=fileSuffix,
                featureList=regionNames, subprocessCount=parallelTaskCount)

            self.add_subtask(masksSubtask)

            years = range(startYear, endYear + 1)

            # in the end, we'll combine all the time series into one, but we
            # create this task first so it's easier to tell it to run after all
            # the compute tasks
            combineSubtask = CombineRegionalProfileTimeSeriesSubtask(
                self, startYears=years, endYears=years,
                regionGroup=regionGroup)

            # run one subtask per year
            for year in years:
                computeSubtask = ComputeRegionTimeSeriesSubtask(
                    self, startYear=year, endYear=year,
                    masksSubtask=masksSubtask, regionGroup=regionGroup,
                    regionNames=regionNames)
                self.add_subtask(computeSubtask)
                computeSubtask.run_after(masksSubtask)
                combineSubtask.run_after(computeSubtask)

            self.add_subtask(combineSubtask)

            for index, regionName in enumerate(regionNames):

                fullSuffix = sectionSuffix + '_' + regionName[0].lower() + \
                    regionName[1:].replace(' ', '')

                plotRegionSubtask = PlotRegionTimeSeriesSubtask(
                    self, regionGroup, regionName, index, controlConfig,
                    sectionName, fullSuffix)
                plotRegionSubtask.run_after(combineSubtask)
                self.add_subtask(plotRegionSubtask)

        # }}}

    # }}}


class ComputeRegionTimeSeriesSubtask(AnalysisTask):  # {{{
    '''
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
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, startYear, endYear, masksSubtask,
                 regionGroup, regionNames):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``TimeSeriesOceanRegions``
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
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        suffix = regionGroup[0].upper() + \
            regionGroup[1:].replace(' ', '')

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

        parallelTaskCount = self.config.getint('execute', 'parallelTaskCount')
        self.subprocessCount = min(parallelTaskCount,
                                   self.config.getint(self.taskName,
                                                      'subprocessCount'))
        self.daskThreads = min(
            multiprocessing.cpu_count(),
            self.config.getint(self.taskName, 'daskThreads'))
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
        super(ComputeRegionTimeSeriesSubtask, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        # }}}

    def run_task(self):  # {{{
        '''
        Compute the regional-mean time series
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.logger.info("\nCompute time series of regional means...")

        startDate = '{:04d}-01-01_00:00:00'.format(self.startYear)
        endDate = '{:04d}-12-31_23:59:59'.format(self.endYear)

        regionGroup = self.regionGroup
        sectionSuffix = regionGroup[0].upper() + \
            regionGroup[1:].replace(' ', '')
        timeSeriesName = sectionSuffix[0].lower() + sectionSuffix[1:]
        sectionName = 'timeSeries{}'.format(sectionSuffix)

        outputDirectory = '{}/{}/'.format(
            build_config_full_path(self.config, 'output',
                                   'timeseriesSubdirectory'),
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

        variables = self.config.getExpression(sectionName, 'variables')

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

        # Load mesh related variables
        try:
            restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for Antarctic melt calculations')

        cellsChunk = 32768
        timeChunk = 1

        datasets = []
        for timeIndex, fileName in enumerate(inputFiles):

            dsTimeSlice = open_mpas_dataset(
                fileName=fileName,
                calendar=self.calendar,
                variableList=variableList,
                startDate=startDate,
                endDate=endDate)
            datasets.append(dsTimeSlice)

        chunk = {'Time': timeChunk, 'nCells': cellsChunk}

        with dask.config.set(schedular='threads',
                             pool=ThreadPool(self.daskThreads)):
            # combine data sets into a single data set
            dsIn = xarray.concat(datasets, 'Time').chunk(chunk)

            chunk = {'nCells': cellsChunk}
            dsRestart = xarray.open_dataset(restartFileName)
            dsRestart = dsRestart.isel(Time=0).chunk(chunk)
            dsIn['areaCell'] = dsRestart.areaCell
            if 'landIceMask' in dsRestart:
                # only the region outside of ice-shelf cavities
                dsIn['openOceanMask'] = dsRestart.landIceMask == 0

            dsIn['zMid'] = compute_zmid(dsRestart.bottomDepth,
                                        dsRestart.maxLevelCell,
                                        dsRestart.layerThickness)

            regionMaskFileName = self.masksSubtask.maskFileName

            dsRegionMask = xarray.open_dataset(regionMaskFileName)

            maskRegionNames = [bytes.decode(name) for name in
                               dsRegionMask.regionNames.values]

            datasets = []
            regionIndices = []
            for regionName in self.regionNames:

                self.logger.info('    region: {}'.format(regionName))
                regionIndex = maskRegionNames.index(regionName)
                regionIndices.append(regionIndex)

                chunk = {'nCells': cellsChunk}
                dsMask = dsRegionMask.isel(nRegions=regionIndex).chunk(chunk)

                cellMask = dsMask.regionCellMasks == 1
                if 'openOceanMask' in dsIn:
                    cellMask = numpy.logical_and(cellMask, dsIn.openOceanMask)
                dsRegion = dsIn.where(cellMask, drop=True)

                totalArea = dsRegion['areaCell'].sum()
                self.logger.info('      totalArea: {} mil. km^2'.format(
                    1e-12*totalArea.values))

                self.logger.info("Don't worry about the following dask "
                                 "warnings.")
                depthMask = numpy.logical_and(dsRegion.zMid >= dsMask.zmin,
                                              dsRegion.zMid <= dsMask.zmax)
                depthMask.compute()
                self.logger.info("Dask warnings should be done.")
                dsRegion['depthMask'] = depthMask

                layerThickness = dsRegion.timeMonthly_avg_layerThickness
                dsRegion['volCell'] = (dsRegion.areaCell*layerThickness).where(
                    depthMask)
                totalVol = dsRegion.volCell.sum(dim='nVertLevels').sum(
                    dim='nCells')
                totalVol.compute()
                self.logger.info('      totalVol (mil. km^3): {}'.format(
                    1e-15*totalVol.values))

                dsRegion = dsRegion.transpose('Time', 'nCells', 'nVertLevels')

                dsOut = xarray.Dataset()
                dsOut['totalVol'] = totalVol
                dsOut.totalVol.attrs['units'] = 'm^3'
                dsOut['totalArea'] = totalArea
                dsOut.totalArea.attrs['units'] = 'm^2'

                for var in variables:
                    outName = var['name']
                    self.logger.info('      {}'.format(outName))
                    mpasVarName = var['mpas']
                    timeSeries = dsRegion[mpasVarName]
                    units = timeSeries.units
                    description = timeSeries.long_name

                    if 'nVertLevels' in timeSeries.dims:
                        timeSeries = \
                            (dsRegion.volCell*timeSeries.where(depthMask)).sum(
                                dim='nVertLevels').sum(dim='nCells') / totalVol
                    else:
                        timeSeries = \
                            (dsRegion.areaCell*timeSeries.where(cellMask)).sum(
                                dim='nCells') / totalArea

                    timeSeries.compute()

                    dsOut[outName] = timeSeries
                    dsOut[outName].attrs['units'] = units
                    dsOut[outName].attrs['description'] = description

                datasets.append(dsOut)

            # combine data sets into a single data set
            dsOut = xarray.concat(datasets, 'nRegions')

            dsOut.coords['regionNames'] = dsRegionMask.regionNames.isel(
                nRegions=regionIndices)
            dsOut.coords['year'] = (('Time'), years)
            dsOut['year'].attrs['units'] = 'years'
            dsOut.coords['month'] = (('Time'), months)
            dsOut['month'].attrs['units'] = 'months'

            write_netcdf(dsOut, outFileName)
        # }}}
    # }}}


class CombineRegionalProfileTimeSeriesSubtask(AnalysisTask):  # {{{
    '''
    Combine individual time series into a single data set
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, startYears, endYears, regionGroup):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``TimeSeriesOceanRegions``
            The main task of which this is a subtask

        startYears, endYears : list of int
            The beginning and end of each time series to combine

        regionGroup : str
            The name of the region group being computed (e.g. "Antarctic
            Basins")

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

        self.startYears = startYears
        self.endYears = endYears
        self.regionGroup = regionGroup
        # }}}

    def run_task(self):  # {{{
        '''
        Combine the time series
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        regionGroup = self.regionGroup
        timeSeriesName = regionGroup[0].lower() + \
            regionGroup[1:].replace(' ', '')

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

            ds = xarray.open_mfdataset(inFileNames, concat_dim='Time',
                                       decode_times=False)

            write_netcdf(ds, outFileName)
        # }}}
    # }}}


class PlotRegionTimeSeriesSubtask(AnalysisTask):
    """
    Plots time-series output of Antarctic sub-ice-shelf melt rates.

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
                 controlConfig, sectionName, fullSuffix):
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

        regionIndex : int
            The index into the dimension ``nRegions`` of the region to plot

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)

        sectionName : str
            The config section with options for this regionGroup

        fullSuffix : str
            The regionGroup and regionName combined and modified to be
            appropriate as a task or file suffix
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
        self.prefix = fullSuffix[0].lower() + fullSuffix[1:]

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

        baseDirectory = build_config_full_path(
            config, 'output', 'timeSeriesSubdirectory')

        inFileName = '{}/{}.nc'.format(baseDirectory, self.prefix)

        startYear = config.getint('timeSeries', 'startYear')
        endYear = config.getint('timeSeries', 'endYear')
        regionGroup = self.regionGroup
        timeSeriesName = regionGroup[0].lower() + \
            regionGroup[1:].replace(' ', '')

        inFileName = '{}/{}/{}_{:04d}-{:04d}.nc'.format(
            baseDirectory, timeSeriesName, timeSeriesName, startYear, endYear)

        dsIn = xarray.open_dataset(inFileName).isel(nRegions=self.regionIndex)

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
            dsRef = xarray.open_dataset()

        mainRunName = config.get('runs', 'mainRunName')
        movingAverageMonths = 1

        self.logger.info('  Make plots...')

        groupLink = self.regionGroup[0].lower() + \
            self.regionGroup[1:].replace(' ', '')

        for var in self.variables:
            varName = var['name']
            title = '{} in {}'.format(var['title'],  self.regionName)
            mainArray = dsIn[varName]
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

            timeseries_analysis_plot(config, fields, movingAverageMonths,
                                     title, xLabel, yLabel, outFileName,
                                     calendar=calendar,
                                     lineColors=lineColors,
                                     lineWidths=lineWidths,
                                     legendText=legendText)

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
