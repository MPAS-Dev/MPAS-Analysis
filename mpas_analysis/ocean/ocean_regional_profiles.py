# -*- coding: utf-8 -*-
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
#

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr
import numpy as np
import os
import matplotlib.pyplot as plt

from geometric_features import FeatureCollection, read_feature_collection

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.io.utility import build_config_full_path, \
    get_files_year_month, make_directories, decode_strings, get_region_mask
from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf
from mpas_analysis.shared.timekeeping.utility import days_to_datetime
from mpas_analysis.shared.regions import get_feature_list
from mpas_analysis.shared.climatology import compute_climatology
from mpas_analysis.shared.constants import constants
from mpas_analysis.ocean.plot_hovmoller_subtask import PlotHovmollerSubtask
from mpas_analysis.shared.html import write_image_xml
from mpas_analysis.shared.plot import savefig, add_inset


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

    def __init__(self, config, regionMasksTask, controlConfig=None):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        regionMasksTask : ``ComputeRegionMasks``
            A task for computing region masks

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

        self.startYear = config.getint('climatology', 'startYear')
        self.endYear = config.getint('climatology', 'endYear')

        self.fields = config.getExpression('oceanRegionalProfiles', 'fields')

        self.seasons = config.getExpression('oceanRegionalProfiles', 'seasons')

        self.regionMaskSuffix = config.get('oceanRegionalProfiles',
                                           'regionMaskSuffix')

        self.regionNames = config.getExpression('oceanRegionalProfiles',
                                                'regionNames')

        plotHovmoller = config.getboolean('oceanRegionalProfiles',
                                          'plotHovmoller')

        self.regionMaskSuffix = config.get('oceanRegionalProfiles',
                                           'regionMaskSuffix')

        hovmollerGalleryGroup = config.get('oceanRegionalProfiles',
                                           'hovmollerGalleryGroup')

        masksFile = get_region_mask(config,
                                    '{}.geojson'.format(self.regionMaskSuffix))

        masksSubtask = regionMasksTask.add_mask_subtask(
            masksFile, outFileSuffix=self.regionMaskSuffix)

        if 'all' in self.regionNames:
            self.regionNames = get_feature_list(masksFile)

        self.masksSubtask = masksSubtask

        years = range(self.startYear, self.endYear + 1)

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

        if plotHovmoller:
            for field in self.fields:
                prefix = field['prefix']
                for regionName in self.regionNames:
                    subtaskName = 'plotHovmoller_{}_{}'.format(
                        prefix, regionName.replace(' ', '_'))
                    inFileName = '{}/{}_{:04d}-{:04d}.nc'.format(
                        self.regionMaskSuffix, self.regionMaskSuffix,
                        self.startYear, self.endYear)
                    titleName = field['titleName']
                    caption = 'Time series of {} {} vs ' \
                              'depth'.format(regionName.replace('_', ' '),
                                             titleName)
                    hovmollerSubtask = PlotHovmollerSubtask(
                        parentTask=self,
                        regionName=regionName,
                        inFileName=inFileName,
                        outFileLabel='{}_hovmoller'.format(prefix),
                        fieldNameInTitle=titleName,
                        mpasFieldName='{}_mean'.format(prefix),
                        unitsLabel=field['units'],
                        sectionName='{}OceanRegionalHovmoller'.format(prefix),
                        thumbnailSuffix='',
                        imageCaption=caption,
                        galleryGroup=hovmollerGalleryGroup,
                        groupSubtitle=None,
                        groupLink='ocnreghovs',
                        galleryName=titleName,
                        subtaskName=subtaskName,
                        controlConfig=controlConfig,
                        regionMaskFile=masksFile)
                    hovmollerSubtask.run_after(combineSubtask)
                    self.add_subtask(hovmollerSubtask)

        for field in self.fields:
            prefix = field['prefix']
            for regionName in self.regionNames:
                for season in self.seasons:
                    plotSubtask = PlotRegionalProfileTimeSeriesSubtask(
                        self, season, regionName, field, controlConfig)
                    plotSubtask.run_after(combineSubtask)
                    self.add_subtask(plotSubtask)

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

    startYear, endYear : int
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

        startYear, endYear : int
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
        Compute time series of regional profiles
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
        regionNames = decode_strings(dsRegionMask.regionNames)

        regionIndices = []
        for regionToPlot in self.parentTask.regionNames:
            for index, regionName in enumerate(regionNames):
                if regionToPlot == regionName:
                    regionIndices.append(index)
                    break

        # select only those regions we want to plot
        dsRegionMask = dsRegionMask.isel(nRegions=regionIndices)
        cellMasks = dsRegionMask.regionCellMasks
        regionNamesVar = dsRegionMask.regionNames

        totalArea = (cellMasks * areaCell * vertMask).sum('nCells')

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
                    (cellMasks * areaCell * var).sum('nCells') / totalArea

                meanSquaredName = '{}_meanSquared'.format(prefix)
                dsLocal[meanSquaredName] = \
                    (cellMasks * areaCell * var**2).sum('nCells') / totalArea

            # drop the original variables
            dsLocal = dsLocal.drop_vars(variableList)

            datasets.append(dsLocal)

        # combine data sets into a single data set
        dsOut = xr.concat(datasets, 'Time')

        dsOut.coords['regionNames'] = regionNamesVar
        dsOut['totalArea'] = totalArea
        dsOut.coords['year'] = (('Time'), years)
        dsOut['year'].attrs['units'] = 'years'
        dsOut.coords['month'] = (('Time'), months)
        dsOut['month'].attrs['units'] = 'months'

        # Note: restart file, not a mesh file because we need refBottomDepth,
        # not in a mesh file
        try:
            restartFile = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for plotting time series vs. depth')

        with xr.open_dataset(restartFile) as dsRestart:
            depths = dsRestart.refBottomDepth.values
            z = np.zeros(depths.shape)
            z[0] = -0.5 * depths[0]
            z[1:] = -0.5 * (depths[0:-1] + depths[1:])

        dsOut.coords['z'] = (('nVertLevels'), z)
        dsOut['z'].attrs['units'] = 'meters'

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
        parentTask : ``OceanRegionalProfiles``
            The main task of which this is a subtask

        startYear, endYear : list of int
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

        useExisting = False
        if os.path.exists(outputFileName):
            ds = xr.open_dataset(outputFileName, decode_times=False)
            if ds.sizes['Time'] > 0:
                useExisting = True
            else:
                ds.close()

        if not useExisting:

            inFileNames = []
            for startYear, endYear in zip(self.startYears, self.endYears):
                inFileName = '{}/{}_{:04d}-{:04d}.nc'.format(
                    outputDirectory, timeSeriesName, startYear, endYear)
                inFileNames.append(inFileName)

            ds = xr.open_mfdataset(inFileNames, combine='nested',
                                   concat_dim='Time', decode_times=False)

            ds.load()

            ds['totalArea'] = ds['totalArea'].isel(Time=0)

            write_netcdf(ds, outputFileName)

        regionNames = ds['regionNames']
        ds = ds.drop('regionNames')

        profileMask = ds['totalArea'] > 0

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
                write_netcdf(dsSeason, outputFileName)

        # }}}
    # }}}


class PlotRegionalProfileTimeSeriesSubtask(AnalysisTask):  # {{{
    '''
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

    controlConfig :  ``MpasAnalysisConfigParser``
        Configuration options for a control run (if any)
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, season, regionName, field, controlConfig):
        # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``AnalysisTask``
            The parent task of which this is a subtask

        season : str
            The season being plotted

        regionName : str
            The region being plotted

        field : dict
            Information about the field (e.g. temperature) being plotted

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        '''
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

        self.parentTask = parentTask
        self.controlConfig = controlConfig

        self.season = season
        self.regionName = regionName
        self.field = field

        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(PlotRegionalProfileTimeSeriesSubtask, self).setup_and_check()

        self.xmlFileNames = []
        self.filePrefixes = {}

        self.filePrefix = '{}_{}_{}_years{:04d}-{:04d}'.format(
            self.field['prefix'], self.regionName.replace(' ', '_'),
            self.season, self.parentTask.startYear,
            self.parentTask.endYear)
        self.xmlFileNames = ['{}/{}.xml'.format(self.plotsDirectory,
                                                self.filePrefix)]
        # }}}

    def run_task(self):  # {{{
        """
        Plot a depth profile with variability
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config
        startYear = self.parentTask.startYear
        endYear = self.parentTask.endYear

        regionMaskSuffix = config.get('oceanRegionalProfiles',
                                      'regionMaskSuffix')
        regionMaskFile = get_region_mask(config,
                                         '{}.geojson'.format(regionMaskSuffix))

        fcAll = read_feature_collection(regionMaskFile)

        fc = FeatureCollection()
        for feature in fcAll.features:
            if feature['properties']['name'] == self.regionName:
                fc.add_feature(feature)
                break

        inDirectory = build_config_full_path(config, 'output',
                                             'profilesSubdirectory')
        timeSeriesName = self.parentTask.regionMaskSuffix
        inFileName = '{}/{}_{}_{:04d}-{:04d}.nc'.format(
            inDirectory, timeSeriesName, self.season,
            self.parentTask.startYear, self.parentTask.endYear)

        ds = xr.open_dataset(inFileName)
        allRegionNames = decode_strings(ds.regionNames)

        regionIndex = allRegionNames.index(self.regionName)
        ds = ds.isel(nRegions=regionIndex)
        meanFieldName = '{}_mean'.format(self.field['prefix'])
        stdFieldName = '{}_std'.format(self.field['prefix'])

        mainRunName = config.get('runs', 'mainRunName')
        profileGalleryGroup = config.get('oceanRegionalProfiles',
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

            controlFileName = '{}/{}_{}_{:04d}-{:04d}.nc'.format(
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

        depthRange = config.getExpression('oceanRegionalProfiles',
                                          'depthRange')
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

        savefig(outFileName, tight=False)

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
        # }}}

    def plot(self, zArrays, fieldArrays, errArrays, lineColors, lineWidths,
             legendText, title, xLabel, yLabel, xLim=None, yLim=None,
             figureSize=(10, 4), dpi=None):  # {{{
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
            plt.legend()

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
        return fig  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
