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

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.plot.plotting import timeseries_analysis_plot

from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf

from mpas_analysis.shared.io.utility import build_config_full_path

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

    def __init__(self, config, mpasTimeSeriesTask, controlConfig=None):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

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

        regionMaskDirectory = build_config_full_path(config,
                                                     'diagnostics',
                                                     'regionMaskSubdirectory')

        regionGroups = config.getExpression(self.taskName, 'regionGroups')

        for regionGroup in regionGroups:
            sectionSuffix = regionGroup[0].upper() + \
                regionGroup[1:].replace(' ', '')
            fileSuffix = sectionSuffix[0].lower() + sectionSuffix[1:]
            sectionName = 'timeSeries{}'.format(sectionSuffix)

            regionMaskFile = config.getExpression(sectionName, 'regionMask')

            regionMaskFile = '{}/{}'.format(regionMaskDirectory,
                                            regionMaskFile)
            regionNames = config.getExpression(sectionName, 'regionNames')

            if 'all' in regionNames:
                regionNames = get_feature_list(regionMaskFile)

            masksSubtask = ComputeRegionMasksSubtask(
                self, regionMaskFile, outFileSuffix=fileSuffix,
                featureList=regionNames)

            self.add_subtask(masksSubtask)

            for index, regionName in enumerate(regionNames):

                fullSuffix = sectionSuffix + '_' + regionName[0].lower() + \
                    regionName[1:].replace(' ', '')

                computeTimeSeriesTask = ComputeRegionTimeSeriesSubtask(
                        self, mpasTimeSeriesTask, masksSubtask, regionGroup,
                        regionName, sectionName, fullSuffix)
                self.add_subtask(computeTimeSeriesTask)
                plotRegionSubtask = PlotRegionTimeSeriesSubtask(
                    self, regionGroup, regionName, index,  controlConfig,
                    sectionName, fullSuffix)
                plotRegionSubtask.run_after(computeTimeSeriesTask)
                self.add_subtask(plotRegionSubtask)

        # }}}

    # }}}


class ComputeRegionTimeSeriesSubtask(AnalysisTask):  # {{{
    """
    Computes time-series of temperature, salinity, etc. in an ocean region

    Attributes
    ----------
    mpasTimeSeriesTask : ``MpasTimeSeriesTask``
        The task that extracts the time series from MPAS monthly output

    masksSubtask : ``ComputeRegionMasksSubtask``
        A task for creating mask files for each region to plot

    regionGroup : str
        The name of the region group being computed (e.g. "Antarctic Basins")

    regionName : str
        A region to compute time series for
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, mpasTimeSeriesTask, masksSubtask,
                 regionGroup, regionName, sectionName, fullSuffix):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

        masksSubtask : ``ComputeRegionMasksSubtask``
            A task for creating mask files for each region to plot

        regionGroup : str
            The name of the region group being computed (e.g. "Antarctic
            Basins")

        regionName : str
            A region to compute time series for

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
        super(ComputeRegionTimeSeriesSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='compute{}'.format(fullSuffix))

        self.mpasTimeSeriesTask = mpasTimeSeriesTask
        self.run_after(mpasTimeSeriesTask)

        self.masksSubtask = masksSubtask
        self.run_after(masksSubtask)

        self.regionGroup = regionGroup

        self.regionName = regionName

        baseDirectory = build_config_full_path(
            self.config, 'output', 'timeSeriesSubdirectory')

        self.outFileName = '{}/{}.nc'.format(
            baseDirectory, fullSuffix[0].lower() + fullSuffix[1:])
        self.sectionName = sectionName

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        IOError
            If a restart file is not present

        ValueError
            If ``config_land_ice_flux_mode`` is not one of ``standalone`` or
            ``coupled``
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #   self.inDirectory, self.plotsDirectory, self.namelist, self.streams
        #   self.calendar
        super(ComputeRegionTimeSeriesSubtask, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        config = self.config

        # Load mesh related variables
        try:
            self.restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for Antarctic melt calculations')

        # get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        self.startDate = config.get('timeSeries', 'startDate')
        self.endDate = config.get('timeSeries', 'endDate')

        self.variables = config.getExpression(self.sectionName, 'variables')

        self.variableList = [var['mpas'] for var in self.variables] + \
            ['timeMonthly_avg_layerThickness']

        self.mpasTimeSeriesTask.add_variables(variableList=self.variableList)

        return  # }}}

    def run_task(self):  # {{{
        """
        Computes time-series of Antarctic sub-ice-shelf melt rates.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.logger.info(r"\Computing {} regional-mean time series...")

        self.logger.info('  Load melt rate data...')

        mpasTimeSeriesTask = self.mpasTimeSeriesTask

        # Load data:
        inputFile = mpasTimeSeriesTask.outputFile
        dsIn = open_mpas_dataset(fileName=inputFile,
                                 calendar=self.calendar,
                                 variableList=self.variableList,
                                 startDate=self.startDate,
                                 endDate=self.endDate)

        try:
            if os.path.exists(self.outFileName):
                # The file already exists so load it
                dsOut = xarray.open_dataset(self.outFileName)
                if numpy.all(dsOut.Time.values == dsIn.Time.values):
                    return
                else:
                    self.logger.warning('File {} is incomplete. Deleting '
                                        'it.'.format(self.outFileName))
                    os.remove(self.outFileName)
        except OSError:
            # something is potentailly wrong with the file, so let's delete
            # it and try again
            self.logger.warning('Problems reading file {}. Deleting '
                                'it.'.format(self.outFileName))
            os.remove(self.outFileName)

        anyVarsWithDepth = False
        for var in self.variables:
            mpasVarName = var['mpas']
            if 'nVertLevels' in dsIn[mpasVarName].dims:
                anyVarsWithDepth = True

        with xarray.open_dataset(self.restartFileName) as dsRestart:
            dsRestart = dsRestart.isel(Time=0)
            areaCell = dsRestart.areaCell
            if 'landIceMask' in dsRestart:
                # only the region outside of ice-shelf cavities
                openOceanMask = dsRestart.landIceMask == 0
            else:
                openOceanMask = None

            if anyVarsWithDepth:
                zMid = compute_zmid(dsRestart.bottomDepth,
                                    dsRestart.maxLevelCell,
                                    dsRestart.layerThickness)

            dsRestart.close()

        regionMaskFileName = self.masksSubtask.maskFileName

        dsRegionMask = xarray.open_dataset(regionMaskFileName)

        # figure out the indices of the regions to plot
        maskRegionNames = [bytes.decode(name) for name in
                           dsRegionMask.regionNames.values]

        regionIndex = maskRegionNames.index(self.regionName)

        dsRegionMask = dsRegionMask.isel(nRegions=regionIndex)
        cellMask = dsRegionMask.regionCellMasks == 1
        if openOceanMask is not None:
            cellMask = numpy.logical_and(cellMask, openOceanMask)

        areaCell = areaCell.where(cellMask)
        totalArea = areaCell.sum()

        if anyVarsWithDepth:
            zMin = dsRegionMask.zmin.values
            zMax = dsRegionMask.zmax.values

            depthMask = numpy.logical_and(zMid >= zMin, zMid <= zMax)
            depthMask = numpy.logical_and(depthMask, cellMask)

        dsOut = xarray.Dataset()

        for var in self.variables:
            mpasVarName = var['mpas']
            hasDepth = 'nVertLevels' in dsIn[mpasVarName].dims
            if hasDepth:
                chunk = {'Time': 1}
            else:
                chunk = {'Time': 12}
            timeSeries = dsIn[mpasVarName].chunk(chunk)
            print(dsIn[mpasVarName])
            units = timeSeries.units
            description = timeSeries.long_name

            if hasDepth:
                layerThickness = dsIn.timeMonthly_avg_layerThickness.where(
                    depthMask)
                volCell = areaCell*layerThickness
                timeSeries = (volCell*timeSeries.where(depthMask)).sum(
                    dim='nVertLevels').sum(dim='nCells')
                totalVol = volCell.sum(dim='nVertLevels').sum(dim='nCells')
                timeSeries = timeSeries/totalVol
            else:
                timeSeries = (areaCell*timeSeries.where(cellMask)).sum(
                    dim='nCells') / totalArea

            timeSeries.compute()

            outName = var['name']
            dsOut[outName] = timeSeries
            dsOut[outName].attrs['units'] = units
            dsOut[outName].attrs['description'] = description

        write_netcdf(dsOut, self.outFileName)

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

        dsIn = xarray.open_dataset(inFileName)

        plotControl = self.controlConfig is not None
        if plotControl:
            controlRunName = self.controlConfig.get('runs', 'mainRunName')
            baseDirectory = build_config_full_path(
                self.controlConfig, 'output', 'timeSeriesSubdirectory')

            inFileName = '{}/{}.nc'.format(baseDirectory, self.prefix)
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
