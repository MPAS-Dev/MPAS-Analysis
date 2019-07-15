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
import csv
import dask
import multiprocessing
from multiprocessing.pool import ThreadPool
import matplotlib.pyplot as plt

from geometric_features import FeatureCollection, read_feature_collection

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.plot import timeseries_analysis_plot, savefig, \
    add_inset

from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories, build_obs_path, decode_strings, get_region_mask

from mpas_analysis.shared.html import write_image_xml

from mpas_analysis.shared.regions import ComputeRegionMasksSubtask, \
    get_feature_list


class TimeSeriesAntarcticMelt(AnalysisTask):  # {{{
    """
    Performs analysis of the time-series output of Antarctic sub-ice-shelf
    melt rates.
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Stephen Price

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
        super(TimeSeriesAntarcticMelt, self).__init__(
            config=config,
            taskName='timeSeriesAntarcticMelt',
            componentName='ocean',
            tags=['timeSeries', 'melt', 'landIceCavities', 'antarctic'])

        self.iceShelfMasksFile = get_region_mask(config, 'iceShelves.geojson')

        iceShelvesToPlot = config.getExpression('timeSeriesAntarcticMelt',
                                                'iceShelvesToPlot')
        if 'all' in iceShelvesToPlot:
            iceShelvesToPlot = get_feature_list(self.iceShelfMasksFile)

        parallelTaskCount = config.getWithDefault('execute',
                                                  'parallelTaskCount',
                                                  default=1)

        masksSubtask = ComputeRegionMasksSubtask(
            self, self.iceShelfMasksFile, outFileSuffix='iceShelfMasks',
            featureList=iceShelvesToPlot, subprocessCount=parallelTaskCount)

        self.add_subtask(masksSubtask)

        computeMeltSubtask = ComputeMeltSubtask(self, mpasTimeSeriesTask,
                                                masksSubtask, iceShelvesToPlot)
        self.add_subtask(computeMeltSubtask)

        for index, iceShelf in enumerate(iceShelvesToPlot):
            plotMeltSubtask = PlotMeltSubtask(self, iceShelf, index,
                                              controlConfig)
            plotMeltSubtask.run_after(computeMeltSubtask)
            self.add_subtask(plotMeltSubtask)

        # }}}

    # }}}


class ComputeMeltSubtask(AnalysisTask):  # {{{
    """
    Computes time-series of Antarctic sub-ice-shelf melt rates.

    Attributes
    ----------
    mpasTimeSeriesTask : ``MpasTimeSeriesTask``
        The task that extracts the time series from MPAS monthly output

    masksSubtask : ``ComputeRegionMasksSubtask``
        A task for creating mask files for each ice shelf to plot

    iceShelvesToPlot : list of str
        A list of ice shelves to plot
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Stephen Price

    def __init__(self, parentTask, mpasTimeSeriesTask, masksSubtask,
                 iceShelvesToPlot):  # {{{
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
            A task for creating mask files for each ice shelf to plot

        iceShelvesToPlot : list of str
            A list of ice shelves to plot
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(ComputeMeltSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='computeMeltRates')

        self.mpasTimeSeriesTask = mpasTimeSeriesTask
        self.run_after(mpasTimeSeriesTask)

        self.masksSubtask = masksSubtask
        self.run_after(masksSubtask)

        self.iceShelvesToPlot = iceShelvesToPlot

        parallelTaskCount = self.config.getint('execute', 'parallelTaskCount')
        self.subprocessCount = min(parallelTaskCount,
                                   self.config.getint(self.taskName,
                                                      'subprocessCount'))
        self.daskThreads = min(
            multiprocessing.cpu_count(),
            self.config.getint(self.taskName, 'daskThreads'))
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
        super(ComputeMeltSubtask, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        config = self.config

        landIceFluxMode = self.namelist.get('config_land_ice_flux_mode')
        if landIceFluxMode not in ['standalone', 'coupled']:
            raise ValueError('*** timeSeriesAntarcticMelt requires '
                             'config_land_ice_flux_mode \n'
                             '    to be standalone or coupled.  Otherwise, no '
                             'melt rates are available \n'
                             '    for plotting.')

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

        self.variableList = \
            ['timeMonthly_avg_landIceFreshwaterFlux']
        self.mpasTimeSeriesTask.add_variables(variableList=self.variableList)

        return  # }}}

    def run_task(self):  # {{{
        """
        Computes time-series of Antarctic sub-ice-shelf melt rates.
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Stephen Price

        self.logger.info(r"\Computing Antarctic melt rate time series...")

        self.logger.info('  Load melt rate data...')

        mpasTimeSeriesTask = self.mpasTimeSeriesTask
        config = self.config

        baseDirectory = build_config_full_path(
            config, 'output', 'timeSeriesSubdirectory')

        outFileName = '{}/iceShelfAggregatedFluxes.nc'.format(baseDirectory)

        # Load data:
        inputFile = mpasTimeSeriesTask.outputFile
        dsIn = open_mpas_dataset(fileName=inputFile,
                                 calendar=self.calendar,
                                 variableList=self.variableList,
                                 startDate=self.startDate,
                                 endDate=self.endDate)
        try:
            if os.path.exists(outFileName):
                # The file already exists so load it
                dsOut = xarray.open_dataset(outFileName)
                if numpy.all(dsOut.Time.values == dsIn.Time.values):
                    return
                else:
                    self.logger.warning('File {} is incomplete. Deleting '
                                        'it.'.format(outFileName))
                    os.remove(outFileName)
        except OSError:
            # something is potentailly wrong with the file, so let's delete
            # it and try again
            self.logger.warning('Problems reading file {}. Deleting '
                                'it.'.format(outFileName))
            os.remove(outFileName)

        with dask.config.set(schedular='threads',
                             pool=ThreadPool(self.daskThreads)):
            # work on data from simulations
            freshwaterFlux = dsIn.timeMonthly_avg_landIceFreshwaterFlux.chunk(
                {'Time': 12})

            restartFileName = \
                mpasTimeSeriesTask.runStreams.readpath('restart')[0]

            dsRestart = xarray.open_dataset(restartFileName)
            areaCell = \
                dsRestart.landIceFraction.isel(Time=0) * dsRestart.areaCell

            regionMaskFileName = self.masksSubtask.maskFileName

            dsRegionMask = xarray.open_dataset(regionMaskFileName)

            # figure out the indices of the regions to plot
            regionNames = decode_strings(dsRegionMask.regionNames)

            regionIndices = []
            for iceShelf in self.iceShelvesToPlot:
                for index, regionName in enumerate(regionNames):
                    if iceShelf == regionName:
                        regionIndices.append(index)
                        break

            # select only those regions we want to plot
            dsRegionMask = dsRegionMask.isel(nRegions=regionIndices)
            cellMasks = dsRegionMask.regionCellMasks.chunk({'nRegions': 10})

            # convert from kg/s to kg/yr
            totalMeltFlux = constants.sec_per_year * \
                (cellMasks * areaCell * freshwaterFlux).sum(dim='nCells')
            totalMeltFlux.compute()

            totalArea = (cellMasks * areaCell).sum(dim='nCells')

            # from kg/m^2/yr to m/yr
            meltRates = (1. / constants.rho_fw) * (totalMeltFlux / totalArea)
            meltRates.compute()

            # convert from kg/yr to GT/yr
            totalMeltFlux /= constants.kg_per_GT

            dsOut = xarray.Dataset()
            dsOut['totalMeltFlux'] = totalMeltFlux
            dsOut.totalMeltFlux.attrs['units'] = 'GT a$^{-1}$'
            dsOut.totalMeltFlux.attrs['description'] = \
                'Total melt flux summed over each ice shelf or region'
            dsOut['meltRates'] = meltRates
            dsOut.meltRates.attrs['units'] = 'm a$^{-1}$'
            dsOut.meltRates.attrs['description'] = \
                'Melt rate averaged over each ice shelf or region'

            write_netcdf(dsOut, outFileName)

        # }}}

    # }}}


class PlotMeltSubtask(AnalysisTask):
    """
    Plots time-series output of Antarctic sub-ice-shelf melt rates.

    Attributes
    ----------
    iceShelf : str
        Name of the ice shelf to plot

    regionIndex : int
        The index into the dimension ``nRegions`` of the ice shelf to plot

    controlConfig : ``MpasAnalysisConfigParser``
        The configuration options for the control run (if any)

    """
    # Authors
    # -------
    # Xylar Asay-Davis, Stephen Price

    def __init__(self, parentTask, iceShelf, regionIndex, controlConfig):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        iceShelf : str
            Name of the ice shelf to plot

        regionIndex : int
            The index into the dimension ``nRegions`` of the ice shelf to plot

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(PlotMeltSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='plotMeltRates_{}'.format(iceShelf.replace(' ', '_')))

        self.iceShelfMasksFile = parentTask.iceShelfMasksFile
        self.iceShelf = iceShelf
        self.regionIndex = regionIndex
        self.controlConfig = controlConfig

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
        super(PlotMeltSubtask, self).setup_and_check()

        self.xmlFileNames = []

        for prefix in ['melt_flux', 'melt_rate']:
            self.xmlFileNames.append(
                '{}/{}_{}.xml'.format(self.plotsDirectory, prefix,
                                      self.iceShelf.replace(' ', '_')))
        return  # }}}

    def run_task(self):  # {{{
        """
        Plots time-series output of Antarctic sub-ice-shelf melt rates.
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Stephen Price

        self.logger.info("\nPlotting Antarctic melt rate time series for "
                         "{}...".format(self.iceShelf))

        self.logger.info('  Load melt rate data...')

        config = self.config
        calendar = self.calendar

        iceShelfMasksFile = self.iceShelfMasksFile

        fcAll = read_feature_collection(iceShelfMasksFile)

        fc = FeatureCollection()
        for feature in fcAll.features:
            if feature['properties']['name'] == self.iceShelf:
                fc.add_feature(feature)
                break

        totalMeltFlux, meltRates = self._load_ice_shelf_fluxes(config)

        plotControl = self.controlConfig is not None
        if plotControl:
            controlRunName = self.controlConfig.get('runs', 'mainRunName')

            refTotalMeltFlux, refMeltRates = \
                self._load_ice_shelf_fluxes(self.controlConfig)

        # Load observations from multiple files and put in dictionary based
        # on shelf keyname
        observationsDirectory = build_obs_path(config, 'ocean',
                                               'meltSubdirectory')
        obsFileNameDict = {'Rignot et al. (2013)':
                           'Rignot_2013_melt_rates.csv',
                           'Rignot et al. (2013) SS':
                           'Rignot_2013_melt_rates_SS.csv'}

        obsDict = {}  # dict for storing dict of obs data
        for obsName in obsFileNameDict:
            obsFileName = '{}/{}'.format(observationsDirectory,
                                         obsFileNameDict[obsName])
            obsDict[obsName] = {}
            obsFile = csv.reader(open(obsFileName, 'rU'))
            next(obsFile, None)  # skip the header line
            for line in obsFile:  # some later useful values commented out
                shelfName = line[0]
                if shelfName != self.iceShelf:
                    continue

                # surveyArea = line[1]
                meltFlux = float(line[2])
                meltFluxUncertainty = float(line[3])
                meltRate = float(line[4])
                meltRateUncertainty = float(line[5])
                # actualArea = float( line[6] )  # actual area here is in sq km

                # build dict of obs. keyed to filename description
                # (which will be used for plotting)
                obsDict[obsName] = {
                    'meltFlux': meltFlux,
                    'meltFluxUncertainty': meltFluxUncertainty,
                    'meltRate': meltRate,
                    'meltRateUncertainty': meltRateUncertainty}
                break

        # If areas from obs file used need to be converted from sq km to sq m

        mainRunName = config.get('runs', 'mainRunName')
        movingAverageMonths = config.getint('timeSeriesAntarcticMelt',
                                            'movingAverageMonths')

        outputDirectory = build_config_full_path(config, 'output',
                                                 'timeseriesSubdirectory')

        make_directories(outputDirectory)

        self.logger.info('  Make plots...')

        # get obs melt flux and unc. for shelf (similar for rates)
        obsMeltFlux = []
        obsMeltFluxUnc = []
        obsMeltRate = []
        obsMeltRateUnc = []
        for obsName in obsDict:
            if len(obsDict[obsName]) > 0:
                obsMeltFlux.append(
                    obsDict[obsName]['meltFlux'])
                obsMeltFluxUnc.append(
                    obsDict[obsName]['meltFluxUncertainty'])
                obsMeltRate.append(
                    obsDict[obsName]['meltRate'])
                obsMeltRateUnc.append(
                    obsDict[obsName]['meltRateUncertainty'])
            else:
                # append NaN so this particular obs won't plot
                self.logger.warning('{} observations not available for '
                                    '{}'.format(obsName, self.iceShelf))
                obsMeltFlux.append(None)
                obsMeltFluxUnc.append(None)
                obsMeltRate.append(None)
                obsMeltRateUnc.append(None)

        title = self.iceShelf.replace('_', ' ')

        xLabel = 'Time (yr)'
        yLabel = 'Melt Flux (GT/yr)'

        timeSeries = totalMeltFlux.isel(nRegions=self.regionIndex)

        filePrefix = 'melt_flux_{}'.format(self.iceShelf.replace(' ', '_'))
        outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

        fields = [timeSeries]
        lineColors = ['k']
        lineWidths = [2.5]
        legendText = [mainRunName]
        if plotControl:
            fields.append(refTotalMeltFlux.isel(nRegions=self.regionIndex))
            lineColors.append('r')
            lineWidths.append(1.2)
            legendText.append(controlRunName)

        fig = timeseries_analysis_plot(config, fields, movingAverageMonths,
                                       title, xLabel, yLabel,
                                       calendar=calendar,
                                       lineColors=lineColors,
                                       lineWidths=lineWidths,
                                       legendText=legendText,
                                       obsMean=obsMeltFlux,
                                       obsUncertainty=obsMeltFluxUnc,
                                       obsLegend=list(obsDict.keys()))

        # do this before the inset because otherwise it moves the inset
        # and cartopy doesn't play too well with tight_layout anyway
        plt.tight_layout()

        add_inset(fig, fc, width=2.0, height=2.0)

        savefig(outFileName)

        caption = 'Running Mean of Total Melt Flux  under Ice ' \
                  'Shelves in the {} Region'.format(title)
        write_image_xml(
            config=config,
            filePrefix=filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='Antarctic Melt Time Series',
            groupLink='antmelttime',
            gallery='Total Melt Flux',
            thumbnailDescription=title,
            imageDescription=caption,
            imageCaption=caption)

        xLabel = 'Time (yr)'
        yLabel = 'Melt Rate (m/yr)'

        timeSeries = meltRates.isel(nRegions=self.regionIndex)

        filePrefix = 'melt_rate_{}'.format(self.iceShelf.replace(' ', '_'))
        outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

        fields = [timeSeries]
        lineColors = ['k']
        lineWidths = [2.5]
        legendText = [mainRunName]
        if plotControl:
            fields.append(refMeltRates.isel(nRegions=self.regionIndex))
            lineColors.append('r')
            lineWidths.append(1.2)
            legendText.append(controlRunName)

        if config.has_option(self.taskName, 'firstYearXTicks'):
            firstYearXTicks = config.getint(self.taskName,
                                            'firstYearXTicks')
        else:
            firstYearXTicks = None

        if config.has_option(self.taskName, 'yearStrideXTicks'):
            yearStrideXTicks = config.getint(self.taskName,
                                             'yearStrideXTicks')
        else:
            yearStrideXTicks = None

        fig = timeseries_analysis_plot(config, fields, movingAverageMonths,
                                       title, xLabel, yLabel,
                                       calendar=calendar,
                                       lineColors=lineColors,
                                       lineWidths=lineWidths,
                                       legendText=legendText,
                                       obsMean=obsMeltRate,
                                       obsUncertainty=obsMeltRateUnc,
                                       obsLegend=list(obsDict.keys()),
                                       firstYearXTicks=firstYearXTicks,
                                       yearStrideXTicks=yearStrideXTicks)

        # do this before the inset because otherwise it moves the inset
        # and cartopy doesn't play too well with tight_layout anyway
        plt.tight_layout()

        add_inset(fig, fc, width=2.0, height=2.0)

        savefig(outFileName)

        caption = 'Running Mean of Area-averaged Melt Rate under Ice ' \
                  'Shelves in the {} Region'.format(title)
        write_image_xml(
            config=config,
            filePrefix=filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='Antarctic Melt Time Series',
            groupLink='antmelttime',
            gallery='Area-averaged Melt Rate',
            thumbnailDescription=title,
            imageDescription=caption,
            imageCaption=caption)
        # }}}

    def _load_ice_shelf_fluxes(self, config):  # {{{
        """
        Reads melt flux time series and computes regional total melt flux and
        mean melt rate.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        baseDirectory = build_config_full_path(
            config, 'output', 'timeSeriesSubdirectory')

        outFileName = '{}/iceShelfAggregatedFluxes.nc'.format(baseDirectory)

        dsOut = xarray.open_dataset(outFileName)
        return dsOut.totalMeltFlux, dsOut.meltRates
        # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
