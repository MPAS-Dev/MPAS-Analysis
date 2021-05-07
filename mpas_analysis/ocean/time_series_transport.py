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

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.plot import timeseries_analysis_plot, savefig, \
    add_inset

from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf

from mpas_analysis.shared.io.utility import build_config_full_path, \
    get_files_year_month, decode_strings

from mpas_analysis.shared.html import write_image_xml

from mpas_analysis.shared.transects import ComputeTransectMasksSubtask


class TimeSeriesTransport(AnalysisTask):  # {{{
    """
    Extract and plot time series of transport through transects on the MPAS
    mesh.
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Stephen Price

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
        super(TimeSeriesTransport, self).__init__(
            config=config,
            taskName='timeSeriesTransport',
            componentName='ocean',
            tags=['timeSeries', 'transport'])

        startYear = config.getint('timeSeries', 'startYear')
        endYear = config.get('timeSeries', 'endYear')
        if endYear == 'end':
            # a valid end year wasn't found, so likely the run was not found,
            # perhaps because we're just listing analysis tasks
            endYear = startYear
        else:
            endYear = int(endYear)

        years = [year for year in range(startYear, endYear + 1)]

        transectsToPlot = config.getExpression('timeSeriesTransport',
                                               'transectsToPlot')
        if len(transectsToPlot) == 0:
            return

        masksSubtask = ComputeTransectMasksSubtask(
            parentTask=self, transectGroup='Transport Transects')

        transectsToPlot = masksSubtask.expand_transect_names(transectsToPlot)
        transportTransectFileName = masksSubtask.geojsonFileName

        self.add_subtask(masksSubtask)

        # in the end, we'll combine all the time series into one, but we
        # create this task first so it's easier to tell it to run after all
        # the compute tasks
        combineSubtask = CombineTransportSubtask(
            self, startYears=years, endYears=years)

        # run one subtask per year
        for year in years:
            computeSubtask = ComputeTransportSubtask(
                self, startYear=year, endYear=year, masksSubtask=masksSubtask,
                transectsToPlot=transectsToPlot)
            self.add_subtask(computeSubtask)
            computeSubtask.run_after(masksSubtask)
            combineSubtask.run_after(computeSubtask)

        for index, transect in enumerate(transectsToPlot):
            plotTransportSubtask = PlotTransportSubtask(
                self, transect, index, controlConfig, transportTransectFileName)
            plotTransportSubtask.run_after(combineSubtask)
            self.add_subtask(plotTransportSubtask)

        # }}}

    # }}}


class ComputeTransportSubtask(AnalysisTask):  # {{{
    """
    Computes time-series of transport through transects.

    Attributes
    ----------
    startYear, endYear : int
        The beginning and end of the time series to compute

    masksSubtask : ``ComputeRegionMasksSubtask``
        A task for creating mask files for each ice shelf to plot

    transectsToPlot : list of str
        A list of transects to plot
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Stephen Price

    def __init__(self, parentTask, startYear, endYear,
                 masksSubtask, transectsToPlot):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        startYear, endYear : int
            The beginning and end of the time series to compute

        masksSubtask : ``ComputeRegionMasksSubtask``
            A task for creating mask files for each ice shelf to plot

        transectsToPlot : list of str
            A list of transects to plot
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(ComputeTransportSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='computeTransport_{:04d}-{:04d}'.format(startYear,
                                                                endYear))

        self.subprocessCount = self.config.getint('timeSeriesTransport',
                                                  'subprocessCount')
        self.startYear = startYear
        self.endYear = endYear

        self.masksSubtask = masksSubtask
        self.run_after(masksSubtask)

        self.transectsToPlot = transectsToPlot
        self.restartFileName = None
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
        super(ComputeTransportSubtask, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        # Load mesh related variables
        try:
            self.restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for transport calculations')

        # }}}

    def run_task(self):  # {{{
        """
        Computes time-series of transport through transects.
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Stephen Price

        self.logger.info("Computing time series of transport through "
                         "transects...")

        config = self.config

        startDate = '{:04d}-01-01_00:00:00'.format(self.startYear)
        endDate = '{:04d}-12-31_23:59:59'.format(self.endYear)

        outputDirectory = '{}/transport/'.format(
            build_config_full_path(config, 'output', 'timeseriesSubdirectory'))
        try:
            os.makedirs(outputDirectory)
        except OSError:
            pass

        outFileName = '{}/transport_{:04d}-{:04d}.nc'.format(
            outputDirectory, self.startYear, self.endYear)

        inputFiles = sorted(self.historyStreams.readpath(
            'timeSeriesStatsMonthlyOutput', startDate=startDate,
            endDate=endDate, calendar=self.calendar))

        years, months = get_files_year_month(inputFiles,
                                             self.historyStreams,
                                             'timeSeriesStatsMonthlyOutput')

        variableList = ['timeMonthly_avg_layerThickness']
        with open_mpas_dataset(fileName=inputFiles[0],
                               calendar=self.calendar,
                               startDate=startDate,
                               endDate=endDate) as dsIn:
            if 'timeMonthly_avg_normalTransportVelocity' in dsIn:
                variableList.append('timeMonthly_avg_normalTransportVelocity')
            elif 'timeMonthly_avg_normalGMBolusVelocity' in dsIn:
                variableList = variableList + \
                    ['timeMonthly_avg_normalVelocity',
                     'timeMonthly_avg_normalGMBolusVelocity']
            else:
                self.logger.warning('Cannot compute transport velocity. '
                                    'Using advection velocity.')
                variableList.append('timeMonthly_avg_normalVelocity')

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

        transectMaskFileName = self.masksSubtask.maskFileName

        dsTransectMask = xarray.open_dataset(transectMaskFileName)

        # figure out the indices of the transects to plot
        maskTransectNames = decode_strings(dsTransectMask.transectNames)

        dsMesh = xarray.open_dataset(self.restartFileName)
        dsMesh = dsMesh[['dvEdge', 'cellsOnEdge']]
        dsMesh.load()
        dvEdge = dsMesh.dvEdge
        cellsOnEdge = dsMesh.cellsOnEdge - 1

        transectIndices = []
        transectData = []
        self.logger.info('  Caching transect data...')
        for transect in self.transectsToPlot:
            self.logger.info('    transect: {}'.format(transect))
            try:
                transectIndex = maskTransectNames.index(transect)
            except ValueError:
                self.logger.warning('      Not found in masks. Skipping.')
                continue
            transectIndices.append(transectIndex)

            # select the current transect
            dsMask = dsTransectMask.isel(nTransects=[transectIndex])
            dsMask.load()
            edgeIndices = numpy.flatnonzero(dsMask.transectEdgeMasks.values)
            edgeIndices = edgeIndices[edgeIndices >= 0].astype(int)
            edgeSign = dsMask.transectEdgeMaskSigns.isel(
                nEdges=edgeIndices)

            dv = dvEdge.isel(nEdges=edgeIndices)
            coe = cellsOnEdge.isel(nEdges=edgeIndices)
            transectData.append({'edgeIndices': edgeIndices,
                                 'edgeSign': edgeSign,
                                 'dv': dv,
                                 'coe': coe,
                                 'transect': transect})

        timeDatasets = []
        self.logger.info('  Computing transport...')
        for fileName in inputFiles:
            self.logger.info('    input file: {}'.format(fileName))
            dsTimeSlice = open_mpas_dataset(
                fileName=fileName,
                calendar=self.calendar,
                variableList=variableList,
                startDate=startDate,
                endDate=endDate)
            dsTimeSlice.load()

            transectDatasets = []
            for data in transectData:
                self.logger.info('    transect: {}'.format(data['transect']))

                edgeIndices = data['edgeIndices']
                coe = data['coe']
                edgeSign = data['edgeSign']
                dv = data['dv']
                dsIn = dsTimeSlice.isel(nEdges=edgeIndices)

                # work on data from simulations
                if 'timeMonthly_avg_normalTransportVelocity' in dsIn:
                    vel = dsIn.timeMonthly_avg_normalTransportVelocity
                elif 'timeMonthly_avg_normalGMBolusVelocity' in dsIn:
                    vel = (dsIn.timeMonthly_avg_normalVelocity +
                           dsIn.timeMonthly_avg_normalGMBolusVelocity)
                else:
                    vel = dsIn.timeMonthly_avg_normalVelocity

                # get layer thickness on edges by averaging adjacent cells
                h = 0.5 * dsIn.timeMonthly_avg_layerThickness.isel(
                    nCells=coe).sum(dim='TWO')

                edgeTransport = edgeSign * vel * h * dv

                # convert from m^3/s to Sv
                transport = (constants.m3ps_to_Sv * edgeTransport.sum(
                    dim=['nEdges', 'nVertLevels']))

                dsOut = xarray.Dataset()
                dsOut['transport'] = transport
                dsOut.transport.attrs['units'] = 'Sv'
                dsOut.transport.attrs['description'] = \
                    'Transport through transects'
                transectDatasets.append(dsOut)

            dsOut = xarray.concat(transectDatasets, 'nTransects')
            timeDatasets.append(dsOut)

        # combine data sets into a single data set
        dsOut = xarray.concat(timeDatasets, 'Time')
        dsOut.coords['transectNames'] = dsTransectMask.transectNames.isel(
            nTransects=transectIndices)
        dsOut.coords['year'] = (('Time',), years)
        dsOut['year'].attrs['units'] = 'years'
        dsOut.coords['month'] = (('Time',), months)
        dsOut['month'].attrs['units'] = 'months'
        write_netcdf(dsOut, outFileName)

        # }}}

    # }}}


class CombineTransportSubtask(AnalysisTask):  # {{{
    """
    Combine individual time series into a single data set
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, startYears, endYears):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``TimeSeriesOceanRegions``
            The main task of which this is a subtask

        startYears, endYears : list of int
            The beginning and end of each time series to combine

        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(CombineTransportSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='combineTimeSeries')

        self.startYears = startYears
        self.endYears = endYears
        # }}}

    def run_task(self):  # {{{
        """
        Combine the time series
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        outputDirectory = '{}/transport/'.format(
            build_config_full_path(self.config, 'output',
                                   'timeseriesSubdirectory'))

        outFileName = '{}/transport_{:04d}-{:04d}.nc'.format(
            outputDirectory, self.startYears[0], self.endYears[-1])

        if not os.path.exists(outFileName):
            inFileNames = []
            for startYear, endYear in zip(self.startYears, self.endYears):
                inFileName = '{}/transport_{:04d}-{:04d}.nc'.format(
                    outputDirectory, startYear, endYear)
                inFileNames.append(inFileName)

            ds = xarray.open_mfdataset(inFileNames, combine='nested',
                                       concat_dim='Time', decode_times=False)
            ds.load()
            write_netcdf(ds, outFileName)
        # }}}
    # }}}


class PlotTransportSubtask(AnalysisTask):
    """
    Plots time-series output of transport through transects.

    Attributes
    ----------
    transect : str
        Name of the transect to plot

    transectIndex : int
        The index into the dimension ``nTransects`` of the transect to plot

    controlConfig : ``MpasAnalysisConfigParser``
        The configuration options for the control run (if any)

    """
    # Authors
    # -------
    # Xylar Asay-Davis, Stephen Price

    def __init__(self, parentTask, transect, transectIndex, controlConfig,
                 transportTransectFileName):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        transect : str
            Name of the transect to plot

        transectIndex : int
            The index into the dimension ``nTransects`` of the transect to plot

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(PlotTransportSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='plotTransport_{}'.format(transect.replace(' ', '_')))

        self.transportTransectFileName = transportTransectFileName
        self.transect = transect
        self.transectIndex = transectIndex
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
        super(PlotTransportSubtask, self).setup_and_check()

        self.xmlFileNames = ['{}/transport_{}.xml'.format(
            self.plotsDirectory, self.transect.replace(' ', '_'))]
        # }}}

    def run_task(self):  # {{{
        """
        Plots time-series output of transport through transects.
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Stephen Price

        self.logger.info("\nPlotting time series of transport through "
                         "{}...".format(self.transect))

        self.logger.info('  Load transport data...')

        obsDict = {'Drake Passage': [120, 175],
                   'Tasmania-Ant': [147, 167],
                   'Africa-Ant': None,
                   'Antilles Inflow': [-23.1, -13.7],
                   'Mona Passage': [-3.8, -1.4],
                   'Windward Passage': [-7.2, -6.8],
                   'Florida-Cuba': [30, 33],
                   'Florida-Bahamas': [30, 33],
                   'Indonesian Throughflow': [-21, -11],
                   'Agulhas': [-90, -50],
                   'Mozambique Channel': [-20, -8],
                   'Bering Strait': [0.6, 1.0],
                   'Lancaster Sound': [-1.0, -0.5],
                   'Fram Strait': [-4.7, 0.7],
                   'Davis Strait': [-1.6, -3.6],
                   'Barents Sea Opening': [1.4, 2.6],
                   'Nares Strait': [-1.8, 0.2],
                   'Denmark Strait': None,
                   'Iceland-Faroe-Scotland': None}

        config = self.config
        calendar = self.calendar

        fcAll = read_feature_collection(self.transportTransectFileName)

        fc = FeatureCollection()
        for feature in fcAll.features:
            if feature['properties']['name'] == self.transect:
                fc.add_feature(feature)
                break

        transport, trans_mean, trans_std = self._load_transport(config)

        if self.transect in obsDict:
            bounds = obsDict[self.transect]
        else:
            bounds = None

        plotControl = self.controlConfig is not None

        mainRunName = config.get('runs', 'mainRunName')
        movingAverageMonths = config.getint('timeSeriesTransport',
                                            'movingAverageMonths')

        self.logger.info('  Plotting...')

        transectName = self.transect.replace('_', ' ')
        title = transectName
        thumbnailDescription = transectName

        xLabel = 'Time (yr)'
        yLabel = 'Transport (Sv)'

        filePrefix = 'transport_{}'.format(self.transect.replace(' ', '_'))
        outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

        fields = [transport]
        lineColors = ['k']
        lineWidths = [2.5]
        meanString = 'mean={:.2f} $\pm$ {:.2f}'.format(trans_mean, trans_std)
        if plotControl:
            controlRunName = self.controlConfig.get('runs', 'mainRunName')
            ref_transport, ref_mean, ref_std = \
                self._load_transport(self.controlConfig)
            refMeanString = 'mean={:.2f} $\pm$ {:.2f}'.format(ref_mean, ref_std)
            fields.append(ref_transport)
            lineColors.append('r')
            lineWidths.append(1.2)
            legendText = ['{} ({})'.format(mainRunName, meanString),
                          '{} ({})'.format(controlRunName, refMeanString)]

        else:
            legendText = [mainRunName]
            title = '{} ({})'.format(title, meanString)

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

        fig = timeseries_analysis_plot(config, fields, calendar=calendar,
                                       title=title, xlabel=xLabel,
                                       ylabel=yLabel,
                                       movingAveragePoints=movingAverageMonths,
                                       lineColors=lineColors,
                                       lineWidths=lineWidths,
                                       legendText=legendText,
                                       firstYearXTicks=firstYearXTicks,
                                       yearStrideXTicks=yearStrideXTicks)

        if bounds is not None:
            t = transport.Time.values
            plt.gca().fill_between(t, bounds[0]*numpy.ones_like(t),
                                   bounds[1]*numpy.ones_like(t), alpha=0.3,
                                   label='observations')
            plt.legend(loc='lower left')


        # do this before the inset because otherwise it moves the inset
        # and cartopy doesn't play too well with tight_layout anyway
        plt.tight_layout()

        add_inset(fig, fc, width=2.0, height=2.0)

        savefig(outFileName, config)

        caption = 'Transport through the {} Transect'.format(transectName)
        write_image_xml(
            config=config,
            filePrefix=filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='Transport Time Series',
            groupLink='transporttime',
            thumbnailDescription=thumbnailDescription,
            imageDescription=caption,
            imageCaption=caption)
        # }}}

    def _load_transport(self, config):  # {{{
        """
        Reads transport time series for this transect
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        baseDirectory = build_config_full_path(
            config, 'output', 'timeSeriesSubdirectory')

        startYear = config.getint('timeSeries', 'startYear')
        endYear = config.getint('timeSeries', 'endYear')

        inFileName = '{}/transport/transport_{:04d}-{:04d}.nc'.format(
            baseDirectory, startYear, endYear)

        dsIn = xarray.open_dataset(inFileName)
        transport = dsIn.transport.isel(nTransects=self.transectIndex)
        mean = transport.mean().values
        std = transport.std().values
        return transport, mean, std  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
