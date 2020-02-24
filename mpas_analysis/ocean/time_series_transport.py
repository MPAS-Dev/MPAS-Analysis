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
    decode_strings, get_region_mask

from mpas_analysis.shared.html import write_image_xml

from mpas_analysis.shared.transects import ComputeTransectMasksSubtask

from mpas_analysis.shared.regions import get_feature_list


class TimeSeriesTransport(AnalysisTask):  # {{{
    """
    Extract and plot time series of transport through transects on the MPAS
    mesh.
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
        super(TimeSeriesTransport, self).__init__(
            config=config,
            taskName='timeSeriesTransport',
            componentName='ocean',
            tags=['timeSeries', 'transport'])

        transportTransectFileName = \
            get_region_mask(config, 'transportTransects.geojson')

        transectsToPlot = config.getExpression('timeSeriesTransport',
                                               'transectsToPlot')
        if 'all' in transectsToPlot:
            transectsToPlot = get_feature_list(transportTransectFileName)

        masksSubtask = ComputeTransectMasksSubtask(
            self, transportTransectFileName, outFileSuffix='transportMasks')

        self.add_subtask(masksSubtask)

        computeTransportSubtask = ComputeTransportSubtask(
            self, mpasTimeSeriesTask, masksSubtask, transectsToPlot)
        self.add_subtask(computeTransportSubtask)

        for index, transect in enumerate(transectsToPlot):
            plotTransportSubtask = PlotTransportSubtask(
                self, transect, index, controlConfig, transportTransectFileName)
            plotTransportSubtask.run_after(computeTransportSubtask)
            self.add_subtask(plotTransportSubtask)

        # }}}

    # }}}


class ComputeTransportSubtask(AnalysisTask):  # {{{
    """
    Computes time-series of transport through transects.

    Attributes
    ----------
    mpasTimeSeriesTask : ``MpasTimeSeriesTask``
        The task that extracts the time series from MPAS monthly output

    masksSubtask : ``ComputeRegionMasksSubtask``
        A task for creating mask files for each ice shelf to plot

    transectsToPlot : list of str
        A list of transects to plot
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Stephen Price

    def __init__(self, parentTask, mpasTimeSeriesTask, masksSubtask,
                 transectsToPlot):  # {{{
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
            subtaskName='computeTransport')

        self.mpasTimeSeriesTask = mpasTimeSeriesTask
        self.run_after(mpasTimeSeriesTask)

        self.masksSubtask = masksSubtask
        self.run_after(masksSubtask)

        self.transectsToPlot = transectsToPlot

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
        super(ComputeTransportSubtask, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        config = self.config

        # Load mesh related variables
        try:
            self.restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for transport calculations')

        # get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        self.startDate = config.get('timeSeries', 'startDate')
        self.endDate = config.get('timeSeries', 'endDate')

        try:
            self.variableList = \
                ['timeMonthly_avg_normalTransportVelocity',
                 'timeMonthly_avg_layerThickness']
            self.mpasTimeSeriesTask.add_variables(
                variableList=self.variableList)
        except ValueError:
            try:
                self.variableList = \
                    ['timeMonthly_avg_normalVelocity',
                     'timeMonthly_avg_normalGMBolusVelocity'
                     'timeMonthly_avg_layerThickness']
                self.mpasTimeSeriesTask.add_variables(
                    variableList=self.variableList)
            except ValueError:
                print('Warning: computing transport with advection velocity '
                      'without GM contribution')
                self.variableList = \
                    ['timeMonthly_avg_normalVelocity',
                     'timeMonthly_avg_layerThickness']
                self.mpasTimeSeriesTask.add_variables(
                    variableList=self.variableList)

        return  # }}}

    def run_task(self):  # {{{
        """
        Computes time-series of transport through transects.
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Stephen Price

        self.logger.info("Computing time series of transport through "
                         "transects...")

        self.logger.info('  Load transport velocity data...')

        mpasTimeSeriesTask = self.mpasTimeSeriesTask
        config = self.config

        baseDirectory = build_config_full_path(
            config, 'output', 'timeSeriesSubdirectory')

        outFileName = '{}/transectTransport.nc'.format(baseDirectory)

        # Load data:
        inputFile = mpasTimeSeriesTask.outputFile
        dsIn = open_mpas_dataset(fileName=inputFile,
                                 calendar=self.calendar,
                                 variableList=self.variableList,
                                 startDate=self.startDate,
                                 endDate=self.endDate)

        dsIn = dsIn.chunk({'Time': 12})
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
            # something is potentially wrong with the file, so let's delete
            # it and try again
            self.logger.warning('Problems reading file {}. Deleting '
                                'it.'.format(outFileName))
            os.remove(outFileName)

        with dask.config.set(schedular='threads',
                             pool=ThreadPool(self.daskThreads)):

            dsMesh = xarray.open_dataset(self.restartFileName)
            dvEdge = dsMesh.dvEdge

            # work on data from simulations
            if 'timeMonthly_avg_normalTransportVelocity' in dsIn:
                vel = dsIn.timeMonthly_avg_normalTransportVelocity
            elif 'timeMonthly_avg_normalGMBolusVelocity' in dsIn:
                vel = (dsIn.timeMonthly_avg_normalVelocity +
                       dsIn.timeMonthly_avg_normalGMBolusVelocity)
            else:
                vel = dsIn.timeMonthly_avg_normalVelocity

            # get layer thickness on edges by averaging adjacent cells
            h = 0.5*dsIn.timeMonthly_avg_layerThickness.isel(
                nCells=(dsMesh.cellsOnEdge-1)).sum(dim='TWO')

            transectMaskFileName = self.masksSubtask.maskFileName

            dsTransectMask = xarray.open_dataset(transectMaskFileName)

            # figure out the indices of the transects to plot
            transectNames = decode_strings(dsTransectMask.transectNames)

            transectIndices = []
            outTransectNames = []
            for transect in self.transectsToPlot:
                found = False
                for index, otherName in enumerate(transectNames):
                    if transect == otherName:
                        transectIndices.append(index)
                        outTransectNames.append(transect)
                        found = True
                        break
                if not found:
                    self.logger.warning('transect {} was not found in transect '
                                        'masks'.format(transect))

            # select only those transects we want to plot
            dsTransectMask = dsTransectMask.isel(nTransects=transectIndices)
            edgeSign = dsTransectMask.transectEdgeMaskSigns.chunk(
                {'nTransects': 5})

            # convert from m^3/s to Sv
            transport = (constants.m3ps_to_Sv *
                         (edgeSign * vel * h * dvEdge).sum(
                             dim=['nEdges', 'nVertLevels']))
            self.logger.info(transport)
            transport.compute()

            dsOut = xarray.Dataset()
            dsOut['transport'] = transport
            dsOut.transport.attrs['units'] = 'Sv'
            dsOut.transport.attrs['description'] = \
                'Transport through transects'
            dsOut['transectNames'] = ('nTransects', outTransectNames)

            write_netcdf(dsOut, outFileName)

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

        config = self.config
        calendar = self.calendar

        fcAll = read_feature_collection(self.transportTransectFileName)

        fc = FeatureCollection()
        for feature in fcAll.features:
            if feature['properties']['name'] == self.transect:
                fc.add_feature(feature)
                break

        transport = self._load_transport(config)

        plotControl = self.controlConfig is not None

        mainRunName = config.get('runs', 'mainRunName')
        movingAverageMonths = config.getint('timeSeriesTransport',
                                            'movingAverageMonths')

        self.logger.info('  Plotting...')

        title = self.transect.replace('_', ' ')

        xLabel = 'Time (yr)'
        yLabel = 'Transport (Sv)'

        filePrefix = 'transport_{}'.format(self.transect.replace(' ', '_'))
        outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

        fields = [transport]
        lineColors = ['k']
        lineWidths = [2.5]
        legendText = [mainRunName]
        if plotControl:
            controlRunName = self.controlConfig.get('runs', 'mainRunName')
            refTransport = self._load_transport(self.controlConfig)
            fields.append(refTransport)
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
                                       firstYearXTicks=firstYearXTicks,
                                       yearStrideXTicks=yearStrideXTicks)

        # do this before the inset because otherwise it moves the inset
        # and cartopy doesn't play too well with tight_layout anyway
        plt.tight_layout()

        add_inset(fig, fc, width=2.0, height=2.0)

        savefig(outFileName)

        caption = 'Transport through the {} Transect'.format(title)
        write_image_xml(
            config=config,
            filePrefix=filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='Transport Time Series',
            groupLink='transporttime',
            thumbnailDescription=title,
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

        inFileName = '{}/transectTransport.nc'.format(baseDirectory)

        dsIn = xarray.open_dataset(inFileName)
        return dsIn.transport.isel(nTransects=self.transectIndex)  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
