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
import numpy
import matplotlib.pyplot as plt

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.ocean.compute_anomaly_subtask import ComputeAnomalySubtask
from mpas_analysis.shared.io import open_mpas_dataset
from mpas_analysis.shared.io.utility import build_config_full_path
from mpas_analysis.shared.plot import timeseries_analysis_plot, savefig
from mpas_analysis.shared.html import write_image_xml


class TimeSeriesSSHAnomaly(AnalysisTask):
    """
    Plots a time series of the global mean sea surface height relative to a
    reference time.

    Attributes
    ----------
    variableDict : dict
        A dictionary of variables from the time series stats monthly output
        (keys), together with shorter, more convenient names (values)

    timeSeriesFileName : str
        The name of the file where the ssh anomaly is stored

    controlConfig :  mpas_analysis.configuration.MpasAnalysisConfigParser
        Configuration options for a control run (if one is provided)

    filePrefix : str
        The basename (without extension) of the PNG and XML files to write out
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasTimeSeriesTask, controlConfig):

        """
        Construct the analysis task.

        Parameters
        ----------
        config :  mpas_analysis.configuration.MpasAnalysisConfigParser
            Configuration options

        mpasTimeSeriesTask : mpas_analysis.shared.time_series.MpasTimeSeriesTask
            The task that extracts the time series from MPAS monthly output

        controlConfig :  mpas_analysis.configuration.MpasAnalysisConfigParser
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super().__init__(
            config=config,
            taskName='timeSeriesSSHAnomaly',
            componentName='ocean',
            tags=['timeSeries', 'ssh', 'publicObs', 'anomaly'])

        self.controlConfig = controlConfig
        self.filePrefix = None

        movingAveragePoints = config.getint('timeSeriesSSHAnomaly',
                                            'movingAveragePoints')

        self.variableDict = {}
        for var in ['volumeCellGlobal', 'areaCellGlobal']:
            key = 'timeMonthly_avg_{}'.format(var)
            self.variableDict[key] = var

        self.timeSeriesFileName = 'globalMeanSSHAnomaly.nc'

        anomalyTask = ComputeAnomalySubtask(
            parentTask=self,
            mpasTimeSeriesTask=mpasTimeSeriesTask,
            outFileName=self.timeSeriesFileName,
            variableList=list(self.variableDict.keys()),
            movingAveragePoints=movingAveragePoints,
            alter_dataset=self._compute_mean_water_column)
        self.add_subtask(anomalyTask)

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        OSError
            If files are not present
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #   self.inDirectory, self.plotsDirectory, self.namelist, self.streams
        #   self.calendar
        super().setup_and_check()

        config = self.config

        mainRunName = config.get('runs', 'mainRunName')

        self.xmlFileNames = []

        self.filePrefix = 'ssh_anomaly_global_{}'.format(mainRunName)
        self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                    self.filePrefix))

    def run_task(self):
        """
        Performs analysis of the time-series output of sea-surface temperature
        (SST).
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Milena Veneziani

        self.logger.info("\nPlotting time series of SSH anomaly...")

        self.logger.info('  Load SSH anomaly data...')

        config = self.config
        calendar = self.calendar

        startDate = self.config.get('timeSeries', 'startDate')
        endDate = self.config.get('timeSeries', 'endDate')

        mainRunName = config.get('runs', 'mainRunName')

        variableList = ['meanWaterColumnThickness']

        baseDirectory = build_config_full_path(
            config, 'output', 'timeSeriesSubdirectory')
        fileName = '{}/{}'.format(baseDirectory, self.timeSeriesFileName)
        ds = open_mpas_dataset(fileName=fileName,
                               calendar=calendar,
                               variableList=variableList,
                               timeVariableNames=None,
                               startDate=startDate,
                               endDate=endDate)

        if self.controlConfig is not None:
            baseDirectory = build_config_full_path(
                self.controlConfig, 'output', 'timeSeriesSubdirectory')

            controlFileName = '{}/{}'.format(baseDirectory,
                                             self.timeSeriesFileName)

            controlStartYear = self.controlConfig.getint(
                'timeSeries', 'startYear')
            controlEndYear = self.controlConfig.getint('timeSeries', 'endYear')
            controlStartDate = '{:04d}-01-01_00:00:00'.format(controlStartYear)
            controlEndDate = '{:04d}-12-31_23:59:59'.format(controlEndYear)

            dsRef = open_mpas_dataset(
                fileName=controlFileName,
                calendar=calendar,
                variableList=variableList,
                timeVariableNames=None,
                startDate=controlStartDate,
                endDate=controlEndDate)
        else:
            dsRef = None

        self.logger.info('  Make plots...')
        title = 'Global Mean SSH Anomaly'
        xLabel = 'Time (years)'
        yLabel = 'SSH Anomaly (m)'

        varName = variableList[0]
        deltaSSH = ds[varName]

        outFileName = '{}/{}.png'.format(self.plotsDirectory, self.filePrefix)

        lineColors = [config.get('timeSeries', 'mainColor')]
        lineWidths = [3]

        fields = [deltaSSH]
        legendText = [mainRunName]

        if dsRef is not None:
            refDeltaSSH = dsRef[varName]
            fields.append(refDeltaSSH)
            lineColors.append(config.get('timeSeries', 'controlColor'))
            lineWidths.append(1.5)
            controlRunName = self.controlConfig.get('runs', 'mainRunName')
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

        timeseries_analysis_plot(config, fields, calendar=calendar,
                                 title=title, xlabel=xLabel, ylabel=yLabel,
                                 lineColors=lineColors,
                                 lineWidths=lineWidths,
                                 legendText=legendText,
                                 firstYearXTicks=firstYearXTicks,
                                 yearStrideXTicks=yearStrideXTicks)

        if config.has_option(self.taskName, 'fitStartYear') and \
                config.has_option(self.taskName, 'fitEndYear'):
            fitStartYear = config.getint(self.taskName, 'fitStartYear')
            fitEndYear = config.getint(self.taskName, 'fitEndYear')
            fitStartDate = '{:04d}-01-01_00:00:00'.format(fitStartYear)
            fitEndDate = '{:04d}-12-31_23:59:59'.format(fitEndYear)

            dsFit = open_mpas_dataset(fileName=fileName, calendar=calendar,
                                      variableList=variableList,
                                      timeVariableNames=None,
                                      startDate=fitStartDate,
                                      endDate=fitEndDate)

            time = dsFit.Time.values
            deltaSSH = dsFit[varName].values

            coeffs = numpy.polyfit(time, deltaSSH, deg=1)
            print(coeffs)
            line = numpy.poly1d(coeffs)
            ylim = plt.gca().get_ylim()

            # m/day --> m/yr
            slope = coeffs[0]*365

            label = r'linear fit {:04d}-{:04d}: slope = ${}$ m/yr'.format(
                fitStartYear, fitEndYear, _as_si(slope, 2))

            fitColor = config.get('timeSeries', 'fitColor1')
            full_time = ds.Time.values
            plt.plot(full_time, line(full_time), color=fitColor,
                     linewidth=0.6*lineWidths[0],
                     label=label)

            plt.plot([time[0], time[-1]], line([time[0], time[-1]]), '.',
                     color=fitColor, markersize=16)
            plt.gca().set_ylim(ylim)
            plt.legend(loc='best')

        savefig(outFileName, config)

        caption = title
        write_image_xml(
            config=config,
            filePrefix=self.filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='Time Series',
            groupLink='timeseries',
            thumbnailDescription='SSH Anomaly',
            imageDescription=caption,
            imageCaption=caption)

    def _compute_mean_water_column(self, ds):
        """
        Compute a time series of the global mean water-column thickness.
        """

        # for convenience, rename the variables to simpler, shorter names
        ds = ds.rename(self.variableDict)

        ds['meanWaterColumnThickness'] = \
            ds['volumeCellGlobal'] / ds['areaCellGlobal']
        ds.meanWaterColumnThickness.attrs['units'] = 'm'
        ds.meanWaterColumnThickness.attrs['description'] = \
            'Global mean thickness of the water column'

        return ds


def _as_si(x, ndp):
    """
    https://stackoverflow.com/a/31453961/7728169
    """
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))
