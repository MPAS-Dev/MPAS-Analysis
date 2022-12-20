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

import os

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.io import write_netcdf

from mpas_analysis.shared.timekeeping.utility import \
    get_simulation_start_time, string_to_datetime

from mpas_analysis.shared.timekeeping.MpasRelativeDelta import \
    MpasRelativeDelta

from mpas_analysis.shared.io.utility import build_config_full_path

from mpas_analysis.shared.time_series import \
    compute_moving_avg_anomaly_from_start


class ComputeAnomalySubtask(AnalysisTask):
    """
    A subtask for computing anomalies of moving averages and writing them out.

    Attributes
    ----------

    mpasTimeSeriesTask : ``MpasTimeSeriesTask``
        The task that extracts the time series from MPAS monthly output

    outFileName : str
        The file name (usually without full path) where the resulting
        data set should be written

    variableList : list of str
        Variables to be included in the data set

    movingAveragePoints : int
        The number of points (months) used in the moving average used to
        smooth the data set

    alter_dataset : function
        A function that takes an ``xarray.Dataset`` and returns an
        ``xarray.Dataset`` for manipulating the data set (e.g. adding a new
        variable computed from others).  This operation is performed before
        computing moving averages and anomalies, so that these operations are
        also performed on any new variables added to the data set.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, mpasTimeSeriesTask, outFileName,
                 variableList, movingAveragePoints,
                 subtaskName='computeAnomaly', alter_dataset=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``AnalysisTask``
            The parent task of which this is a subtask

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

        outFileName : str
            The file name (usually without full path) where the resulting
            data set should be written

        variableList : list of str
            Variables to be included in the data set

        movingAveragePoints : int
            The number of points (months) used in the moving average used to
            smooth the data set

        subtaskName :  str, optional
            The name of the subtask

        alter_dataset : function
            A function that takes an ``xarray.Dataset`` and returns an
            ``xarray.Dataset`` for manipulating the data set (e.g. adding a new
            variable computed from others).  This operation is performed before
            computing moving averages and anomalies, so that these operations
            are also performed on any new variables added to the data set.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(ComputeAnomalySubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName='ocean',
            tags=parentTask.tags,
            subtaskName=subtaskName)

        self.mpasTimeSeriesTask = mpasTimeSeriesTask

        self.run_after(mpasTimeSeriesTask)

        self.outFileName = outFileName
        self.variableList = variableList
        self.movingAveragePoints = movingAveragePoints

        self.alter_dataset = alter_dataset

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
        super(ComputeAnomalySubtask, self).setup_and_check()

        startDate = self.config.get('timeSeries', 'startDate')
        endDate = self.config.get('timeSeries', 'endDate')

        delta = MpasRelativeDelta(string_to_datetime(endDate),
                                  string_to_datetime(startDate),
                                  calendar=self.calendar)

        months = delta.months + 12*delta.years

        if months <= self.movingAveragePoints:
            raise ValueError('Cannot meaninfully perform a rolling mean '
                             'because the time series is too short.')

        self.mpasTimeSeriesTask.add_variables(variableList=self.variableList)

        self.inputFile = self.mpasTimeSeriesTask.outputFile

    def run_task(self):
        """
        Performs analysis of ocean heat content (OHC) from time-series output.
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Milena Veneziani, Greg Streletz

        self.logger.info("\nComputing anomalies...")

        config = self.config
        startDate = config.get('timeSeries', 'startDate')
        endDate = config.get('timeSeries', 'endDate')

        if config.has_option('timeSeries', 'anomalyRefYear'):
            anomalyYear = config.getint('timeSeries', 'anomalyRefYear')
            anomalyRefDate = '{:04d}-01-01_00:00:00'.format(anomalyYear)
            anomalyEndDate = '{:04d}-12-31_23:59:59'.format(anomalyYear)
        else:
            anomalyRefDate = get_simulation_start_time(self.runStreams)
            anomalyYear = int(anomalyRefDate[0:4])
            anomalyEndDate = '{:04d}-12-31_23:59:59'.format(anomalyYear)

        ds = compute_moving_avg_anomaly_from_start(
            timeSeriesFileName=self.inputFile,
            variableList=self.variableList,
            anomalyStartTime=anomalyRefDate,
            anomalyEndTime=anomalyEndDate,
            startDate=startDate,
            endDate=endDate,
            calendar=self.calendar,
            movingAveragePoints=self.movingAveragePoints,
            alter_dataset=self.alter_dataset)

        outFileName = self.outFileName
        if not os.path.isabs(outFileName):
            baseDirectory = build_config_full_path(
                config, 'output', 'timeSeriesSubdirectory')

            outFileName = '{}/{}'.format(baseDirectory,
                                         outFileName)

        write_netcdf(ds, outFileName)
