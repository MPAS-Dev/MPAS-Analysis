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
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
#

from mpas_analysis.shared.io import open_mpas_dataset
from mpas_analysis.shared.time_series.moving_average import compute_moving_avg


def compute_moving_avg_anomaly_from_start(timeSeriesFileName, variableList,
                                          anomalyStartTime, anomalyEndTime,
                                          startDate, endDate, calendar,
                                          movingAveragePoints=12,
                                          alter_dataset=None):
    """
    Compute the rolling mean of the anomaly of a quantity from the beginning
    of the simulation (such that the rolling mean starts at zero by definition)

    Parameters
    ----------
    timeSeriesFileName :  str
        a file produced by ``MpasTimeSeriesTask`` containing variables, the
        anomaly and rolling mean of which is to be computed

    variableList :  list of str
        variable names to include in the resulting data set

    anomalyStartTime, anomalyEndTime :  str
        the start and end times of the reference point for the anomaly

    startDate, endDate :  str
        the start and end dates of the time series

    calendar : {'gregorian', 'gregoraian_noleap'}
        The calendar used in the MPAS run

    movingAveragePoints : int, optional
        The number of points (months) over which to perform the rolling average
        of the data set

    alter_dataset : function
        A function for manipulating the data set (e.g. computing new
        variables), taking an ``xarray.Dataset`` as input argument and
        returning an ``xarray.Dataset``

    Returns
    -------
    ds : ``xarray.Dataset``
        The anomaly of the rolling time mean from the start of the simulation
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    ds = open_mpas_dataset(fileName=timeSeriesFileName,
                           calendar=calendar,
                           variableList=variableList,
                           startDate=startDate,
                           endDate=endDate)

    if alter_dataset is not None:
        ds = alter_dataset(ds)

    dsStart = open_mpas_dataset(
        fileName=timeSeriesFileName,
        calendar=calendar,
        variableList=variableList,
        startDate=anomalyStartTime,
        endDate=anomalyEndTime)

    if alter_dataset is not None:
        dsStart = alter_dataset(dsStart)

    dsStart = dsStart.isel(Time=slice(0, movingAveragePoints)).mean('Time')

    for variable in ds.data_vars:
        ds[variable] = ds[variable] - dsStart[variable]

    ds = compute_moving_avg(ds, movingAveragePoints=movingAveragePoints)

    return ds

