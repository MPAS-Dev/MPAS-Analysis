"""
Functions for creating climatologies from monthly time series data

Authors
-------
Xylar Asay-Davis

Last Modified
-------------
02/26/2017
"""

import numpy as np
import netCDF4

from ..interpolation.interpolate import interp_fields, init_tree
from ..constants import constants

from ..timekeeping.utility import days_to_datetime


def compute_monthly_climatology(ds, calendar):
    """
    Compute a monthly climatology data set from a data set with Time expressed
    as days since 0001-01-01 with the given calendar.

    Parameters
    ----------
    ds : an xarray data set with a 'Time' coordinate expressed as days since
        0001-01-01

    calendar: {'gregorian', 'gregorian_noleap'}
        The name of one of the calendars supported by MPAS cores

    Returns
    -------
    monthlyClimatology : an xarray data set with a new 'month' coordinate,
        containing monthly climatologies of all variables in ds

    Authors
    -------
    Luke Van Roekel, Milena Veneziani, Xylar Asay-Davis

    Last Modified
    -------------
    02/25/2017
    """
    months = [date.month for date in days_to_datetime(ds.Time,
                                                      calendar=calendar)]

    ds.coords['month'] = ('Time', months)
    monthlyClimatology = ds.groupby('month').mean('Time')
    return monthlyClimatology


def init_model_interpolation(mpasMeshFileName):
    """
    Given an MPAS mesh file, computes weight information needed to perform
    climatology interpolation to a comparison grid.  The resolution of the
    comparison grid is determined via the constants module.

    Parameters
    ----------
    mpasMeshFileName : an MPAS file contining mesh information

    Returns
    -------
    d : the distance to the nearest MPAS cell center from a given comparison
        grid point

    inds : the index of the nearest MPAS cell center to a given comparison
        grid point

    lonTarg, latTarg : the longitude and latitude coordinates (as 2D arrays)
        of points on the default comparison grid

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    02/26/2017
    """

    ncFile = netCDF4.Dataset(mpasMeshFileName, mode='r')
    lonCell = ncFile.variables["lonCell"][:]
    latCell = ncFile.variables["latCell"][:]
    ncFile.close()

    return init_tree(np.rad2deg(lonCell),
                     np.rad2deg(latCell),
                     constants.lonmin,
                     constants.lonmax,
                     constants.latmin,
                     constants.latmax,
                     constants.dLongitude,
                     constants.dLatitude)


def init_observations_interpolation(dsObs):
    """
    Given monthly climatology from observations on a lat/lon grid, interpolates
    the desired field over the desired months to a the default comparison grid.

    Parameters
    ----------
    dsObs : an xarray data set containing longitude and latitude information
        for the observations

    Returns
    -------
    d : the distance to the nearest MPAS cell center from a given comparison
        grid point

    inds : the index of the nearest MPAS cell center to a given comparison
        grid point

    lonTarg, latTarg : the longitude and latitude coordinates (as 2D arrays)
        of points on the default comparison grid

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    02/26/2017
    """

    latData, lonData = np.meshgrid(dsObs.lat.values,
                                   dsObs.lon.values)
    latData = latData.flatten()
    lonData = lonData.flatten()

    # initialize interpolation variables
    return init_tree(lonData, latData,
                     constants.lonmin,
                     constants.lonmax,
                     constants.latmin,
                     constants.latmax,
                     constants.dLongitude,
                     constants.dLatitude)


def interpolate_model_climatology(interpolationData, monthlyClimatology, field,
                                  monthsValue):
    """
    Given a monthly climatology on the mpas grid, together with corresponding
    interpolation weighting data, interpolates desired field over the desired
    months to a the default comparison grid.

    Parameters
    ----------
    interpolationData : a tuple returned by init_model_interpolation that
        contains the information needed to interpolate from the MPAS mesh
        to the comparison grid

    monthlyClimatology : an xarray data set containing a monthly climatology

    field : the name of a field in monthlyClimatology to be interpolated

    monthsValue : a single month or an array of months to be averaged together
        before interpolation

    Returns
    -------
    modelOutput : a numpy array containing the field interpolated to the
        default comparison grid

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    02/26/2017
    """

    d, inds, lonTarg, latTarg = interpolationData

    nLon = lonTarg.shape[0]
    nLat = lonTarg.shape[1]

    modelOutput = np.zeros((nLon, nLat))

    if isinstance(monthsValue, (int, long)):
        modelData = monthlyClimatology.sel(month=monthsValue)[field].values
    else:

        modelData = (np.sum(
            constants.daysInMonth[monthsValue-1] *
            monthlyClimatology.sel(month=monthsValue)[field].values.T,
            axis=1) /
            np.sum(constants.daysInMonth[monthsValue-1]))

    modelOutput = interp_fields(modelData, d, inds, lonTarg)

    return modelOutput


def interpolate_observation_climatology(interpolationData, dsObs,
                                        obsFieldName, monthsValue):
    """
    Given monthly climatology from observations on a lat/lon grid, interpolates
    the desired field over the desired months to a the default comparison grid.

    Parameters
    ----------
    interpolationData : a tuple returned by init_observations_interpolation
        that contains the information needed to interpolate from the
        observation grid to the comparison grid

    dsObs : an xarray data set containing a monthly climatology from
        observations

    obsFieldName : the name of a field in dsObs to be interpolated

    monthsValue : a single month or an array of months to be averaged together
        before interpolation

    Returns
    -------
    observations : a numpy array containing the field interpolated to the
        default comparison grid

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    02/26/2017
    """

    latData, lonData = np.meshgrid(dsObs.lat.values,
                                   dsObs.lon.values)
    latData = latData.flatten()
    lonData = lonData.flatten()

    daysarray = np.ones((12,
                         dsObs[obsFieldName].values.shape[1],
                         dsObs[obsFieldName].values.shape[2]))

    for monthIndex, dval in enumerate(constants.daysInMonth):
        daysarray[monthIndex, :, :] = dval
        inds = np.where(np.isnan(
                dsObs[obsFieldName][monthIndex, :, :].values))
        daysarray[monthIndex, inds[0], inds[1]] = np.NaN

    d, inds, lonTarg, latTarg = interpolationData

    nLon = lonTarg.shape[0]
    nLat = latTarg.shape[1]

    observations = np.zeros((nLon, nLat))

    # Interpolate and compute biases

    if isinstance(monthsValue, (int, long)):
        obsData = dsObs.sel(month=monthsValue)[obsFieldName].values
    else:
        obsData = \
            (np.nansum(
                    daysarray[monthsValue-1, :, :] *
                    dsObs.sel(month=monthsValue)[obsFieldName].values,
                    axis=0) /
             np.nansum(daysarray[monthsValue-1, :, :], axis=0))

    observations = interp_fields(obsData.flatten(), d, inds, lonTarg)

    return observations

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
