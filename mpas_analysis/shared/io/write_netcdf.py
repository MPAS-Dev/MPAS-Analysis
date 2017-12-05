'''
Functions for writing data sets

Functions
---------
write_netcdf - write an xarray data set to a NetCDF file using finite fill
    values

Authors
-------
Xylar Asay-Davis

'''

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import netCDF4
import numpy


def write_netcdf(ds, fileName, fillValues=netCDF4.default_fillvals):  # {{{
    '''
    Write an xarray data set to a NetCDF file using finite fill values

    Parameters
    ----------
    ds : xarray.Dataset object
        The xarray data set to be written to a file

    fileName : str
        The fileName to write the data set to

    fillValues : dict
        A dictionary of fill values for each supported data type.  By default,
        this is the dictionary used by the netCDF4 package.  Key entries should
        be of the form 'f8' (for float64), 'i4' (for int32), etc.

    Authors
    -------
    Xylar Asay-Davis

    '''
    encodingDict = {}
    variableNames = list(ds.data_vars.keys()) + list(ds.coords.keys())
    for variableName in variableNames:
        dtype = ds[variableName].dtype
        for fillType in fillValues:
            if dtype == numpy.dtype(fillType):
                encodingDict[variableName] = \
                    {'_FillValue': fillValues[fillType]}
                break

    ds.to_netcdf(fileName, encoding=encodingDict)

    # }}}

# vim: ai ts=4 sts=4 et sw=4 ft=python
