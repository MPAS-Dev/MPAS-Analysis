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
import netCDF4
import numpy


def write_netcdf(ds, fileName, fillValues=netCDF4.default_fillvals):
    """
    Write an xarray data set to a NetCDF file using finite fill values and
    unicode strings

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
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    encodingDict = {}
    variableNames = list(ds.data_vars.keys()) + list(ds.coords.keys())
    for variableName in variableNames:
        dtype = ds[variableName].dtype

        # add fill values for types that have them
        for fillType in fillValues:
            if dtype == numpy.dtype(fillType):
                encodingDict[variableName] = \
                    {'_FillValue': fillValues[fillType]}
                break

        # make strings write as unicode instead
        if dtype.type is numpy.string_:
            encodingDict[variableName] = {'dtype': str}

    ds.to_netcdf(fileName, encoding=encodingDict)


# vim: ai ts=4 sts=4 et sw=4 ft=python
