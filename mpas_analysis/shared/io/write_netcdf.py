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
import netCDF4
import numpy


def write_netcdf_with_fill(ds, fileName, fillValues=netCDF4.default_fillvals):
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
        isNumeric = numpy.issubdtype(ds[variableName].dtype,
                                     numpy.number)
        if isNumeric and numpy.any(numpy.isnan(ds[variableName])):
            dtype = ds[variableName].dtype
            for fillType in fillValues:
                if dtype == numpy.dtype(fillType):
                    encodingDict[variableName] = \
                        {'_FillValue': fillValues[fillType]}
                    break
        else:
            encodingDict[variableName] = {'_FillValue': None}

        # make strings write as unicode instead
        if dtype.type is numpy.bytes_:
            encodingDict[variableName] = {'dtype': str}

    unlimited_dims = ds.encoding.get('unlimited_dims', None)
    if unlimited_dims is not None:
        if isinstance(unlimited_dims, str):
            unlimited_dims = {unlimited_dims}
        unlimited_dims = [dim for dim in unlimited_dims if dim in ds.dims]
        ds.encoding['unlimited_dims'] = set(unlimited_dims)
    ds.to_netcdf(fileName, encoding=encodingDict)
