
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
import numpy
import xarray


def add_periodic_lon(ds, lonDim, degrees=True):
    """
    Add a single grid point that is a periodic image to the end of the data set
    in the lon dimension
    """

    lonRange = ds[lonDim][-1].values - ds[lonDim][0].values
    if degrees:
        period = 360.
    else:
        period = 2. * numpy.pi

    if numpy.abs(lonRange - period) < 1e-10:
        # already periodic
        return

    nLon = ds.sizes[lonDim]
    lonIndices = xarray.DataArray(numpy.append(numpy.arange(nLon), [0]),
                                  dims=('newLon',))
    ds.load()
    ds = ds.isel({lonDim: lonIndices})
    ds = ds.rename({'newLon': lonDim})

    # fix the last entry in lon
    lon = ds[lonDim].values
    lon[-1] += period
    dims = ds[lonDim].dims
    attrs = ds[lonDim].attrs
    ds[lonDim] = (dims, lon)
    ds[lonDim].attrs = attrs

    return ds
