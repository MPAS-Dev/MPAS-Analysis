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


def interp_1d(ds, inInterpDim, inInterpCoord, outInterpDim,
              outInterpCoord):
    """
    Interpolate 1D or 2D fields in 1D

    Parameters
    ----------
    ds : ``xarray.Dataset```
        The data set containing variables before 1D interpolation

    inInterpDim, outInterpDim : str
        The name of the dimensions to interpolate before and after
        interpolation

    inInterpCoord, outInterpCoord : str
        The name of the coordinates to interpolate from and to.  Each
        of these can be 1D (vertical) or 2D fields
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    indices, weight0 = _compute_weights_and_indices(
        ds, inInterpDim, inInterpCoord, outInterpDim, outInterpCoord)

    # conert coords to normal data variables
    coords = list(ds.coords)
    ds = ds.reset_coords(coords)

    ds = ds.map(_interp_1d_array, args=(indices, weight0, inInterpDim))

    # conert back to coords
    ds = ds.set_coords(coords)

    return ds


def _compute_weights_and_indices(ds, inInterpDim, inInterpCoord, outInterpDim,
                                 outInterpCoord):
    """
    add interpolation weights and indices to the data set
    """

    xIn = ds[inInterpCoord]
    xOut = ds[outInterpCoord]

    outDims = xOut.dims
    inDims = xIn.dims

    inAxis = inDims.index(inInterpDim)
    outAxis = outDims.index(outInterpDim)

    inInterpSize = ds.sizes[inInterpDim]
    outInterpSize = ds.sizes[outInterpDim]

    xIn = xIn.values
    xOut = xOut.values

    # we need to add size-zero dimensions to xIn and/or xOut for any dimensions
    # that need to be "broadcast" (i.e. if they are in xIn but not xOut or
    # visa versa).
    allOutDims = [outDims[axis] for axis in range(outAxis)]
    for axis in range(inAxis):
        dim = inDims[axis]
        if dim not in allOutDims:
            allOutDims.append(dim)

    allOutDims.append(outInterpDim)

    allOutDims.extend(outDims[outAxis + 1:])

    for axis in range(inAxis + 1, len(inDims)):
        dim = inDims[axis]
        if dim not in allOutDims:
            allOutDims.append(dim)

    outSizes = [ds.sizes[d] for d in allOutDims]

    allInDims = list(allOutDims)
    index = allInDims.index(outInterpDim)
    allInDims.pop(index)
    allInDims.insert(index, inInterpDim)

    shape = list(outSizes)
    for index, dim in enumerate(allOutDims):
        if dim not in outDims:
            shape[index] = 1

    xOut = xOut.reshape(shape)

    shape = [ds.sizes[d] for d in allInDims]
    for index, dim in enumerate(allInDims):
        if dim not in inDims:
            shape[index] = 1

    xIn = xIn.reshape(shape)

    inAxis = allInDims.index(inInterpDim)
    outAxis = allOutDims.index(outInterpDim)

    indexArrays = numpy.indices(outSizes, int)
    indices = {}
    for index, dim in enumerate(allInDims):
        indices[dim] = indexArrays[index]
    index0 = indices[inInterpDim]
    index0[:] = -1

    weight0 = numpy.nan * numpy.ones(outSizes)

    for outIndex in range(outInterpSize):
        outInd = [slice(None)] * xOut.ndim
        outInd[outAxis] = outIndex
        outInd = tuple(outInd)
        x = numpy.array(xOut[outInd])
        for inIndex in range(inInterpSize - 1):
            ind0 = [slice(None)] * xIn.ndim
            ind0[inAxis] = inIndex
            ind0 = tuple(ind0)
            x0 = numpy.array(xIn[ind0])
            ind1 = [slice(None)] * xIn.ndim
            ind1[inAxis] = inIndex + 1
            ind1 = tuple(ind1)
            x1 = numpy.array(xIn[ind1])
            dx = x1 - x0
            frac = (x - x0) / dx
            valid = numpy.isfinite(frac)
            mask = numpy.zeros(valid.shape, bool)
            mask[valid] = numpy.logical_and(frac[valid] >= 0.,
                                            frac[valid] < 1.)
            if inIndex == inInterpSize - 2:
                mask = numpy.logical_or(mask, x == x1)
            if numpy.count_nonzero(mask) == 0:
                continue

            localIndex = numpy.array(index0[outInd])
            localIndex[mask] = inIndex
            index0[outInd] = localIndex

            localWeight = numpy.array(weight0[outInd])
            localWeight[mask] = 1. - frac[mask]
            weight0[outInd] = localWeight

    for dim in indices:
        indices[dim] = xarray.DataArray(indices[dim], dims=allOutDims)
    weight0 = xarray.DataArray(weight0, dims=allOutDims)

    return indices, weight0


def _interp_1d_array(da, indices, weight0, inInterpDim):
    """
    iterpolate a data array in 1D with the given indices and weights
    """
    if inInterpDim not in da.dims:
        return da

    indices0 = {}
    indices1 = {}
    for dim in da.dims:
        if dim in indices.keys():
            indices0[dim] = indices[dim]
            if dim == inInterpDim:
                indices1[dim] = indices[dim] + 1
            else:
                indices1[dim] = indices[dim]

    da0 = da.isel(**indices0)
    da1 = da.isel(**indices1)
    da = weight0 * da0 + (1. - weight0) * da1
    return da
