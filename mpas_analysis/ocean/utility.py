# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#
"""
A utility for computing common ocean fields (e.g. zMid) from datasets
"""
# Authors
# -------
# Xylar Asay-Davis

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
import xarray


def compute_zmid(bottomDepth, maxLevelCell, layerThickness):  # {{{
    """
    Computes zMid given data arrays for bottomDepth, maxLevelCell and
    layerThickness

    Parameters
    ----------
    bottomDepth : ``xarray.DataArray``
        the depth of the ocean bottom (positive)

    maxLevelCell : ``xarray.DataArray``
        the 1-based vertical index of the bottom of the ocean

    layerThickness : ``xarray.DataArray``
        the thickness of MPAS-Ocean layers (possibly as a function of time)

    Returns
    -------
    zMid : ``xarray.DataArray``
        the vertical coordinate defining the middle of each layer, masked below
        the bathymetry
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    nVertLevels = layerThickness.sizes['nVertLevels']

    vertIndex = \
        xarray.DataArray.from_dict({'dims': ('nVertLevels',),
                                    'data': numpy.arange(nVertLevels)})

    layerThickness = layerThickness.where(vertIndex < maxLevelCell)

    thicknessSum = layerThickness.sum(dim='nVertLevels')
    thicknessCumSum = layerThickness.cumsum(dim='nVertLevels')
    zSurface = -bottomDepth+thicknessSum

    zLayerBot = zSurface - thicknessCumSum

    zMid = zLayerBot + 0.5*layerThickness

    return zMid  # }}}


def nans_to_numpy_mask(field):  # {{{
    """
    Convert a numpy array with NaNs to a masked numpy array
    """
    field = numpy.ma.masked_array(
        field, numpy.isnan(field))
    return field  # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
