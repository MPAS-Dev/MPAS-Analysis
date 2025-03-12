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
"""
A utility for computing common ocean fields (e.g. zMid) from datasets
"""
# Authors
# -------
# Xylar Asay-Davis

import numpy
import xarray


def add_standard_regions_and_subset(ds, config, regionShortNames=None):
    """
    Add standard region names (``regionNames`` coordinate) to a dataset and
    rename ``nOceanRegionsTmp`` dimension to ``nOceanRegions`` (if present).
    Shorter standard region names are in ``regionNamesShort``.

    Parameters
    ----------
    ds : xarray.Dataset
        the dataset to which region names should be added

    config : mpas_tools.config.MpasConfigParser
        Configuration options

    regionShortNames : list of str, optional
        A list of a subset of the short region names to use to subset the
        dataset

    Returns
    -------
    ds : xarray.Dataset
        the dataset with region names added and possibly subsetted
    """
    ds = ds.copy()
    if 'nOceanRegionsTmp' in ds.dims:
        ds = ds.rename({'nOceanRegionsTmp': 'nOceanRegions'})

    allShortNames = config.getexpression('regions', 'regionShortNames')
    regionNames = config.getexpression('regions', 'regionNames')
    ds.coords['regionShortNames'] = ('nOceanRegions', allShortNames)
    ds.coords['regionNames'] = ('nOceanRegions', regionNames)
    if regionShortNames is not None:
        regionIndices = \
            [allShortNames.index(name) for name in regionShortNames]
        ds = ds.isel(nOceanRegions=regionIndices)
    return ds


def get_standard_region_names(config, regionShortNames):
    """
    Add standard region names from the short names

    Parameters
    ----------
    config : mpas_tools.config.MpasConfigParser
        Configuration options

    regionShortNames : list of str
        A list of short region names

    Returns
    -------
    regionNames : list of str
        A list of full standard region names
    """
    allShortNames = config.getexpression('regions', 'regionShortNames')
    regionNames = config.getexpression('regions', 'regionNames')
    regionNameMap = {shortName: regionName for shortName, regionName in
                     zip(allShortNames, regionNames)}
    regionNames = [regionNameMap[shortName] for shortName in regionShortNames]

    return regionNames


def compute_zmid(bottomDepth, maxLevelCell, layerThickness):
    """
    Computes zMid given data arrays for bottomDepth, maxLevelCell and
    layerThickness

    Parameters
    ----------
    bottomDepth : ``xarray.DataArray``
        the depth of the ocean bottom (positive)

    maxLevelCell : ``xarray.DataArray``
        the 0-based vertical index of the bottom of the ocean

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

    layerThickness = layerThickness.where(vertIndex <= maxLevelCell)

    thicknessSum = layerThickness.sum(dim='nVertLevels')
    thicknessCumSum = layerThickness.cumsum(dim='nVertLevels')
    zSurface = -bottomDepth + thicknessSum

    zLayerBot = zSurface - thicknessCumSum

    zMid = zLayerBot + 0.5 * layerThickness

    return zMid


def compute_zinterface(bottomDepth, maxLevelCell, layerThickness):
    """
    Computes zInterface given data arrays for bottomDepth, maxLevelCell and
    layerThickness

    Parameters
    ----------
    bottomDepth : ``xarray.DataArray``
        the depth of the ocean bottom (positive)

    maxLevelCell : ``xarray.DataArray``
        the 0-based vertical index of the bottom of the ocean

    layerThickness : ``xarray.DataArray``
        the thickness of MPAS-Ocean layers (possibly as a function of time)

    Returns
    -------
    zInterface : ``xarray.DataArray``
        the vertical coordinate defining the interfaces between layers, masked
        below the bathymetry
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    nVertLevels = layerThickness.sizes['nVertLevels']

    vertIndex = \
        xarray.DataArray.from_dict({'dims': ('nVertLevels',),
                                    'data': numpy.arange(nVertLevels)})

    layerThickness = layerThickness.where(vertIndex <= maxLevelCell)
    thicknessSum = layerThickness.sum(dim='nVertLevels')

    zSurface = -bottomDepth + thicknessSum

    zInterfaceList = [zSurface]

    zTop = zSurface

    for zIndex in range(nVertLevels):
        zBot = zTop - layerThickness.isel(nVertLevels=zIndex)
        zInterfaceList.append(zBot)
        zTop = zBot

    zInterface = xarray.concat(zInterfaceList, dim='nVertLevelsP1').transpose(
        'nCells', 'nVertLevelsP1')
    return zInterface
