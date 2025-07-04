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

import numpy as np
import xarray as xr


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

    vertIndex = xr.DataArray.from_dict(
        {
            'dims': ('nVertLevels',),
            'data': np.arange(nVertLevels)
        }
    )

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

    vertIndex = xr.DataArray.from_dict(
        {
            'dims': ('nVertLevels',),
            'data': np.arange(nVertLevels)
        }
    )

    layerThickness = layerThickness.where(vertIndex <= maxLevelCell)
    thicknessSum = layerThickness.sum(dim='nVertLevels')

    zSurface = -bottomDepth + thicknessSum

    zInterfaceList = [zSurface]

    zTop = zSurface

    for zIndex in range(nVertLevels):
        zBot = zTop - layerThickness.isel(nVertLevels=zIndex)
        zInterfaceList.append(zBot)
        zTop = zBot

    zInterface = xr.concat(zInterfaceList, dim='nVertLevelsP1').transpose(
        'nCells', 'nVertLevelsP1')
    return zInterface


def vector_cell_to_edge_isotropic(ds_mesh, zonal_cell, meridional_cell):
    """
    Compute the zonal and meridional components of a vector at edges from
    cell-centered components using isotropic area-weighted averaging.

    Parameters
    ----------
    ds_mesh : xarray.Dataset
        MPAS mesh variables, must include:
        - verticesOnEdge
        - cellsOnVertex
        - kiteAreasOnVertex

    zonal_cell : xarray.DataArray
        Zonal component at cell centers (nCells,)

    meridional_cell : xarray.DataArray
        Meridional component at cell centers (nCells,)

    Returns
    -------
    zonal_edge : xarray.DataArray
        Zonal component at edges (nEdges,)

    meridional_edge : xarray.DataArray
        Meridional component at edges (nEdges,)
    """
    vertices_on_edge = ds_mesh.verticesOnEdge - 1
    cells_on_vertex = ds_mesh.cellsOnVertex - 1
    kite_areas_on_vertex = ds_mesh.kiteAreasOnVertex

    n_edges = vertices_on_edge.sizes['nEdges']
    vertex_degree = cells_on_vertex.sizes['vertexDegree']

    zonal_edge = np.zeros(n_edges, dtype=float)
    meridional_edge = np.zeros(n_edges, dtype=float)
    area_sum = np.zeros(n_edges, dtype=float)

    for v in range(2):
        # all valid edges have 2 valid vertices on that edge
        voe = vertices_on_edge.isel(TWO=v)
        for c in range(vertex_degree):
            # cells on vertices on edge
            covoe = cells_on_vertex.isel(
                vertexDegree=c,
                nVertices=voe
            )
            valid = covoe >= 0
            valid_covoe = covoe.isel(nEdges=valid)
            valid_voe = voe.isel(nEdges=valid)
            area = kite_areas_on_vertex.isel(
                vertexDegree=c,
                nVertices=valid_voe
            ).values
            if np.any(area == 0):
                raise ValueError(
                    "Some kite areas of valid cells on vertex have zero area. "
                    "This seems to be a bug in the mesh or "
                    "vector_cell_to_edge_isotropic()."
                )
            zcell = zonal_cell.isel(nCells=valid_covoe).values
            mcell = meridional_cell.isel(nCells=valid_covoe).values
            zonal_edge[valid] += zcell * area
            meridional_edge[valid] += mcell * area
            area_sum[valid] += area

    if np.any(area_sum == 0):
        raise ValueError(
            "Some edges have zero area.  This seems to be a bug in the mesh "
            "or vector_cell_to_edge_isotropic()."
        )

    # Normalize by the area sum to get the average
    zonal_edge /= area_sum
    meridional_edge /= area_sum

    # Wrap as xarray DataArrays
    zonal_edge = xr.DataArray(zonal_edge, dims=('nEdges',))
    meridional_edge = xr.DataArray(meridional_edge, dims=('nEdges',))
    return zonal_edge, meridional_edge


def vector_to_edge_normal(ds_mesh, zonal_edge, meridional_edge):
    """
    Compute the normal component of a vector at an edge from
    the zonal and meridional components.

    Parameters
    ----------
    ds_mesh : xarray.Dataset
        MPAS mesh variables, must include:
        - angleEdge

    zonal_edge : xarray.DataArray
        Zonal component at edges (nEdges,)

    meridional_edge : xarray.DataArray
        Meridional component at edges (nEdges,)

    Returns
    -------
    normal_edge : xarray.DataArray
        Normal component at edges (nEdges,)
    """

    angle_edge = ds_mesh.angleEdge
    normal_edge = (
        np.cos(angle_edge) * zonal_edge + np.sin(angle_edge) * meridional_edge
    )
    return normal_edge
