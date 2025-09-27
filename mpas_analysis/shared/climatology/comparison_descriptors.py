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
Functions for creating climatologies from monthly time series data
"""
# Authors
# -------
# Xylar Asay-Davis

import numpy

from mpas_analysis.shared.constants import constants
from mpas_analysis.shared.projection import (
    comparison_grid_file_suffixes,
    comparison_grid_option_suffixes,
    get_pyproj_projection,
    known_comparison_grids,
)

from pyremap import LatLonGridDescriptor, ProjectionGridDescriptor


def get_comparison_descriptor(config, comparison_grid_name):
    """
    Get the comparison grid descriptor from the comparison_grid_name.

    Parameters
    ----------
    config : tranche.Tranche
        Contains configuration options

    comparison_grid_name : {'latlon', 'antarctic', 'arctic', 'north_atlantic',
                            'north_pacific', 'subpolar_north_atlantic', 'fris'}
        The name of the comparison grid to use for remapping.

    Raises
    ------
    ValueError
        If comparison_grid_name does not describe a known comparison grid
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani

    if comparison_grid_name not in known_comparison_grids:
        raise ValueError(
            f'Unknown comparison grid type {comparison_grid_name}')

    if comparison_grid_name == 'latlon':
        comparison_descriptor = \
            _get_lat_lon_comparison_descriptor(config)
    else:
        comparison_descriptor = \
            _get_projection_comparison_descriptor(config, comparison_grid_name)

    return comparison_descriptor


def _get_lat_lon_comparison_descriptor(config):
    """
    Get a descriptor of the lat/lon comparison grid, used for remapping and
    determining the grid name

    Parameters
    ----------
    config : tranche.Tranche
        Contains configuration options

    Returns
    -------
    descriptor : LatLonGridDescriptor
        A descriptor of the lat/lon grid
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    section = 'climatology'

    lat_res = config.getfloat(section, 'comparisonLatResolution')
    lon_res = config.getfloat(section, 'comparisonLatResolution')

    nlat = int((constants.latmax - constants.latmin) / lat_res) + 1
    nlon = int((constants.lonmax - constants.lonmin) / lon_res) + 1
    lat = numpy.linspace(constants.latmin, constants.latmax, nlat)
    lon = numpy.linspace(constants.lonmin, constants.lonmax, nlon)

    descriptor = LatLonGridDescriptor.create(lat, lon, units='degrees')

    return descriptor


def _get_projection_comparison_descriptor(config, comparison_grid_name):
    """
    Get a descriptor of any comparison grid base on a projection, used for
    remapping and determining the grid name

    Parameters
    ----------
    config : tranche.Tranche
        Contains configuration options

    comparison_grid_name : str
        One of the projections

    Returns
    -------
    descriptor : pyremap.ProjectionGridDescriptor
        A descriptor of the comparison grid
        (eg. - Arctic, North Atlantic)
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    section = 'climatology'

    option_suffixes = comparison_grid_option_suffixes

    grid_suffixes = comparison_grid_file_suffixes

    if comparison_grid_name not in option_suffixes:
        raise ValueError(f'{comparison_grid_name} is not one of the supported '
                         f'projection grids')

    projection = get_pyproj_projection(comparison_grid_name)

    option_suffix = option_suffixes[comparison_grid_name]
    grid_suffix = grid_suffixes[comparison_grid_name]
    option = f'comparison{option_suffix}Bounds'
    if config.has_option(section, option):
        bounds = config.getexpression(section, option)
        bounds = [1e3 * bound for bound in bounds]
    else:
        width = config.getfloat(
            section, f'comparison{option_suffix}Width')
        option = f'comparison{option_suffix}Height'

        if config.has_option(section, option):
            height = config.getfloat(section, option)
        else:
            height = width
        xmax = 0.5 * width * 1e3
        ymax = 0.5 * height * 1e3
        bounds = [-xmax, xmax, -ymax, ymax]
    width = (bounds[1] - bounds[0]) / 1e3
    height = (bounds[3] - bounds[2]) / 1e3
    res = config.getfloat(
        section, f'comparison{option_suffix}Resolution')

    nx = int(width / res) + 1
    x = numpy.linspace(bounds[0], bounds[1], nx)

    ny = int(height / res) + 1
    y = numpy.linspace(bounds[2], bounds[3], ny)

    mesh_name = f'{width}x{height}km_{res}km_{grid_suffix}'
    descriptor = ProjectionGridDescriptor.create(projection, x, y, mesh_name)

    return descriptor
