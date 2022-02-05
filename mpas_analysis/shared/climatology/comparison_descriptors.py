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
"""
Functions for creating climatologies from monthly time series data
"""
# Authors
# -------
# Xylar Asay-Davis

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
import pyproj

from mpas_analysis.shared.constants import constants

from pyremap import LatLonGridDescriptor, ProjectionGridDescriptor


known_comparison_grids = ['latlon', 'antarctic', 'arctic', 'north_atlantic',
                          'north_pacific']


def get_comparison_descriptor(config, comparison_grid_name):  # {{{
    """
    Get the comparison grid descriptor from the comparison_grid_name.

    Parameters
    ----------
    config :  mpas_analysis.configuration.MpasAnalysisConfigParser
        Contains configuration options

    comparison_grid_name : {'latlon', 'antarctic', 'arctic', 'north_atlantic',
                            'north_pacific'}
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

    return comparison_descriptor  # }}}


def get_projection(comparison_grid_name):  # {{{
    """
    Get the projection from the comparison_grid_name.

    Parameters
    ----------
    comparison_grid_name : {'antarctic', 'arctic', 'north_atlantic',
                            'north_pacific'}
        The name of the comparison grid to use for remapping.

    Raises
    ------
    ValueError
        If comparison_grid_name does not describe a known comparison grid
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    if comparison_grid_name not in known_comparison_grids:
        raise ValueError(
            f'Unknown comparison grid type {comparison_grid_name}')

    if comparison_grid_name == 'latlon':
        raise ValueError('latlon is not a projection grid.')
    elif comparison_grid_name == 'antarctic':
        projection = get_antarctic_stereographic_projection()
    elif comparison_grid_name == 'arctic':
        projection = get_arctic_stereographic_projection()
    elif comparison_grid_name == 'north_atlantic':
        projection = get_north_atlantic_projection()
    elif comparison_grid_name == 'north_atlantic':
        projection = get_north_atlantic_projection()
    else:
        raise ValueError(f'We missed one of the known comparison grids: '
                         f'{comparison_grid_name}')

    return projection  # }}}


def get_antarctic_stereographic_projection():  # {{{
    """
    Get a projection for an Antarctic stereographic comparison grid

    Returns
    -------
    projection : pyproj.Proj
        The projection
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    projection = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 '
                             '+k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')

    return projection  # }}}


def get_arctic_stereographic_projection():  # {{{
    """
    Get a projection for an Arctic stereographic comparison grid

    Returns
    -------
    projection : pyproj.Proj
        The projection
    """
    # Authors
    # -------
    # Milena Veneziani

    projection = pyproj.Proj('+proj=stere +lat_ts=75.0 +lat_0=90 +lon_0=0.0 '
                             '+k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')

    return projection  # }}}


def get_north_atlantic_projection():  # {{{
    """
    Get a projection for a North Atlantic Lambert Conformal Conic (LCC)
    comparison grid

    Returns
    -------
    projection : pyproj.Proj
        The projection
    """
    # Authors
    # -------
    # Xylar Asay-Davis
    # Yohei Takano

    projection = pyproj.Proj('+proj=lcc +lon_0=-45 +lat_1=20 +lat_2=75 '
                             '+x_0=0.0 +y_0=0.0 +ellps=WGS84')

    return projection  # }}}


def get_north_pacific_projection():  # {{{
    """
    Get a projection for a North Pacific Lambert Conformal Conic (LCC)
    comparison grid

    Returns
    -------
    projection : pyproj.Proj
        The projection
    """
    # Authors
    # -------
    # Xylar Asay-Davis
    # Yohei Takano

    projection = pyproj.Proj('+proj=lcc +lon_0=180 +lat_1=20 +lat_2=60 '
                             '+x_0=0.0 +y_0=0.0 +ellps=WGS84')

    return projection  # }}}


def _get_lat_lon_comparison_descriptor(config):  # {{{
    """
    Get a descriptor of the lat/lon comparison grid, used for remapping and
    determining the grid name

    Parameters
    ----------
    config :  mpas_analysis.configuration.MpasAnalysisConfigParser
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

    lat_res = config.getWithDefault(section, 'comparisonLatResolution',
                                    constants.dLatitude)
    lon_res = config.getWithDefault(section, 'comparisonLatResolution',
                                    constants.dLongitude)

    nlat = int((constants.latmax - constants.latmin) / lat_res) + 1
    nlon = int((constants.lonmax - constants.lonmin) / lon_res) + 1
    lat = numpy.linspace(constants.latmin, constants.latmax, nlat)
    lon = numpy.linspace(constants.lonmin, constants.lonmax, nlon)

    descriptor = LatLonGridDescriptor.create(lat, lon, units='degrees')

    return descriptor  # }}}


def _get_projection_comparison_descriptor(config, comparison_grid_name):  # {{{
    """
    Get a descriptor of any comparison grid base on a projection, used for
    remapping and determining the grid name

    Parameters
    ----------
    config :  mpas_analysis.configuration.MpasAnalysisConfigParser
        Contains configuration options

    comparison_grid_name : str
        One of the projections

    Returns
    -------
    descriptor : pyremap.ProjectionGridDescriptor
        A descriptor of the Arctic comparison grid
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    section = 'climatology'

    option_suffixes = {'antarctic': 'AntarcticStereo',
                       'arctic': 'ArcticStereo',
                       'north_atlantic': 'NorthAtlantic',
                       'north_pacific': 'NorthPacific'}

    grid_suffixes = {'antarctic': 'Antarctic_stereo',
                     'arctic': 'Arctic_stereo',
                     'north_atlantic': 'North_Atlantic',
                     'north_pacific': 'North_Pacific'}

    if comparison_grid_name not in option_suffixes:
        raise ValueError(f'{comparison_grid_name} is not one of the supported '
                         f'projection grids')

    projection = get_projection(comparison_grid_name)

    option_suffix = option_suffixes[comparison_grid_name]
    grid_suffix = grid_suffixes[comparison_grid_name]
    width = config.getfloat(
        section, f'comparison{option_suffix}Width')
    option = f'comparison{option_suffix}Height'
    if config.has_option(section, option):
        height = config.getfloat(section, option)
    else:
        height = width
    res = config.getfloat(
        section, f'comparison{option_suffix}Resolution')

    xmax = 0.5 * width * 1e3
    nx = int(width / res) + 1
    x = numpy.linspace(-xmax, xmax, nx)

    ymax = 0.5 * height * 1e3
    ny = int(height / res) + 1
    y = numpy.linspace(-ymax, ymax, ny)

    mesh_name = f'{width}x{height}km_{res}km_{grid_suffix}'
    descriptor = ProjectionGridDescriptor.create(projection, x, y, mesh_name)

    return descriptor  # }}}
