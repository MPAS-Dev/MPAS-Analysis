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

import pyproj
import cartopy


known_comparison_grids = ['latlon', 'antarctic', 'arctic', 'north_atlantic',
                          'north_pacific']


def get_pyproj_projection(comparison_grid_name):
    """
    Get the projection from the comparison_grid_name.

    Parameters
    ----------
    comparison_grid_name : {'antarctic', 'arctic', 'north_atlantic',
                            'north_pacific'}
        The name of the projection comparison grid to use for remapping

    Returns
    -------
    projection : pyproj.Proj
        The projection

    Raises
    ------
    ValueError
        If comparison_grid_name does not describe a known comparison grid
    """
    # Authors
    # -------
    # Xylar Asay-Davis
    # Milena Veneziani
    # Yohei Takano

    if comparison_grid_name not in known_comparison_grids:
        raise ValueError(
            f'Unknown comparison grid type {comparison_grid_name}')

    if comparison_grid_name == 'latlon':
        raise ValueError('latlon is not a projection grid.')
    elif comparison_grid_name == 'antarctic':
        projection = pyproj.Proj(
            '+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0  +k_0=1.0 '
            '+x_0=0.0 +y_0=0.0 +ellps=WGS84')
    elif comparison_grid_name == 'arctic':
        projection = pyproj.Proj(
            '+proj=stere +lat_ts=75.0 +lat_0=90 +lon_0=0.0  +k_0=1.0 '
            '+x_0=0.0 +y_0=0.0 +ellps=WGS84')
    elif comparison_grid_name == 'north_atlantic':
        projection = pyproj.Proj('+proj=lcc +lon_0=-45 +lat_0=45 +lat_1=39 '
                                 '+lat_2=51 +x_0=0.0 +y_0=0.0 +ellps=WGS84')
    elif comparison_grid_name == 'north_pacific':
        projection = pyproj.Proj('+proj=lcc +lon_0=180 +lat_0=40 +lat_1=34 '
                                 '+lat_2=46 +x_0=0.0 +y_0=0.0 +ellps=WGS84')
    else:
        raise ValueError(f'We missed one of the known comparison grids: '
                         f'{comparison_grid_name}')

    return projection


def get_cartopy_projection(comparison_grid_name):
    """
    Get the projection from the comparison_grid_name.

    Parameters
    ----------
    comparison_grid_name : {'antarctic', 'arctic', 'north_atlantic',
                            'north_pacific'}
        The name of the projection comparison grid to use for remapping

    Returns
    -------
    projection : cartopy.crs.Projection
        The projection

    Raises
    ------
    ValueError
        If comparison_grid_name does not describe a known comparison grid
    """
    # Authors
    # -------
    # Xylar Asay-Davis
    # Milena Veneziani
    # Yohei Takano

    if comparison_grid_name not in known_comparison_grids:
        raise ValueError(
            f'Unknown comparison grid type {comparison_grid_name}')

    if comparison_grid_name == 'latlon':
        raise ValueError('latlon is not a projection grid.')

    elif comparison_grid_name == 'antarctic':
        projection = cartopy.crs.Stereographic(
            central_latitude=-90., central_longitude=0.0,
            true_scale_latitude=-71.0)
    elif comparison_grid_name == 'arctic':
        projection = cartopy.crs.Stereographic(
            central_latitude=90., central_longitude=0.0,
            true_scale_latitude=75.0)
    elif comparison_grid_name == 'north_atlantic':
        projection = cartopy.crs.LambertConformal(
            central_latitude=45., central_longitude=-45.,
            standard_parallels=(39., 51.))
    elif comparison_grid_name == 'north_pacific':
        projection = cartopy.crs.LambertConformal(
            central_latitude=40., central_longitude=180.,
            standard_parallels=(34., 46.))
    else:
        raise ValueError(f'We missed one of the known comparison grids: '
                         f'{comparison_grid_name}')

    return projection
