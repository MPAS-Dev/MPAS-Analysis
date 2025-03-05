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

import pyproj
import cartopy


known_comparison_grids = ['latlon', 'antarctic', 'antarctic_extended',
                          'arctic', 'arctic_extended', 'north_atlantic',
                          'north_pacific', 'subpolar_north_atlantic', 'fris']


def get_pyproj_projection(comparison_grid_name):
    """
    Get the projection from the comparison_grid_name.

    Parameters
    ----------
    comparison_grid_name : str
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
    elif comparison_grid_name in ['antarctic', 'antarctic_extended', 'fris']:
        projection = pyproj.Proj(
            '+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0  +k_0=1.0 '
            '+x_0=0.0 +y_0=0.0 +ellps=WGS84')
    elif comparison_grid_name in ['arctic', 'arctic_extended']:
        projection = pyproj.Proj(
            '+proj=stere +lat_ts=75.0 +lat_0=90 +lon_0=0.0  +k_0=1.0 '
            '+x_0=0.0 +y_0=0.0 +ellps=WGS84')
    elif comparison_grid_name == 'north_atlantic':
        projection = pyproj.Proj('+proj=lcc +lon_0=-45 +lat_0=45 +lat_1=39 '
                                 '+lat_2=51 +x_0=0.0 +y_0=0.0 +ellps=WGS84')
    elif comparison_grid_name == 'north_pacific':
        projection = pyproj.Proj('+proj=lcc +lon_0=180 +lat_0=40 +lat_1=34 '
                                 '+lat_2=46 +x_0=0.0 +y_0=0.0 +ellps=WGS84')
    elif comparison_grid_name == 'subpolar_north_atlantic':
        projection = pyproj.Proj('+proj=lcc +lon_0=-40 +lat_0=54 +lat_1=40 '
                                 '+lat_2=68 +x_0=0.0 +y_0=0.0 +ellps=WGS84')
    else:
        raise ValueError(f'We missed one of the known comparison grids: '
                         f'{comparison_grid_name}')

    return projection


def get_cartopy_projection(comparison_grid_name):
    """
    Get the projection from the comparison_grid_name.

    Parameters
    ----------
    comparison_grid_name : str
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

    elif comparison_grid_name in ['antarctic', 'antarctic_extended', 'fris']:
        projection = cartopy.crs.Stereographic(
            central_latitude=-90., central_longitude=0.0,
            true_scale_latitude=-71.0)
    elif comparison_grid_name in ['arctic', 'arctic_extended']:
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
    elif comparison_grid_name == 'subpolar_north_atlantic':
        projection = cartopy.crs.LambertConformal(
            central_latitude=54., central_longitude=-40.,
            standard_parallels=(40., 68.))
    else:
        raise ValueError(f'We missed one of the known comparison grids: '
                         f'{comparison_grid_name}')

    return projection
