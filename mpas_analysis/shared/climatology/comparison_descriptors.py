# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#
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

from mpas_analysis.shared.grid import LatLonGridDescriptor, \
    ProjectionGridDescriptor


def get_comparison_descriptor(config, comparisonGridName):  # {{{
    """
    Get the comparison grid descriptor from the comparisonGridName.

    Parameters
    ----------
    config :  MpasAnalysisConfigParser object
        Contains configuration options

    comparisonGridName : {'latlon', 'antarctic'}
        The name of the comparison grid to use for remapping.

    Raises
    ------
    ValueError
        If comparisonGridName does not describe a known comparions grid
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    if comparisonGridName == 'latlon':
        comparisonDescriptor = \
            _get_lat_lon_comparison_descriptor(config)
    elif comparisonGridName == 'antarctic':
        comparisonDescriptor = \
            _get_antarctic_stereographic_comparison_descriptor(config)
    else:
        raise ValueError('Unknown comaprison grid type {}'.format(
            comparisonGridName))
    return comparisonDescriptor  # }}}


def get_antarctic_stereographic_projection():  # {{{
    """
    Get a projection for an Antarctic steregraphic comparison grid

    Returns
    -------
    projection : ``pyproj.Proj`` object
        The projection
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    projection = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 '
                             '+k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')

    return projection  # }}}


def _get_lat_lon_comparison_descriptor(config):  # {{{
    """
    Get a descriptor of the lat/lon comparison grid, used for remapping and
    determining the grid name

    Parameters
    ----------
    config :  instance of ``MpasAnalysisConfigParser``
        Contains configuration options

    Returns
    -------
    descriptor : ``LatLonGridDescriptor`` object
        A descriptor of the lat/lon grid
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    climSection = 'climatology'

    comparisonLatRes = config.getWithDefault(climSection,
                                             'comparisonLatResolution',
                                             constants.dLatitude)
    comparisonLonRes = config.getWithDefault(climSection,
                                             'comparisonLatResolution',
                                             constants.dLongitude)

    nLat = int((constants.latmax-constants.latmin)/comparisonLatRes)+1
    nLon = int((constants.lonmax-constants.lonmin)/comparisonLonRes)+1
    lat = numpy.linspace(constants.latmin, constants.latmax, nLat)
    lon = numpy.linspace(constants.lonmin, constants.lonmax, nLon)

    descriptor = LatLonGridDescriptor.create(lat, lon, units='degrees')

    return descriptor  # }}}


def _get_antarctic_stereographic_comparison_descriptor(config):  # {{{
    """
    Get a descriptor of an Antarctic steregraphic comparison grid, used for
    remapping and determining the grid name

    Parameters
    ----------
    config :  instance of ``MpasAnalysisConfigParser``
        Contains configuration options

    Returns
    -------
    descriptor : ``ProjectionGridDescriptor`` object
        A descriptor of the Antarctic comparison grid
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    climSection = 'climatology'

    comparisonStereoWidth = config.getfloat(climSection,
                                            'comparisonAntarcticStereoWidth')
    comparisonStereoResolution = config.getfloat(
        climSection, 'comparisonAntarcticStereoResolution')

    projection = get_antarctic_stereographic_projection()

    xMax = 0.5*comparisonStereoWidth*1e3
    nx = int(comparisonStereoWidth/comparisonStereoResolution)+1
    x = numpy.linspace(-xMax, xMax, nx)

    meshName = '{}x{}km_{}km_Antarctic_stereo'.format(
        comparisonStereoWidth, comparisonStereoWidth,
        comparisonStereoResolution)
    descriptor = ProjectionGridDescriptor.create(projection, x, x, meshName)

    return descriptor  # }}}
