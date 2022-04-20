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
Functions for plotting inset maps in plots (e.g. for transects)
"""
# Authors
# -------
# Xylar Asay-Davis

import matplotlib.path
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
import numpy
import shapely.geometry

from geometric_features.plot import subdivide_geom


def add_inset(fig, fc, latlonbuffer=45., polarbuffer=5., width=1.0,
              height=1.0, lowerleft=None, xbuffer=None, ybuffer=None,
              maxlength=1.):
    """
    Plots an inset map showing the location of a transect or polygon.  Shapes
    are plotted on a polar grid if they are entirely poleward of +/-50 deg.
    latitude and with a lat/lon grid if not.

    Parameters
    ----------
    fig : ``matplotlib.figure.Figure``
        A matplotlib figure to add the inset to

    fc : ``geometric_features.FeatureCollection``
        A collection of regions, transects and/or points to plot in the inset

    latlonbuffer : float, optional
        The number of degrees lat/lon to use as a buffer around the shape(s)
        to plot if a lat/lon plot is used.

    polarbuffer : float, optional
        The number of degrees latitude to use as a buffer equatorward of the
        shape(s) in polar plots

    width, height : float, optional
        width and height in inches of the inset

    lowerleft : pair of floats, optional
        the location of the lower left corner of the axis in inches, default
        puts the inset in the upper right corner of ``fig``.

    xbuffer, ybuffer : float, optional
        right and top buffers from the top-right corner (in inches) if
        lowerleft is ``None``.

    maxlength : float or ``None``, optional
        Any segments longer than maxlength will be subdivided in the plot to
        ensure curvature.  If ``None``, no subdivision is performed.

    Returns
    -------
    inset : ``matplotlib.axes.Axes``
        The new inset axis
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    minLon, minLat, maxLon, maxLat = _get_bounds(fc)

    figsize = fig.get_size_inches()
    width /= figsize[0]
    height /= figsize[1]
    if lowerleft is None:
        if xbuffer is None:
            xbuffer = 0.1*width
        else:
            xbuffer /= figsize[0]
        if ybuffer is None:
            ybuffer = xbuffer*figsize[0]/figsize[1]
        else:
            ybuffer /= figsize[1]
        lowerleft = [1.0 - width - xbuffer, 1.0 - height - ybuffer]
    else:
        lowerleft = [lowerleft[0]/figsize[0], lowerleft[1]/figsize[1]]
    bounds = [lowerleft[0], lowerleft[1], width, height]

    if maxLat <= -50:
        # an Antarctic-focused map makes the most sense
        inset = fig.add_axes(bounds,
                             projection=ccrs.SouthPolarStereo())
        extent = [-180., 180., -90., max(-65., maxLat+polarbuffer)]
        _set_circular_boundary(inset)
        xlocator = mticker.FixedLocator(numpy.linspace(-180., 180., 9))
        ylocator = mticker.FixedLocator(numpy.linspace(-90., -50., 9))
    elif minLat >= 50:
        # an Arctic-focused map makes the most sense
        inset = fig.add_axes(bounds,
                             projection=ccrs.NorthPolarStereo())
        extent = [-180, 180, min(65., minLat-polarbuffer), 90]
        _set_circular_boundary(inset)
        xlocator = mticker.FixedLocator(numpy.linspace(-180., 180., 9))
        ylocator = mticker.FixedLocator(numpy.linspace(50., 90., 9))
    else:
        inset = fig.add_axes(bounds,
                             projection=ccrs.PlateCarree())
        extent = [max(-180., minLon-latlonbuffer),
                  min(180., maxLon+latlonbuffer),
                  max(-90., minLat-latlonbuffer),
                  min(90., maxLat+latlonbuffer)]
        xlocator = None
        ylocator = None

    # kind of like "top" justified -- graphics are toward the "north" end of
    # the subplot
    inset.set_anchor('N')

    inset.set_extent(extent,  ccrs.PlateCarree())
    inset.add_feature(cartopy.feature.LAND, zorder=1)
    inset.add_feature(cartopy.feature.OCEAN, zorder=0)

    gl = inset.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                         linewidth=0.5, color='gray', alpha=0.5,
                         linestyle='--')

    if xlocator is not None:
        gl.xlocator = xlocator

    if ylocator is not None:
        gl.ylocator = ylocator

    for feature in fc.features:
        geomtype = feature['geometry']['type']
        shape = shapely.geometry.shape(feature['geometry'])
        if maxlength is not None:
            shape = subdivide_geom(shape, shape.geom_type, maxlength)
        if geomtype in ['Polygon', 'MultiPolygon']:
            inset.add_geometries((shape,), crs=ccrs.PlateCarree(),
                                 edgecolor='blue', facecolor='blue', alpha=0.4,
                                 linewidth=1.)
        elif geomtype in ['Point', 'MultiPoint']:
            inset.add_geometries((shape,), crs=ccrs.PlateCarree(),
                                 edgecolor='none', facecolor='none', alpha=1.,
                                 markersize=3., markeredgecolor='k',
                                 markerfacecolor='k')
        else:
            inset.add_geometries((shape,), crs=ccrs.PlateCarree(),
                                 edgecolor='k', facecolor='none', alpha=1.,
                                 linewidth=1.)
            # put a red point at the beginning and a blue point at the end
            # of the transect to help show the orientation
            begin = shape.coords[0]
            end = shape.coords[-1]
            inset.plot(begin[0], begin[1], color='r', marker='o',
                       markersize=3., transform=ccrs.PlateCarree())
            inset.plot(end[0], end[1], color='g', marker='o',
                       markersize=3., transform=ccrs.PlateCarree())

    return inset


def _set_circular_boundary(ax):
    """Set the boundary of the given axis to be circular (for a polar plot)"""

    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    theta = numpy.linspace(0, 2*numpy.pi, 100)
    center = numpy.array([0.5, 0.5])
    radius = 0.5
    verts = numpy.vstack([numpy.sin(theta), numpy.cos(theta)]).T
    circle = matplotlib.path.Path(verts * radius + center)

    ax.set_boundary(circle, transform=ax.transAxes)


def _get_bounds(fc):
    """Compute the lon/lat bounding box for all transects and regions"""

    bounds = shapely.geometry.GeometryCollection()
    for feature in fc.features:
        shape = shapely.geometry.shape(feature['geometry'])
        shape_bounds = shapely.geometry.box(*shape.bounds)
        bounds = shapely.geometry.box(*bounds.union(shape_bounds).bounds)

    return bounds.bounds
