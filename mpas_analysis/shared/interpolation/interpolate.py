"""
    Module that provides basic nearest neighbor interpolation functionality

    Author: Luke Van Roekel
    Modified: 10/24/2016
"""

import numpy as np
from scipy.spatial import cKDTree
import sys

def lon_lat_to_cartesian(lon, lat, R = 6371222.):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R
    """

    if max(abs(lon)) < 3.0*np.pi:
        lon_r = lon
        lat_r = lat
    else:
        lon_r = np.radians(lon)
        lat_r = np.radians(lat)

    x = R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x,y,z

def init_tree(lon_input, lat_input, lonmin, lonmax, latmin, latmax, dLon, dLat):
    """
    Initializes a KD tree for nearest neighbor searching
    """

    lon_input = lon_input.flatten()
    lat_input = lat_input.flatten()

    if max(lon_input) < 2.*np.pi:
        lon_input = np.rad2deg(lon_input)
    if max(lat_input) < np.pi / 2.:
        lat_input = np.rad2deg(lat_input)

    if max(lon_input) > 180.:
        inds = np.where(lon_input > 180.)
        lon_input[inds] -= 360.

    if lonmax > 180.:
        sys.exit("longitude bounds must be between -180 and 180")

    xs, ys, zs = lon_lat_to_cartesian(lon_input,lat_input)
    tree = cKDTree(zip(xs, ys, zs))

    lonVals = np.arange(lonmin + dLon/2., lonmax + dLon/2., dLon)
    latVals = np.arange(latmin + dLat/2., latmax + dLat/2., dLat)

    latTarg, lonTarg = np.meshgrid(latVals,lonVals)
    xt, yt, zt = lon_lat_to_cartesian(lonTarg.flatten(),latTarg.flatten())

    d, inds = tree.query(zip(xt, yt, zt), k = 1)

    return d, inds, lonTarg, latTarg

def interp_fields(field, d, inds, lonTarg):
    """
        performs nearest neighbor interpolation
    """

    return field.flatten()[inds].reshape(lonTarg.shape)
