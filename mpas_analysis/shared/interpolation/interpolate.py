import numpy as np
from scipy.spatial import cKDTree

def lon_lat_to_cartesian(lon, lat, R = 6371222):
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
            
def init_tree(lon_input, lat_input, dLon, dLat):                    
    lon_input = lon_input.flatten()
    lat_input = lat_input.flatten()

    if max(lon_input) > 180:
        inds = np.where(lon_input > 180)
        lon_input[inds] -= 360

    xs, ys, zs = lon_lat_to_cartesian(lon_input,lat_input)
    tree = cKDTree(zip(xs, ys, zs))

    lonVals = np.arange(-180. + dLon/2., 181 - dLon/2., dLon)
    latVals = np.arange(-90. + dLat/2., 91 - dLat/2., dLat)
    
    latTarg, lonTarg = np.meshgrid(latVals,lonVals)
    xt, yt, zt = lon_lat_to_cartesian(lonTarg.flatten(),latTarg.flatten())

    d, inds = tree.query(zip(xt, yt, zt), k = 1)

    return d, inds, lonTarg, latTarg

def interp_fields(field, d, inds, lonTarg):
    """
        performs nearest neighbor interpolation
    """

    return field.flatten()[inds].reshape(lonTarg.shape)
