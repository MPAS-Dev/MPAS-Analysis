#!/usr/bin/env python
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

"""
Creates a region mask file from ISMIP6 regions on the ISMIP6 8km polar
stereographic grid.  Both a geojson file for the regions and the mask file
will be stored in the given output directory.  Because mask computation with
shapely is relatively slow, the computation can be sped up by running several
threads in parallel.

Usage: Symlink the mpas_analysis directory from one directory up.
Modify the mesh and regions names and the local path to region masks and mesh
file. Optionally, change the number of threads.
"""

import xarray
import os

from geometric_features import GeometricFeatures
from geometric_features.aggregation import get_aggregator_by_name

from mpas_tools.parallel import create_pool
from mpas_tools.io import write_netcdf
from mpas_tools.logging import LoggingContext
from mpas_tools.mesh.mask import compute_projection_grid_region_masks


if __name__ == '__main__':
    meshName = 'ISMIP6_8km'

    ismip6_aggregation, prefix, date = get_aggregator_by_name('ISMIP6 Regions')

    regionsName = '{}{}'.format(prefix, date)

    # the number of parallel threads to use
    processCount = 8

    regionMasksPath = '/home/xylar/Desktop/region_masks'

    # replace with the path to the desired mesh or restart file
    meshFileName = '/home/xylar/Desktop/obs_temperature_1995-2017_8km_x_20m.nc'

    try:
        os.makedirs(regionMasksPath)
    except OSError:
        pass

    geojsonFileName = '{}/{}.geojson'.format(regionMasksPath, regionsName)
    maskFileName = '{}/{}_{}.nc'.format(regionMasksPath, meshName, regionsName)

    pool = create_pool(process_count=processCount, method='spawn')

    gf = GeometricFeatures()
    fc = ismip6_aggregation(gf)
    fc.to_geojson(geojsonFileName)

    with xarray.open_dataset(meshFileName) as ds:
        lon = ds.lon
        lat = ds.lat

        ydim, xdim = lon.dims

        with LoggingContext('compute_lon_lat_region_masks') as logger:
            dsMasks = compute_projection_grid_region_masks(
                lon=lon.values, lat=lat.values, fcMask=fc, logger=logger,
                pool=pool, showProgress=True, xdim=xdim, ydim=ydim)

        write_netcdf(dsMasks, maskFileName)
