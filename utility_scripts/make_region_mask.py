#!/usr/bin/env python
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

'''
Creates a region mask file from a geojson file defining the regions and a mesh
file.  The geojson file should be in the diagnostics/mpas_analysis/region_masks
directory and the mask file will be stored in the same directory, saving the
time of computing the mask for each analysis run.  Because mask computation
with shapely is relatively slow, the computation can be sped up by running
several threads in parallel.

Usage: Copy this script into the main MPAS-Analysis directory (up one level).
Modify the mesh and regions names and the local path to region masks.
Optionally, change the number of threads.
'''

import logging
import sys

from mpas_analysis.shared.regions.compute_region_masks_subtask import \
    compute_mpas_region_masks

# replace with the MPAS mesh name
meshName = 'oEC60to30v3'

# the prefix/suffix used for this set of regions.  This should be the name
# of the geojson file (without extension) and will be used as a suffix on the
# resulting region masks file.
regionsName = 'oceanBasins'

# the number of parallel threads to use
processCount = 64

# replace with the path to the region masks (where the appropriate geojson file
# is stored and where the mask file should be written)
regionMasksPath = '/global/project/projectdirs/acme/xylar/' \
                  'diagnostics_not_public/mpas_analysis/region_masks'

# replace with the path to the desired mesh or restart file
meshFileName = '/global/project/projectdirs/acme/inputdata/ocn/mpas-o/' \
               'oEC60to30v3/oEC60to30v3_60layer.170905.nc'

geojsonFileName = '{}/{}.geojson'.format(regionMasksPath, regionsName)
maskFileName = '{}/{}_{}.nc'.format(regionMasksPath, meshName, regionsName)

logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

compute_mpas_region_masks(geojsonFileName, meshFileName, maskFileName,
                          logger=logger, processCount=processCount)
