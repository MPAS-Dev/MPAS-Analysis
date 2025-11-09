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
import numpy
from pyremap import LatLonGridDescriptor, Remapper

from mpas_analysis.shared.constants import constants

inputFileName = '/media/xylar/extra_data/analysis/output/GMPAS-QU240/' \
    'remap_obs/clim/obs/mld_1.0x1.0degree.nc'

obsDescriptor = LatLonGridDescriptor.read(filename=inputFileName,
                                          lat_var_name='lat',
                                          lon_var_name='lon')

comparisonLatRes = 4.
comparisonLonRes = 4.

nLat = int((constants.latmax - constants.latmin) / comparisonLatRes) + 1
nLon = int((constants.lonmax - constants.lonmin) / comparisonLonRes) + 1
lat = numpy.linspace(constants.latmin, constants.latmax, nLat)
lon = numpy.linspace(constants.lonmin, constants.lonmax, nLon)

comparisonDescriptor = LatLonGridDescriptor.create(lat, lon, units='degrees')

remapper = Remapper(
    map_filename='map.nc',
    src_descriptor=obsDescriptor,
    dst_descriptor=comparisonDescriptor,
)

remapper.build_map()

remapper.ncremap(
    inputFileName,
    'mld_4.0x4.0degree.nc',
    ['mld', 'month', 'year'],
    renormalize=0.05
)
