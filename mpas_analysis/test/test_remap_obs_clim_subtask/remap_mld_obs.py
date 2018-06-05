# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#
import numpy

from mpas_analysis.shared.interpolation import Remapper
from mpas_analysis.shared.grid import LatLonGridDescriptor
from mpas_analysis.shared.constants import constants

inputFileName = '/media/xylar/extra_data/analysis/output/GMPAS-QU240/' \
    'remap_obs/clim/obs/mld_1.0x1.0degree.nc'

obsDescriptor = LatLonGridDescriptor.read(fileName=inputFileName,
                                          latVarName='lat',
                                          lonVarName='lon')

comparisonLatRes = 4.
comparisonLonRes = 4.

nLat = int((constants.latmax-constants.latmin)/comparisonLatRes)+1
nLon = int((constants.lonmax-constants.lonmin)/comparisonLonRes)+1
lat = numpy.linspace(constants.latmin, constants.latmax, nLat)
lon = numpy.linspace(constants.lonmin, constants.lonmax, nLon)

comparisonDescriptor = LatLonGridDescriptor.create(lat, lon, units='degrees')

remapper = Remapper(obsDescriptor, comparisonDescriptor,
                    mappingFileName='map.nc')

remapper.build_mapping_file()

remapper.remap_file(inputFileName, 'mld_4.0x4.0degree.nc',
                    ['mld', 'month', 'year'],
                    renormalize=0.05)
