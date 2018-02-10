#!/usr/bin/env python
# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#

'''
Creates a mapping file that can be used with ncremap (NCO) to remap MPAS files
to a latitude/longitude grid.

Usage: Copy this script into the main MPAS-Analysis directory (up one level).
Modify the grid name, the path to the MPAS grid file and the output grid
resolution.
'''

from mpas_analysis.shared.interpolation import Remapper
from mpas_analysis.shared.grid import MpasMeshDescriptor
from mpas_analysis.shared.climatology import get_comparison_descriptor
from mpas_analysis.configuration import MpasAnalysisConfigParser


# replace with the MPAS mesh name
inGridName = 'oQU240'

# replace with the path to the desired mesh or restart file
inGridFileName = '/media/xylar/extra_data/analysis/edison/G-QU240-master-intel/run/mpaso.rst.0001-01-06_00000.nc'

config = MpasAnalysisConfigParser()
config.read('mpas_analysis/config.default')
# replace 1.0 with the desired resolution of the output mesh
config.set('climatology', 'comparisonLatResolution', '1.0')
config.set('climatology', 'comparisonLonResolution', '1.0')

inDescriptor = MpasMeshDescriptor(inGridFileName, inGridName)

outDescriptor = get_comparison_descriptor(config, 'latlon')
outGridName = outDescriptor.meshName

mappingFileName = 'map_{}_to_{}.nc'.format(inGridName, outGridName)

remapper = Remapper(inDescriptor, outDescriptor, mappingFileName)

remapper.build_mapping_file(method='bilinear')
