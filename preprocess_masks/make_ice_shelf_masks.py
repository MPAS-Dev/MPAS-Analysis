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
Make a mask file for a given mesh from ice-shelf geometric features.

The -m flag is used to specify the name of the ACME mesh to which the
masks should be applied.

Requires:
    * a local link to the MPAS mask creator MpasMaskCreator.x
    * a local link to a mesh file named <mesh_name>_mesh.nc describing the
      desired mesh
    * the region file iceShelves.geojson produced by running
      ./driver_scripts/setup_ice_shelves.py in the geometric_features repo

Produces:
    * <mesh_name>_iceShelfMasks.nc, the mask file

Author: Xylar Asay-Davis
'''

import subprocess
import argparse


parser = \
    argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-m', '--mesh_name', dest='mesh_name',
                    help='The ACME name of the mesh', metavar='MESH_NAME',
                    required=True)
args = parser.parse_args()

meshFileName = '{}_mesh.nc'.format(args.mesh_name)
maskFileName = '{}_iceShelfMasks.nc'.format(args.mesh_name)
regionFileName = 'iceShelves.geojson'

tempRegionMaskFile = 'tempRegionMasks.nc'
subprocess.check_call(['./MpasMaskCreator.x', meshFileName, maskFileName,
                       '-f', regionFileName])
