#!/usr/bin/env python

'''
Make a mask file for a given mesh from geometric features describing
MOC regions.

The -m flag is used to specify the name of the ACME mesh to which the
masks should be applied.

Requires:
    * a local link to the MPAS mask creator MpasMaskCreator.x
    * a local link to the MOC southern boundary extractor tool
      moc_southern_boundary_extractor.py
    * a local link to a mesh file named <mesh_name>_mesh.nc describing the
      desired mesh
    * the region file MOCBasins.geojson produced by running
      ./driver_scripts/setup_ocean_region_groups.py in the geometric_features
      repo

Produces:
    * <mesh_name>_MOCBasinsAndTransectMasks.nc, the mask file

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
maskFileName = '{}_MOCBasinsAndTransectMasks.nc'.format(args.mesh_name)
regionFileName = 'MOCBasins.geojson'

tempRegionMaskFile = 'tempRegionMasks.nc'
subprocess.check_call(['./MpasMaskCreator.x', meshFileName,
                       tempRegionMaskFile, '-f', regionFileName])

subprocess.check_call(['./moc_southern_boundary_extractor.py',
                       '-f', tempRegionMaskFile,
                       '-m', meshFileName,
                       '-o', maskFileName])

subprocess.check_call(['rm', tempRegionMaskFile])


