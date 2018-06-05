#!/usr/bin/env python
# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#

"""
A script for downloading and preprocessing data sets from Schmidtko et al.
(2014, DOI: 10.1126/science.1256117) for use in MPAS-Analysis
"""
# Authors
# -------
# Xylar Asay-Davis

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
import xarray
import os
import sys
import argparse
import pandas
import gsw

from mpas_analysis.shared.io.download import download_files

from mpas_analysis.shared.interpolation import Remapper
from mpas_analysis.shared.grid import LatLonGridDescriptor
from mpas_analysis.shared.climatology.comparison_descriptors \
    import get_comparison_descriptor
from mpas_analysis.configuration \
    import MpasAnalysisConfigParser
from mpas_analysis.shared.io import write_netcdf


def text_to_netcdf(inDir, outDir):
    inFileName = '{}/Antarctic_shelf_data.txt'.format(inDir)
    outFileName = '{}/Schmidtko_et_al_2014_bottom_PT_S_PD_' \
                  'SouthernOcean_0.25x0.125degree.nc'.format(outDir)

    if os.path.exists(outFileName):
        return

    # 1/4 x 1/8 degree grid cells
    cellsPerLon = 4
    cellsPerLat = 8

    obsFile = pandas.read_csv(inFileName, delim_whitespace=True)

    inLon = numpy.array(obsFile.iloc[:, 0])
    inLat = numpy.array(obsFile.iloc[:, 1])

    inZ = numpy.array(obsFile.iloc[:, 2])

    inCT = numpy.array(obsFile.iloc[:, 3])
    inCT_std = numpy.array(obsFile.iloc[:, 4])

    inSA = numpy.array(obsFile.iloc[:, 5])
    inSA_std = numpy.array(obsFile.iloc[:, 6])

    pressure = gsw.p_from_z(inZ, inLat)
    inS = gsw.SP_from_SA(inSA, pressure, inLon, inLat)
    inPT = gsw.pt_from_CT(inSA, inCT)
    inPD = gsw.rho(inSA, inCT, 0.)

    minLat = int(numpy.amin(inLat)*cellsPerLat)/cellsPerLat
    maxLat = int(numpy.amax(inLat)*cellsPerLat)/cellsPerLat
    deltaLat = 1./cellsPerLat
    outLat = numpy.arange(minLat-deltaLat, maxLat+2*deltaLat, deltaLat)

    deltaLon = 1./cellsPerLon
    outLon = numpy.arange(0., 360., deltaLon)

    xIndices = numpy.array(cellsPerLon*inLon + 0.5, int)
    yIndices = numpy.array(cellsPerLat*(inLat - outLat[0]) + 0.5, int)

    Lon, Lat = numpy.meshgrid(outLon, outLat)

    ds = xarray.Dataset()
    ds['lon'] = (('lon',), outLon)
    ds.lon.attrs['units'] = 'degrees'
    ds.lon.attrs['description'] = 'longitutude'

    ds['lat'] = (('lat',), outLat)
    ds.lat.attrs['units'] = 'degrees'
    ds.lat.attrs['description'] = 'latitutude'

    z = numpy.ma.masked_all(Lon.shape)
    z[yIndices, xIndices] = inZ
    ds['z'] = (('lat', 'lon'), z)
    ds.z.attrs['units'] = 'meters'
    ds.z.attrs['description'] = 'depth of the seafloor (positive up)'

    PT = numpy.ma.masked_all(Lon.shape)
    PT[yIndices, xIndices] = inPT
    ds['botTheta'] = (('lat', 'lon'), PT)
    ds.botTheta.attrs['units'] = '$\degree$C'
    ds.botTheta.attrs['description'] = \
        'potential temperature at sea floor'

    PT_std = numpy.ma.masked_all(Lon.shape)
    # neglect difference between std of PT and CT
    PT_std[yIndices, xIndices] = inCT_std
    ds['botThetaStd'] = (('lat', 'lon'), PT_std)
    ds.botThetaStd.attrs['units'] = '$\degree$C'
    ds.botThetaStd.attrs['description'] = \
        'standard deviation in potential temperature at sea floor'

    S = numpy.ma.masked_all(Lon.shape)
    S[yIndices, xIndices] = inS
    ds['botSalinity'] = (('lat', 'lon'), S)
    ds.botSalinity.attrs['units'] = 'PSU'
    ds.botSalinity.attrs['description'] = \
        'salinity at sea floor'

    S_std = numpy.ma.masked_all(Lon.shape)
    # neglect difference between std of S and SA
    S_std[yIndices, xIndices] = inSA_std
    ds['botSalinityStd'] = (('lat', 'lon'), S_std)
    ds.botSalinityStd.attrs['units'] = 'PSU'
    ds.botSalinityStd.attrs['description'] = \
        'standard deviation in salinity at sea floor'

    PD = numpy.ma.masked_all(Lon.shape)
    PD[yIndices, xIndices] = inPD
    ds['botPotentialDensity'] = (('lat', 'lon'), PD)
    ds.botPotentialDensity.attrs['units'] = 'kg m$^{-3}$'
    ds.botPotentialDensity.attrs['description'] = \
        'potential desnity at sea floor'

    write_netcdf(ds, outFileName)


def remap(inDir, outDir):

    inGridName = 'SouthernOcean_0.25x0.125degree'
    inFileName = '{}/Schmidtko_et_al_2014_bottom_PT_S_PD_{}.nc'.format(
            inDir, inGridName)

    config = MpasAnalysisConfigParser()
    config.read('mpas_analysis/config.default')

    inDescriptor = LatLonGridDescriptor()

    inDescriptor = LatLonGridDescriptor.read(inFileName,
                                             latVarName='lat',
                                             lonVarName='lon')

    outDescriptor = get_comparison_descriptor(config, 'antarctic')
    outGridName = outDescriptor.meshName

    outFileName = '{}/Schmidtko_et_al_2014_bottom_PT_S_PD_{}.nc'.format(
            outDir, outGridName)

    mappingFileName = '{}/map_{}_to_{}.nc'.format(inDir, inGridName,
                                                  outGridName)

    remapper = Remapper(inDescriptor, outDescriptor, mappingFileName)

    remapper.build_mapping_file(method='bilinear')

    if not os.path.exists(outFileName):
        print('Remapping...')
        with xarray.open_dataset(inFileName) as dsIn:
            with remapper.remap(dsIn, renormalizationThreshold=0.01) \
                    as remappedMLD:
                print('Done.')
                remappedMLD.attrs['history'] = ' '.join(sys.argv)
                write_netcdf(remappedMLD, outFileName)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--inDir", dest="inDir", required=True,
                        help="Directory where intermediate files used in "
                             "processing should be downloaded")
    parser.add_argument("-o", "--outDir", dest="outDir", required=True,
                        help="Directory where final preprocessed observation "
                             "are stored")
    args = parser.parse_args()

    try:
        os.makedirs(args.inDir)
    except OSError:
        pass

    try:
        os.makedirs(args.outDir)
    except OSError:
        pass

    urlBase = 'https://www.geomar.de/fileadmin/personal/fb1/po/sschmidtko/'
    fileName = 'Antarctic_shelf_data.txt'

    download_files([fileName], urlBase, args.inDir)
    text_to_netcdf(args.inDir, args.inDir)
    remap(args.inDir, args.outDir)
