#!/usr/bin/env python
# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#

"""
A script for downloading and preprocessing data sets from SOSE for use in
MPAS-Analysis
"""
# Authors
# -------
# Xylar Asay-Davis

from __future__ import absolute_import, division, print_function, \
    unicode_literals


import numpy
import xarray
import sys
from scipy.io import loadmat
import os
import argparse
import gzip
import shutil
import gsw

from mpas_analysis.shared.io.download import download_files
from mpas_analysis.shared.interpolation import Remapper
from mpas_analysis.shared.grid import LatLonGridDescriptor
from mpas_analysis.shared.climatology.comparison_descriptors \
    import get_comparison_descriptor
from mpas_analysis.configuration \
    import MpasAnalysisConfigParser
from mpas_analysis.shared.io import write_netcdf

from mds import rdmds


def get_bottom_indices(cellFraction):
    nx, ny, nz = cellFraction.shape
    botIndices = -1*numpy.ones((nx, ny), int)
    for zIndex in range(nz):
        mask = cellFraction[:, :, zIndex] > 0.
        botIndices[mask] = zIndex
    return botIndices


def get_monthly_average_3d(filePrefix, cellFraction, botIndices):
    field, itrs, metadata = rdmds(filePrefix, rec=[0], returnmeta=True)
    nz, ny, nx = field.shape
    # print nx, ny, nz
    yearCount = metadata['nrecords'][0]//12
    dims = [12, nx, ny, nz]

    mask3D = cellFraction <= 0.
    mask2D = botIndices == -1
    xIndices, yIndices = numpy.meshgrid(numpy.arange(nx), numpy.arange(ny),
                                        indexing='ij')
    monthlyClimatologies = numpy.ma.masked_all(dims)
    botMonthlyClimatologies = numpy.ma.masked_all((12, nx, ny))
    for month in range(12):
        first = True
        for year in range(yearCount):
            print('{:04d}-{:02d}'.format(year+2005, month+1))
            recordIndex = year*12 + month
            field = rdmds(filePrefix, rec=[recordIndex])
            field = field.transpose(2, 1, 0)

            field = numpy.ma.masked_array(field, mask=mask3D)
            if first:
                monthlyClimatologies[month, :, :, :] = field/float(yearCount)
                first = False
            else:
                monthlyClimatologies[month, :, :, :] = \
                    monthlyClimatologies[month, :, :, :] + \
                    field/float(yearCount)
        botMonthlyClimatologies[month, :, :] = \
            numpy.ma.masked_array(field[xIndices, yIndices, botIndices],
                                  mask=mask2D)

    monthlyClimatologies = monthlyClimatologies.transpose(0, 2, 1, 3)
    botMonthlyClimatologies = botMonthlyClimatologies.transpose(0, 2, 1)

    return monthlyClimatologies, botMonthlyClimatologies


def get_monthly_average_2d(filePrefix, botIndices):
    field, itrs, metadata = rdmds(filePrefix, rec=[0], returnmeta=True)
    ny, nx = field.shape
    # print nx, ny, nz
    yearCount = metadata['nrecords'][0]//12
    dims = [12, nx, ny]

    mask2D = botIndices == -1
    xIndices, yIndices = numpy.meshgrid(numpy.arange(nx), numpy.arange(ny),
                                        indexing='ij')
    monthlyClimatologies = numpy.ma.masked_all(dims)
    for month in range(12):
        first = True
        for year in range(yearCount):
            print('{:04d}-{:02d}'.format(year+2005, month+1))
            recordIndex = year*12 + month
            field = rdmds(filePrefix, rec=[recordIndex])
            field = field.transpose(1, 0)

            field = numpy.ma.masked_array(field, mask=mask2D)
            if first:
                monthlyClimatologies[month, :, :] = field/float(yearCount)
                first = False
            else:
                monthlyClimatologies[month, :, :] = \
                    monthlyClimatologies[month, :, :] + \
                    field/float(yearCount)

    monthlyClimatologies = monthlyClimatologies.transpose(0, 2, 1)
    return monthlyClimatologies


def unzip_sose_data(inPrefixes, outDir):
    # unzip the gzipped data
    for prefix in inPrefixes:
        zippedFileName = '{}/{}.data.gz'.format(outDir, prefix)
        unzippedFileName = '{}/{}.data'.format(outDir, prefix)
        if not os.path.exists(unzippedFileName):
            print('Unzipping {}.data.gz...'.format(prefix))
            with gzip.open(zippedFileName, 'rb') as f_in:
                with open(unzippedFileName, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print('Done.')


def sose_pt_to_nc(inPrefix, outFileName, lon, lat, z, cellFraction,
                  botIndices):
    if os.path.exists(outFileName):
        dsT = xarray.open_dataset(outFileName)
    else:

        print('Building climatology of potential temperature...')
        field, botField = get_monthly_average_3d(inPrefix, cellFraction,
                                                 botIndices)
        print('Done.')
        zBot = numpy.ma.masked_array(z[botIndices], mask=(botIndices == -1))
        zBot = zBot.transpose(1, 0)

        description = 'Monthly potential temperature climatologies from ' \
                      '2005-2010 average of the Southern Ocean State ' \
                      'Estimate (SOSE)'
        botDescription = 'Monthly potential temperature climatologies at ' \
                         'seafloor from 2005-2010 average from SOSE'
        dictonary = {'dims': ['Time', 'lon', 'lat', 'z'],
                     'coords': {'month': {'dims': ('Time'),
                                          'data': range(1, 13),
                                          'attrs': {'units': 'months'}},
                                'year': {'dims': ('Time'),
                                         'data': numpy.ones(12),
                                         'attrs': {'units': 'years'}},
                                'lon': {'dims': ('lon'),
                                        'data': lon,
                                        'attrs': {'units': 'degrees'}},
                                'lat': {'dims': ('lat'),
                                        'data': lat,
                                        'attrs': {'units': 'degrees'}},
                                'z': {'dims': ('z'),
                                      'data': z,
                                      'attrs': {'units': 'm'}},
                                'zBot': {'dims': ('lat', 'lon'),
                                         'data': zBot,
                                         'attrs': {'units': 'm'}}},
                     'data_vars': {'theta':
                                   {'dims': ('Time', 'lat', 'lon', 'z'),
                                    'data': field,
                                    'attrs': {'units': '$\degree$C',
                                              'description': description}},
                                   'botTheta':
                                   {'dims': ('Time', 'lat', 'lon'),
                                    'data': botField,
                                    'attrs': {'units': '$\degree$C',
                                              'description': botDescription}}}}

        dsT = xarray.Dataset.from_dict(dictonary)
        write_netcdf(dsT, outFileName)
    return dsT


def sose_s_to_nc(inPrefix, outFileName, lon, lat, z, cellFraction, botIndices):
    if os.path.exists(outFileName):
        dsS = xarray.open_dataset(outFileName)
    else:
        print('Building climatology of salinity...')
        field, botField = get_monthly_average_3d(inPrefix, cellFraction,
                                                 botIndices)
        print('Done.')
        zBot = numpy.ma.masked_array(z[botIndices], mask=(botIndices == -1))
        zBot = zBot.transpose(1, 0)

        description = 'Monthly salinity climatologies from 2005-2010 ' \
                      'average of the Southern Ocean State Estimate (SOSE)'
        botDescription = 'Monthly salinity climatologies at sea floor ' \
                         'from 2005-2010 average from SOSE'
        dictonary = {'dims': ['Time', 'lon', 'lat', 'z'],
                     'coords': {'month': {'dims': ('Time'),
                                          'data': range(1, 13),
                                          'attrs': {'units': 'months'}},
                                'year': {'dims': ('Time'),
                                         'data': numpy.ones(12),
                                         'attrs': {'units': 'years'}},
                                'lon': {'dims': ('lon'),
                                        'data': lon,
                                        'attrs': {'units': 'degrees'}},
                                'lat': {'dims': ('lat'),
                                        'data': lat,
                                        'attrs': {'units': 'degrees'}},
                                'z': {'dims': ('z'),
                                      'data': z,
                                      'attrs': {'units': 'm'}},
                                'zBot': {'dims': ('lat', 'lon'),
                                         'data': zBot,
                                         'attrs': {'units': 'm'}}},
                     'data_vars': {'salinity':
                                   {'dims': ('Time', 'lat', 'lon', 'z'),
                                    'data': field,
                                    'attrs': {'units': 'PSU',
                                              'description': description}},
                                   'botSalinity':
                                   {'dims': ('Time', 'lat', 'lon'),
                                    'data': botField,
                                    'attrs': {'units': 'PSU',
                                              'description': botDescription}}}}

        dsS = xarray.Dataset.from_dict(dictonary)
        write_netcdf(dsS, outFileName)

    return dsS


def sose_gammaN_to_nc(inPrefix, outFileName, lon, lat, z, cellFraction,
                      botIndices):
    if os.path.exists(outFileName):
        dsS = xarray.open_dataset(outFileName)
    else:
        print('Building climatology of neutral density...')
        field, botField = get_monthly_average_3d(inPrefix, cellFraction,
                                                 botIndices)
        print('Done.')
        zBot = numpy.ma.masked_array(z[botIndices], mask=(botIndices == -1))
        zBot = zBot.transpose(1, 0)

        description = 'Monthly neutral density climatologies from 2005-2010 ' \
                      'average of the Southern Ocean State Estimate (SOSE)'
        botDescription = 'Monthly neutral density climatologies at sea ' \
                         'floor from 2005-2010 average from SOSE'
        dictonary = {'dims': ['Time', 'lon', 'lat', 'z'],
                     'coords': {'month': {'dims': ('Time'),
                                          'data': range(1, 13),
                                          'attrs': {'units': 'months'}},
                                'year': {'dims': ('Time'),
                                         'data': numpy.ones(12),
                                         'attrs': {'units': 'years'}},
                                'lon': {'dims': ('lon'),
                                        'data': lon,
                                        'attrs': {'units': 'degrees'}},
                                'lat': {'dims': ('lat'),
                                        'data': lat,
                                        'attrs': {'units': 'degrees'}},
                                'z': {'dims': ('z'),
                                      'data': z,
                                      'attrs': {'units': 'm'}},
                                'zBot': {'dims': ('lat', 'lon'),
                                         'data': zBot,
                                         'attrs': {'units': 'm'}}},
                     'data_vars': {'neutralDensity':
                                   {'dims': ('Time', 'lat', 'lon', 'z'),
                                    'data': field,
                                    'attrs': {'units': 'kg m$^{-3}$',
                                              'description': description}},
                                   'botNeutralDensity':
                                   {'dims': ('Time', 'lat', 'lon'),
                                    'data': botField,
                                    'attrs': {'units': 'kg m$^{-3}$',
                                              'description': botDescription}}}}

        dsS = xarray.Dataset.from_dict(dictonary)
        write_netcdf(dsS, outFileName)

    return dsS


def sose_mld_to_nc(inPrefix, outFileName, lon, lat, botIndices):
    if os.path.exists(outFileName):
        dsMLD = xarray.open_dataset(outFileName)
    else:
        print('Building climatology of mixed layer depth...')
        field = get_monthly_average_2d(inPrefix, botIndices)
        # make MLD positive
        field = -field
        print('Done.')

        description = 'Monthly mixed layer depth climatologies from ' \
                      '2005-2010 average of the Southern Ocean State ' \
                      'Estimate (SOSE)'
        dictonary = {'dims': ['Time', 'lon', 'lat'],
                     'coords': {'month': {'dims': ('Time'),
                                          'data': range(1, 13),
                                          'attrs': {'units': 'months'}},
                                'year': {'dims': ('Time'),
                                         'data': numpy.ones(12),
                                         'attrs': {'units': 'years'}},
                                'lon': {'dims': ('lon'),
                                        'data': lon,
                                        'attrs': {'units': 'degrees'}},
                                'lat': {'dims': ('lat'),
                                        'data': lat,
                                        'attrs': {'units': 'degrees'}}},
                     'data_vars': {'mld':
                                   {'dims': ('Time', 'lat', 'lon'),
                                    'data': field,
                                    'attrs': {'units': 'm',
                                              'description': description}}}}

        dsMLD = xarray.Dataset.from_dict(dictonary)
        write_netcdf(dsMLD, outFileName)

    return dsMLD


def sose_u_to_nc(inPrefix, outFileName, lon, lat, z, cellFraction, botIndices):
    if os.path.exists(outFileName):
        dsU = xarray.open_dataset(outFileName)
    else:
        print('Building climatology of zonal velocity...')
        field, botField = get_monthly_average_3d(inPrefix, cellFraction,
                                                 botIndices)
        print('Done.')
        zBot = numpy.ma.masked_array(z[botIndices], mask=(botIndices == -1))
        zBot = zBot.transpose(1, 0)

        description = 'Monthly zonal velocity climatologies from 2005-2010 ' \
                      'average of the Southern Ocean State Estimate (SOSE)'
        botDescription = 'Monthly zonal velocity climatologies at sea floor ' \
                         'from 2005-2010 average from SOSE'
        dictonary = {'dims': ['Time', 'lon', 'lat', 'z'],
                     'coords': {'month': {'dims': ('Time'),
                                          'data': range(1, 13),
                                          'attrs': {'units': 'months'}},
                                'year': {'dims': ('Time'),
                                         'data': numpy.ones(12),
                                         'attrs': {'units': 'years'}},
                                'lon': {'dims': ('lon'),
                                        'data': lon,
                                        'attrs': {'units': 'degrees'}},
                                'lat': {'dims': ('lat'),
                                        'data': lat,
                                        'attrs': {'units': 'degrees'}},
                                'z': {'dims': ('z'),
                                      'data': z,
                                      'attrs': {'units': 'm'}},
                                'zBot': {'dims': ('lat', 'lon'),
                                         'data': zBot,
                                         'attrs': {'units': 'm'}}},
                     'data_vars': {'zonalVel':
                                   {'dims': ('Time', 'lat', 'lon', 'z'),
                                    'data': field,
                                    'attrs': {'units': 'm s$^{-1}$',
                                              'description': description}},
                                   'botZonalVel':
                                   {'dims': ('Time', 'lat', 'lon'),
                                    'data': botField,
                                    'attrs': {'units': 'm s$^{-1}$',
                                              'description': botDescription}}}}

        dsU = xarray.Dataset.from_dict(dictonary)
        write_netcdf(dsU, outFileName)

    return dsU


def sose_v_to_nc(inPrefix, outFileName, lon, lat, z, cellFraction, botIndices):
    if os.path.exists(outFileName):
        dsV = xarray.open_dataset(outFileName)
    else:
        print('Building climatology of meridional velocity...')
        field, botField = get_monthly_average_3d(inPrefix, cellFraction,
                                                 botIndices)
        print('Done.')
        zBot = numpy.ma.masked_array(z[botIndices], mask=(botIndices == -1))
        zBot = zBot.transpose(1, 0)

        description = 'Monthly meridional velocity climatologies from ' \
                      '2005-2010 average of the Southern Ocean State ' \
                      'Estimate (SOSE)'
        botDescription = 'Monthly meridional velocity climatologies at sea ' \
                         'floor from 2005-2010 average from SOSE'
        dictonary = {'dims': ['Time', 'lon', 'lat', 'z'],
                     'coords': {'month': {'dims': ('Time'),
                                          'data': range(1, 13),
                                          'attrs': {'units': 'months'}},
                                'year': {'dims': ('Time'),
                                         'data': numpy.ones(12),
                                         'attrs': {'units': 'years'}},
                                'lon': {'dims': ('lon'),
                                        'data': lon,
                                        'attrs': {'units': 'degrees'}},
                                'lat': {'dims': ('lat'),
                                        'data': lat,
                                        'attrs': {'units': 'degrees'}},
                                'z': {'dims': ('z'),
                                      'data': z,
                                      'attrs': {'units': 'm'}},
                                'zBot': {'dims': ('lat', 'lon'),
                                         'data': zBot,
                                         'attrs': {'units': 'm'}}},
                     'data_vars': {'meridVel':
                                   {'dims': ('Time', 'lat', 'lon', 'z'),
                                    'data': field,
                                    'attrs': {'units': 'm s$^{-1}$',
                                              'description': description}},
                                   'botMeridVel':
                                   {'dims': ('Time', 'lat', 'lon'),
                                    'data': botField,
                                    'attrs': {'units': 'm s$^{-1}$',
                                              'description': botDescription}}}}

        dsV = xarray.Dataset.from_dict(dictonary)
        write_netcdf(dsV, outFileName)

    return dsV


def remap_pt_s(prefix, inGridName, inGridFileName, inDir, inTPrefix,
               inSPrefix, inGammaNPrefix):
    cacheTFileName = '{}_pot_temp_{}.nc'.format(prefix, inGridName)
    cacheSFileName = '{}_salinity_{}.nc'.format(prefix, inGridName)
    cacheGammaNFileName = '{}_neut_den_{}.nc'.format(prefix, inGridName)

    config = MpasAnalysisConfigParser()
    config.read('mpas_analysis/config.default')

    matGrid = loadmat(inGridFileName)
    # lat/lon is a tensor grid so we can use 1-D arrays
    lon = matGrid['XC'][:, 0]
    lat = matGrid['YC'][0, :]
    z = matGrid['RC'][:, 0]
    cellFraction = matGrid['hFacC']

    botIndices = get_bottom_indices(cellFraction)

    with sose_pt_to_nc('{}/{}'.format(inDir, inTPrefix),
                       cacheTFileName, lon, lat, z, cellFraction, botIndices) \
            as dsT:
        inDescriptor = LatLonGridDescriptor()

        inDescriptor = LatLonGridDescriptor.read(cacheTFileName,
                                                 latVarName='lat',
                                                 lonVarName='lon')

        outDescriptor = get_comparison_descriptor(config, 'antarctic')
        outGridName = outDescriptor.meshName

        outTFileName = '{}_pot_temp_{}.nc'.format(prefix, outGridName)
        outSFileName = '{}_salinity_{}.nc'.format(prefix, outGridName)
        outGammaNFileName = '{}_neut_den_{}.nc'.format(prefix, outGridName)

        mappingFileName = '{}/map_C_{}_to_{}.nc'.format(inDir, inGridName,
                                                        outGridName)

        remapper = Remapper(inDescriptor, outDescriptor, mappingFileName)

        remapper.build_mapping_file(method='bilinear')

        if not os.path.exists(outTFileName):
            dsT.reset_coords(names='zBot', inplace=True)
            print('Remapping potential temperature...')
            with remapper.remap(dsT, renormalizationThreshold=0.01) \
                    as remappedT:
                print('Done.')
                remappedT.attrs['history'] = ' '.join(sys.argv)
                remappedT.set_coords(names='zBot', inplace=True)
                write_netcdf(remappedT, outTFileName)

    with sose_s_to_nc('{}/{}'.format(inDir, inSPrefix),
                      cacheSFileName, lon, lat, z, cellFraction, botIndices) \
            as dsS:
        if not os.path.exists(outSFileName):
            dsS.reset_coords(names='zBot', inplace=True)
            print('Remapping salinity...')
            with remapper.remap(dsS, renormalizationThreshold=0.01) \
                    as remappedS:
                print('Done.')
                remappedS.attrs['history'] = ' '.join(sys.argv)
                remappedS.set_coords(names='zBot', inplace=True)
                write_netcdf(remappedS, outSFileName)

    with sose_gammaN_to_nc('{}/{}'.format(inDir, inGammaNPrefix),
                           cacheGammaNFileName, lon, lat, z, cellFraction,
                           botIndices) \
            as dsGammaN:
        if not os.path.exists(outGammaNFileName):
            dsGammaN.reset_coords(names='zBot', inplace=True)
            print('Remapping neutral density...')
            with remapper.remap(dsGammaN, renormalizationThreshold=0.01) \
                    as remappedGammaN:
                print('Done.')
                remappedGammaN.attrs['history'] = ' '.join(sys.argv)
                remappedGammaN.set_coords(names='zBot', inplace=True)
                write_netcdf(remappedGammaN, outGammaNFileName)


def remap_mld(prefix, inGridName, inGridFileName, inDir, inMLDPrefix):
    cacheMLDFileName = '{}_mld_{}.nc'.format(prefix, inGridName)

    config = MpasAnalysisConfigParser()
    config.read('mpas_analysis/config.default')

    matGrid = loadmat(inGridFileName)
    # lat/lon is a tensor grid so we can use 1-D arrays
    lon = matGrid['XC'][:, 0]
    lat = matGrid['YC'][0, :]
    cellFraction = matGrid['hFacC']

    botIndices = get_bottom_indices(cellFraction)

    with sose_mld_to_nc('{}/{}'.format(inDir, inMLDPrefix),
                        cacheMLDFileName, lon, lat, botIndices) as dsMLD:
        inDescriptor = LatLonGridDescriptor()

        inDescriptor = LatLonGridDescriptor.read(cacheMLDFileName,
                                                 latVarName='lat',
                                                 lonVarName='lon')

        outDescriptor = get_comparison_descriptor(config, 'antarctic')
        outGridName = outDescriptor.meshName

        outMLDFileName = '{}_mld_{}.nc'.format(prefix, outGridName)

        mappingFileName = '{}/map_C_{}_to_{}.nc'.format(inDir, inGridName,
                                                        outGridName)

        remapper = Remapper(inDescriptor, outDescriptor, mappingFileName)

        remapper.build_mapping_file(method='bilinear')

        if not os.path.exists(outMLDFileName):
            print('Remapping mixed layer depth...')
            with remapper.remap(dsMLD, renormalizationThreshold=0.01) \
                    as remappedMLD:
                print('Done.')
                remappedMLD.attrs['history'] = ' '.join(sys.argv)
                write_netcdf(remappedMLD, outMLDFileName)


def remap_u(prefix, inGridName, inGridFileName, inDir, inUPrefix):
    cacheUFileName = '{}_zonal_vel_{}.nc'.format(prefix, inGridName)

    config = MpasAnalysisConfigParser()
    config.read('mpas_analysis/config.default')

    matGrid = loadmat(inGridFileName)
    # lat/lon is a tensor grid so we can use 1-D arrays
    lon = matGrid['XG'][:, 0]
    lat = matGrid['YC'][0, :]
    z = matGrid['RC'][:, 0]
    cellFraction = matGrid['hFacW']

    botIndices = get_bottom_indices(cellFraction)

    with sose_u_to_nc('{}/{}'.format(inDir, inUPrefix),
                      cacheUFileName, lon, lat, z, cellFraction, botIndices) \
            as dsU:
        inDescriptor = LatLonGridDescriptor()

        inDescriptor = LatLonGridDescriptor.read(cacheUFileName,
                                                 latVarName='lat',
                                                 lonVarName='lon')

        outDescriptor = get_comparison_descriptor(config, 'antarctic')
        outGridName = outDescriptor.meshName

        outUFileName = '{}_zonal_vel_{}.nc'.format(prefix, outGridName)

        mappingFileName = '{}/map_U_{}_to_{}.nc'.format(inDir, inGridName,
                                                        outGridName)

        remapper = Remapper(inDescriptor, outDescriptor, mappingFileName)

        remapper.build_mapping_file(method='bilinear')

        if not os.path.exists(outUFileName):
            print('Remapping zonal velocity...')
            with remapper.remap(dsU, renormalizationThreshold=0.01) \
                    as remappedU:
                print('Done.')
                remappedU.attrs['history'] = ' '.join(sys.argv)
                write_netcdf(remappedU, outUFileName)


def remap_v(prefix, inGridName, inGridFileName, inDir, inVPrefix):
    cacheVFileName = '{}_merid_vel_{}.nc'.format(prefix, inGridName)

    config = MpasAnalysisConfigParser()
    config.read('mpas_analysis/config.default')

    matGrid = loadmat(inGridFileName)
    # lat/lon is a tensor grid so we can use 1-D arrays
    lon = matGrid['XC'][:, 0]
    lat = matGrid['YG'][0, :]
    z = matGrid['RC'][:, 0]
    cellFraction = matGrid['hFacS']

    botIndices = get_bottom_indices(cellFraction)

    with sose_v_to_nc('{}/{}'.format(inDir, inVPrefix),
                      cacheVFileName, lon, lat, z, cellFraction, botIndices) \
            as dsV:
        inDescriptor = LatLonGridDescriptor()

        inDescriptor = LatLonGridDescriptor.read(cacheVFileName,
                                                 latVarName='lat',
                                                 lonVarName='lon')

        outDescriptor = get_comparison_descriptor(config, 'antarctic')
        outGridName = outDescriptor.meshName

        outVFileName = '{}_merid_vel_{}.nc'.format(prefix, outGridName)

        mappingFileName = '{}/map_V_{}_to_{}.nc'.format(inDir, inGridName,
                                                        outGridName)

        remapper = Remapper(inDescriptor, outDescriptor, mappingFileName)

        remapper.build_mapping_file(method='bilinear')

        if not os.path.exists(outVFileName):
            print('Remapping meridional velocity...')
            with remapper.remap(dsV, renormalizationThreshold=0.01) \
                    as remappedV:
                print('Done.')
                remappedV.attrs['history'] = ' '.join(sys.argv)
                write_netcdf(remappedV, outVFileName)


def compute_vel_mag(prefix, inGridName, inDir):
    config = MpasAnalysisConfigParser()
    config.read('mpas_analysis/config.default')

    outDescriptor = get_comparison_descriptor(config, 'antarctic')
    outGridName = outDescriptor.meshName
    description = 'Monthly velocity magnitude climatologies from ' \
                  '2005-2010 average of the Southern Ocean State ' \
                  'Estimate (SOSE)'
    botDescription = 'Monthly velocity magnitude climatologies at sea ' \
                     'floor from 2005-2010 average from SOSE'

    for gridName in [outGridName]:
        outFileName = '{}_vel_mag_{}.nc'.format(prefix, gridName)
        uFileName = '{}_zonal_vel_{}.nc'.format(prefix, gridName)
        vFileName = '{}_merid_vel_{}.nc'.format(prefix, gridName)
        if not os.path.exists(outFileName):
            with xarray.open_dataset(uFileName) as dsU:
                with xarray.open_dataset(vFileName) as dsV:
                    dsVelMag = dsU.drop(['zonalVel', 'botZonalVel'])
                    dsVelMag['velMag'] = xarray.ufuncs.sqrt(
                            dsU.zonalVel**2 + dsV.meridVel**2)
                    dsVelMag.velMag.attrs['units'] = 'm s$^{-1}$'
                    dsVelMag.velMag.attrs['description'] = description

                    dsVelMag['botVelMag'] = xarray.ufuncs.sqrt(
                            dsU.botZonalVel**2 + dsV.botMeridVel**2)
                    dsVelMag.botVelMag.attrs['units'] = 'm s$^{-1}$'
                    dsVelMag.botVelMag.attrs['description'] = botDescription

                    write_netcdf(dsVelMag, outFileName)


def compute_pot_density(prefix, inGridName, inDir):
    config = MpasAnalysisConfigParser()
    config.read('mpas_analysis/config.default')

    outDescriptor = get_comparison_descriptor(config, 'antarctic')
    outGridName = outDescriptor.meshName
    description = 'Monthly potential density climatologies from ' \
                  '2005-2010 average of the Southern Ocean State ' \
                  'Estimate (SOSE)'
    botDescription = 'Monthly potential density climatologies at sea ' \
                     'floor from 2005-2010 average from SOSE'

    for gridName in [inGridName, outGridName]:
        outFileName = '{}_pot_den_{}.nc'.format(prefix, gridName)
        TFileName = '{}_pot_temp_{}.nc'.format(prefix, gridName)
        SFileName = '{}_salinity_{}.nc'.format(prefix, gridName)
        if not os.path.exists(outFileName):
            with xarray.open_dataset(TFileName) as dsT:
                with xarray.open_dataset(SFileName) as dsS:
                    dsPotDensity = dsT.drop(['theta', 'botTheta'])

                    lat, lon, z = xarray.broadcast(dsS.lat, dsS.lon, dsS.z)
                    pressure = gsw.p_from_z(z.values, lat.values)
                    SA = gsw.SA_from_SP(dsS.salinity.values, pressure,
                                        lon.values, lat.values)
                    CT = gsw.CT_from_pt(SA, dsT.theta.values)
                    dsPotDensity['potentialDensity'] = (dsS.salinity.dims,
                                                        gsw.rho(SA, CT, 0.))
                    dsPotDensity.potentialDensity.attrs['units'] = \
                        'kg m$^{-3}$'
                    dsPotDensity.potentialDensity.attrs['description'] = \
                        description

                    lat, lon, z = xarray.broadcast(dsS.lat, dsS.lon, dsS.zBot)
                    pressure = gsw.p_from_z(z.values, lat.values)
                    SA = gsw.SA_from_SP(dsS.botSalinity.values, pressure,
                                        lon.values, lat.values)
                    CT = gsw.CT_from_pt(SA, dsT.botTheta.values)
                    dsPotDensity['botPotentialDensity'] = \
                        (dsS.botSalinity.dims, gsw.rho(SA, CT, 0.))
                    dsPotDensity.botPotentialDensity.attrs['units'] = \
                        'kg m$^{-3}$'
                    dsPotDensity.botPotentialDensity.attrs['description'] = \
                        botDescription

                    write_netcdf(dsPotDensity, outFileName)


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

    inGridName = 'SouthernOcean_0.167x0.167degree'

    inTPrefix = 'THETA_mnthlyBar.0000000100'
    inSPrefix = 'SALT_mnthlyBar.0000000100'
    inMLDPrefix = 'MLD_mnthlyBar.0000000100'
    inUPrefix = 'UVEL_mnthlyBar.0000000100'
    inVPrefix = 'VVEL_mnthlyBar.0000000100'
    inGammaNPrefix = 'GAMMA_mnthlyBar.0000000100'

    inPrefixes = [inTPrefix, inSPrefix, inMLDPrefix, inUPrefix, inVPrefix,
                  inGammaNPrefix]

    inGridFileName = '{}/grid.mat'.format(args.inDir)

    try:
        os.makedirs(args.inDir)
    except OSError:
        pass

    # dowload the desired file
    download_files(['GRID_README.txt'], urlBase='http://sose.ucsd.edu/DATA',
                   outDir=args.inDir)

    urlBase = 'http://sose.ucsd.edu/DATA/SO6_V2'

    fileList = ['grid.mat']
    for prefix in inPrefixes:
        fileList.append('{}.data.gz'.format(prefix))
        fileList.append('{}.meta'.format(prefix))

    download_files(fileList, urlBase, args.inDir)
    unzip_sose_data(inPrefixes, args.inDir)

    prefix = '{}/SOSE_2005-2010_monthly'.format(args.outDir)

    remap_pt_s(prefix, inGridName, inGridFileName, args.inDir, inTPrefix,
               inSPrefix, inGammaNPrefix)

    remap_mld(prefix, inGridName, inGridFileName, args.inDir, inMLDPrefix)

    remap_u(prefix, inGridName, inGridFileName, args.inDir, inUPrefix)
    remap_v(prefix, inGridName, inGridFileName, args.inDir, inVPrefix)

    compute_vel_mag(prefix, inGridName, args.inDir)

    compute_pot_density(prefix, inGridName, args.inDir)
