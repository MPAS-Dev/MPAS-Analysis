# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#
import numpy
import xarray
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
from scipy.io import loadmat
import os

from mpas_analysis.shared.interpolation import Remapper
from mpas_analysis.shared.grid import LatLonGridDescriptor
from mpas_analysis.shared.climatology.climatology \
    import get_antarctic_stereographic_comparison_descriptor
from mpas_analysis.configuration.MpasAnalysisConfigParser \
    import MpasAnalysisConfigParser

from mds import rdmds


def get_bottom_indices(cellFraction):
    nx, ny, nz = cellFraction.shape
    botIndices = -1*numpy.ones((nx, ny), int)
    for zIndex in range(nz):
        mask = cellFraction[:, :, zIndex] > 0.
        botIndices[mask] = zIndex
    return botIndices


def get_monthly_average(filePrefix):
    field, itrs, metadata = rdmds(filePrefix, rec=[0], returnmeta=True)
    nz, ny, nx = field.shape
    # print nx, ny, nz
    yearCount = metadata['nrecords'][0]/12
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
            print '{:04d}-{:02d}'.format(year+2005, month+1)
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


inGridName = 'SouthernOcean_0.167x0.167degree'

inTFileName = '/media/xylar/extra_data/data_overflow/observations/' \
              'SouthernOcean/SOSE/monthly/THETA_mnthlyBar.0000000100'
inSFileName = '/media/xylar/extra_data/data_overflow/observations/' \
              'SouthernOcean/SOSE/monthly/SALT_mnthlyBar.0000000100'
inGridFileName = '/media/xylar/extra_data/data_overflow/observations/' \
                 'SouthernOcean/SOSE/grid.mat'

prefix = 'SOSE_2005-2010_monthly_'

cacheTFileName = '{}_pot_temp_{}.nc'.format(prefix, inGridName)
cacheSFileName = '{}_salinity_{}.nc'.format(prefix, inGridName)
outTFileName = '{}_pot_temp_{}.nc'.format(prefix, outGridName)
outSFileName = '{}_salinity_{}.nc'.format(prefix, outGridName)

config = MpasAnalysisConfigParser()
config.read('config.default')


inDescriptor = LatLonGridDescriptor()

if not os.path.exists(cacheTFileName) or not os.path.exists(cacheSFileName):
    matGrid = loadmat(inGridFileName)
    # lat/lon is a tensor grid so we can use 1-D arrays
    lon = matGrid['XC'][:, 0]
    lat = matGrid['YC'][0, :]
    z = matGrid['RC'][:, 0]
    cellFraction = matGrid['hFacC']

    botIndices = get_bottom_indices(cellFraction)

if os.path.exists(cacheTFileName):
    dsT = xarray.open_dataset(cacheTFileName)
else:
    field, botField = get_monthly_average(inTFileName)

    description = 'Monthly potential temperature climatologies from ' \
                  '2005-2010 average of the Southern Ocean State Estimate ' \
                  '(SOSE)'
    botDescription = 'Monthly potential temperature climatologies at sea ' \
                     'floor from 2005-2010 average from SOSE'
    dictonary = {'dims': ['Time', 'lon', 'lat', 'depth'],
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
                            'depth': {'dims': ('depth'),
                                      'data': z,
                                      'attrs': {'units': 'm'}}},
                 'data_vars': {'theta':
                               {'dims': ('Time', 'lat', 'lon', 'depth'),
                                'data': field,
                                'attrs': {'units': '$^\circ$C',
                                          'description': description}},
                               'botTheta':
                               {'dims': ('Time', 'lat', 'lon'),
                                'data': botField,
                                'attrs': {'units': '$^\circ$C',
                                          'description': botDescription}}}}

    dsT = xarray.Dataset.from_dict(dictonary)
    dsT.to_netcdf(cacheTFileName)

if os.path.exists(cacheSFileName):
    dsS = xarray.open_dataset(cacheSFileName)
else:
    field, botField = get_monthly_average(inSFileName)

    description = 'Monthly salinity climatologies from 2005-2010 ' \
                  'average of the Southern Ocean State Estimate (SOSE)'
    botDescription = 'Monthly salinity climatologies at sea floor ' \
                     'from 2005-2010 average from SOSE'
    dictonary = {'dims': ['Time', 'lon', 'lat', 'depth'],
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
                            'depth': {'dims': ('depth'),
                                      'data': z,
                                      'attrs': {'units': 'm'}}},
                 'data_vars': {'salinity':
                               {'dims': ('Time', 'lat', 'lon', 'depth'),
                                'data': field,
                                'attrs': {'units': 'PSU',
                                          'description': description}},
                               'botSalinity':
                               {'dims': ('Time', 'lat', 'lon'),
                                'data': botField,
                                'attrs': {'units': 'PSU',
                                          'description': botDescription}}}}

    dsS = xarray.Dataset.from_dict(dictonary)
    dsS.to_netcdf(cacheSFileName)

inDescriptor = LatLonGridDescriptor.read(cacheTFileName,  latVarName='lat',
                                         lonVarName='lon')

outDescriptor = get_antarctic_stereographic_comparison_descriptor(config)
outGridName = outDescriptor.meshName

mappingFileName = 'map_{}_to_{}.nc'.format(inGridName, outGridName)

remapper = Remapper(inDescriptor, outDescriptor, mappingFileName)

remapper.build_mapping_file(method='bilinear')

remappedT = remapper.remap(dsT, renormalizationThreshold=0.01)

remappedT.attrs['history'] = ' '.join(sys.argv)
remappedT.to_netcdf(outTFileName)

remappedS = remapper.remap(dsS, renormalizationThreshold=0.01)

remappedS.attrs['history'] = ' '.join(sys.argv)
remappedS.to_netcdf(outSFileName)

normT = colors.Normalize(vmin=-2.0, vmax=2.0)
normS = colors.Normalize(vmin=33.0, vmax=35.0)

plt.figure()
plt.imshow(remappedT.botTheta.values[0, :, :], origin='lower', cmap='RdBu_r',
           norm=normT)
plt.colorbar()
plt.figure()
plt.imshow(remappedS.botSalinity.values[0, :, :], origin='lower',
           cmap='RdBu_r', norm=normS)
plt.colorbar()

plt.show()
