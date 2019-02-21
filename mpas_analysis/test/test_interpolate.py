# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2018 Los Alamos National Security, LLC. All rights reserved.
# Copyright (c) 2018 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2018 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
'''
Unit test infrastructure for horizontal interpolation.

Xylar Asay-Davis
04/06/2017
'''

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import pytest
import shutil
import os
import tempfile
import numpy
import xarray
import pyproj

from mpas_analysis.shared.interpolation import Remapper
from mpas_analysis.shared.grid import MpasMeshDescriptor, \
    LatLonGridDescriptor, ProjectionGridDescriptor
from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.configuration import MpasAnalysisConfigParser


@pytest.mark.usefixtures('loaddatadir')
class TestInterp(TestCase):

    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        self.renormalizationThreshold = 0.01

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def get_mpas_descriptor(self):
        mpasMeshFileName = str(self.datadir.join('mpasMesh.nc'))
        timeSeriesFileName = str(self.datadir.join('timeSeries.0002-01-01.nc'))
        descriptor = MpasMeshDescriptor(mpasMeshFileName, meshName='oQU240')

        return (descriptor, mpasMeshFileName, timeSeriesFileName)

    def get_latlon_file_descriptor(self):
        latLonGridFileName = str(self.datadir.join('SST_annual_1870-1900.nc'))
        descriptor = LatLonGridDescriptor.read(latLonGridFileName,
                                               latVarName='lat',
                                               lonVarName='lon')

        return (descriptor, latLonGridFileName)

    def get_latlon_array_descriptor(self):
        configPath = str(self.datadir.join('config.analysis'))
        config = MpasAnalysisConfigParser()
        config.read(configPath)

        lat = numpy.array(config.getExpression('interpolate', 'lat',
                                               usenumpyfunc=True))
        lon = numpy.array(config.getExpression('interpolate', 'lon',
                                               usenumpyfunc=True))

        descriptor = LatLonGridDescriptor.create(lat, lon, units='degrees')
        return descriptor

    def get_stereographic_array_descriptor(self):

        # projection for BEDMAP2 and other common Antarctic data sets
        projection = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 '
                                 '+lon_0=0.0  +k_0=1.0 +x_0=0.0 +y_0=0.0 '
                                 '+ellps=WGS84')

        # a square 61x61 cell map with 100 km resolution and
        xMax = 3000e3
        res = 100e3
        nx = 2 * int(xMax / res) + 1
        x = numpy.linspace(-xMax, xMax, nx)
        meshName = '{}km_Antarctic_stereo'.format(int(res * 1e-3))
        descriptor = ProjectionGridDescriptor.create(projection, x, x,
                                                     meshName)
        return descriptor

    def get_file_names(self, suffix):
        weightFileName = '{}/weights_{}.nc'.format(self.test_dir, suffix)
        outFileName = '{}/remapped_{}.nc'.format(self.test_dir, suffix)
        refFileName = '{}/ref_{}.nc'.format(self.datadir, suffix)
        return (weightFileName, outFileName, refFileName)

    def build_remapper(self, sourceDescriptor, destinationDescriptor,
                       weightFileName):

        remapper = Remapper(sourceDescriptor, destinationDescriptor,
                            weightFileName)

        remapper.build_mapping_file(method='bilinear')

        assert os.path.exists(remapper.mappingFileName)

        return remapper

    def check_remap(self, inFileName, outFileName, refFileName, remapper,
                    remap_file=True):

        dsRef = xarray.open_dataset(refFileName)
        if remap_file:
            # first, test interpolation with ncremap
            remapper.remap_file(inFileName=inFileName,
                                outFileName=outFileName)

            assert os.path.exists(outFileName)
            dsRemapped = xarray.open_dataset(outFileName)
            # drop some extra vairables added by ncremap that aren't in the
            # reference data set
            dsRemapped = dsRemapped.drop(['lat_bnds', 'lon_bnds', 'gw',
                                          'area'])
            self.assertDatasetApproxEqual(dsRemapped, dsRef)

        # now, try in-memory remapping
        ds = xarray.open_dataset(inFileName)
        dsRemapped = remapper.remap(ds, self.renormalizationThreshold)
        self.assertDatasetApproxEqual(dsRemapped, dsRef)

    def test_mpas_to_latlon_file(self):
        '''
        test horizontal interpolation from an MPAS mesh to a destination
        lat/lon grid determined from a file containing 'lat' and 'lon' coords

        Xylar Asay-Davis
        04/06/2017
        '''

        weightFileName, outFileName, refFileName = \
            self.get_file_names(suffix='mpas_to_latlon_file')

        sourceDescriptor, mpasMeshFileName, timeSeriesFileName = \
            self.get_mpas_descriptor()
        destinationDescriptor, latLonGridFileName = \
            self.get_latlon_file_descriptor()

        remapper = self.build_remapper(sourceDescriptor, destinationDescriptor,
                                       weightFileName)
        self.check_remap(timeSeriesFileName, outFileName, refFileName,
                         remapper, remap_file=True)

    def test_mpas_to_latlon_array(self):
        '''
        test horizontal interpolation from an MPAS mesh to a destination
        lat/lon grid determined from config options 'lat' and 'lon'.

        Xylar Asay-Davis
        04/06/2017
        '''

        weightFileName, outFileName, refFileName = \
            self.get_file_names(suffix='mpas_to_latlon_array')

        sourceDescriptor, mpasMeshFileName, timeSeriesFileName = \
            self.get_mpas_descriptor()
        destinationDescriptor = self.get_latlon_array_descriptor()

        remapper = self.build_remapper(sourceDescriptor, destinationDescriptor,
                                       weightFileName)
        self.check_remap(timeSeriesFileName, outFileName, refFileName,
                         remapper, remap_file=True)

    def test_latlon_file_to_latlon_array(self):
        '''
        test horizontal interpolation from a lat/lon grid to a destination
        lat/lon grid determined from config options 'lat' and 'lon'.

        Xylar Asay-Davis
        04/06/2017
        '''

        weightFileName, outFileName, refFileName = \
            self.get_file_names(suffix='latlon_file_to_latlon_array')

        sourceDescriptor, latLonGridFileName = \
            self.get_latlon_file_descriptor()
        destinationDescriptor = self.get_latlon_array_descriptor()

        remapper = self.build_remapper(sourceDescriptor, destinationDescriptor,
                                       weightFileName)
        self.check_remap(latLonGridFileName, outFileName, refFileName,
                         remapper, remap_file=True)

    def test_mpas_to_stereographic_array(self):
        '''
        test horizontal interpolation from an MPAS mesh to a destination
        stereographic grid.

        Xylar Asay-Davis
        04/06/2017
        '''

        weightFileName, outFileName, refFileName = \
            self.get_file_names(suffix='mpas_to_stereographic_array')

        sourceDescriptor, mpasMeshFileName, timeSeriesFileName = \
            self.get_mpas_descriptor()
        destinationDescriptor = self.get_stereographic_array_descriptor()

        remapper = self.build_remapper(sourceDescriptor, destinationDescriptor,
                                       weightFileName)

        # ncremap doesn't support stereographic grids so just check the
        # Remapper
        self.check_remap(timeSeriesFileName, outFileName, refFileName,
                         remapper, remap_file=False)

    def test_latlon_file_to_stereographic_array(self):
        '''
        test horizontal interpolation from a lat/lon grid to a destination
        stereographic grid.

        Xylar Asay-Davis
        04/06/2017
        '''

        weightFileName, outFileName, refFileName = \
            self.get_file_names(suffix='latlon_file_to_stereographic_array')

        sourceDescriptor, latLonGridFileName = \
            self.get_latlon_file_descriptor()
        destinationDescriptor = self.get_stereographic_array_descriptor()

        remapper = self.build_remapper(sourceDescriptor, destinationDescriptor,
                                       weightFileName)

        # ncremap doesn't support stereographic grids so just check the
        # Remapper
        self.check_remap(latLonGridFileName, outFileName, refFileName,
                         remapper, remap_file=False)

    def test_stereographic_array_to_latlon_array(self):
        '''
        test horizontal interpolation from a stereographic grid to a
        destination lat/lon grid determined from config options 'lat' and
        'lon'.

        Xylar Asay-Davis
        04/06/2017
        '''

        weightFileName, outFileName, refFileName = \
            self.get_file_names(suffix='stereographic_array_to_latlon_array')

        sourceDescriptor = self.get_stereographic_array_descriptor()
        destinationDescriptor = self.get_latlon_array_descriptor()

        Lat = sourceDescriptor.coords['lat']['data']

        # now, let's make a more complicated field with more dimensions to
        # check
        inField = numpy.reshape(Lat, (1, Lat.shape[0], Lat.shape[1], 1))
        inField = inField.repeat(3, axis=0)
        inField = inField.repeat(2, axis=3)

        datasetDict = {'dims': ('dim0', 'x', 'y', 'dim3'),
                       'coords': sourceDescriptor.coords,
                       'data_vars': {'complicated':
                                     {'dims': ('dim0', 'x', 'y', 'dim3'),
                                      'data': inField}}}

        ds = xarray.Dataset.from_dict(datasetDict)
        inFileName = '{}/unmapped_stereographic_array_to_latlon_' \
                     'array.nc'.format(self.test_dir)
        ds.to_netcdf(inFileName)

        remapper = self.build_remapper(sourceDescriptor, destinationDescriptor,
                                       weightFileName)

        self.check_remap(inFileName, outFileName, refFileName,
                         remapper, remap_file=False)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
