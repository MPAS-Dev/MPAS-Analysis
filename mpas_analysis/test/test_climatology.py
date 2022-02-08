# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2020 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2020 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2020 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
"""
Unit test infrastructure for climatologies.

Xylar Asay-Davis
04/11/2017
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import pytest
import tempfile
import shutil
import os
import numpy
import xarray
from pyremap import MpasMeshDescriptor, LatLonGridDescriptor

from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.shared.generalized_reader.generalized_reader \
    import open_multifile_dataset
from mpas_analysis.configuration import MpasAnalysisConfigParser
from mpas_analysis.shared.climatology import \
    get_comparison_descriptor, get_remapper, \
    add_years_months_days_in_month, compute_climatology, \
    compute_monthly_climatology
from mpas_analysis.shared.constants import constants


@pytest.mark.usefixtures("loaddatadir")
class TestClimatology(TestCase):

    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def setup_config(self, maxChunkSize=10000):
        config = MpasAnalysisConfigParser()

        config.add_section('execute')
        config.set('execute', 'mapParallelExec', 'None')

        config.add_section('diagnostics')
        config.set('diagnostics', 'base_path', self.test_dir)
        config.set('diagnostics', 'customDirectory', 'none')
        config.set('diagnostics', 'mappingSubdirectory', 'maps')

        config.add_section('input')
        config.set('input', 'maxChunkSize', str(maxChunkSize))
        config.set('input', 'mpasMeshName', 'QU240')

        config.add_section('output')
        config.set('output', 'baseDirectory', self.test_dir)
        config.set('output', 'mappingSubdirectory', '.')
        config.set('output', 'mpasClimatologySubdirectory', 'clim/mpas')

        config.add_section('climatology')
        config.set('climatology', 'startYear', '2')
        config.set('climatology', 'endYear', '2')
        config.set('climatology', 'comparisonLatResolution', '0.5')
        config.set('climatology', 'comparisonLonResolution', '0.5')

        config.set('climatology', 'mpasInterpolationMethod', 'bilinear')

        config.add_section('oceanObservations')
        config.set('oceanObservations', 'interpolationMethod', 'bilinear')
        config.set('oceanObservations', 'climatologySubdirectory', 'clim/obs')
        config.set('oceanObservations', 'remappedClimSubdirectory',
                   'clim/obs/remap')

        return config

    def setup_mpas_remapper(self, config):
        mpasMeshFileName = '{}/mpasMesh.nc'.format(self.datadir)

        comparisonDescriptor = \
            get_comparison_descriptor(config, comparison_grid_name='latlon')

        mpasDescriptor = MpasMeshDescriptor(
            mpasMeshFileName, meshName=config.get('input', 'mpasMeshName'))

        remapper = get_remapper(
            config=config, sourceDescriptor=mpasDescriptor,
            comparisonDescriptor=comparisonDescriptor,
            mappingFilePrefix='map', method=config.get(
                'climatology', 'mpasInterpolationMethod'))

        return remapper

    def setup_obs_remapper(self, config, fieldName):
        gridFileName = '{}/obsGrid.nc'.format(self.datadir)

        comparisonDescriptor = \
            get_comparison_descriptor(config, comparison_grid_name='latlon')

        obsDescriptor = LatLonGridDescriptor.read(fileName=gridFileName,
                                                  latVarName='lat',
                                                  lonVarName='lon')

        remapper = \
            get_remapper(
                config=config, sourceDescriptor=obsDescriptor,
                comparisonDescriptor=comparisonDescriptor,
                mappingFilePrefix='map_obs_{}'.format(fieldName),
                method=config.get('oceanObservations',
                                  'interpolationMethod'))

        return remapper

    def open_test_ds(self, config, calendar):
        fileNames = ['{}/timeSeries.0002-{:02d}-01.nc'.format(self.datadir,
                                                              month)
                     for month in [1, 2, 3]]

        variableMap = {'mld': ['timeMonthly_avg_tThreshMLD'],
                       'Time': [['xtime_startMonthly', 'xtime_endMonthly']]}
        variableList = ['mld']

        ds = open_multifile_dataset(
            fileNames=fileNames,
            calendar=calendar,
            config=config,
            timeVariableName='Time',
            variableList=variableList,
            variableMap=variableMap)

        assert(len(ds.Time) == 3)
        return ds

    def test_get_mpas_remapper(self):
        config = self.setup_config()

        explicitMappingPath = '{}/maps'.format(self.test_dir)
        os.makedirs(explicitMappingPath)

        fileBase = 'map_QU240_to_0.5x0.5degree_bilinear.nc'

        defaultMappingFileName = '{}/{}'.format(self.test_dir, fileBase)
        explicitMappingFileName = '{}/{}'.format(explicitMappingPath, fileBase)

        for mappingFileName, setName in [(defaultMappingFileName, False),
                                         (explicitMappingFileName, True)]:

            remapper = self.setup_mpas_remapper(config)

            assert (os.path.abspath(mappingFileName) ==
                    os.path.abspath(remapper.mappingFileName))
            assert os.path.exists(mappingFileName)

            assert isinstance(remapper.sourceDescriptor,
                              MpasMeshDescriptor)
            assert isinstance(remapper.destinationDescriptor,
                              LatLonGridDescriptor)

            if not setName:
                # copy the mapping file so it exists in the 'maps' dir
                shutil.copyfile(defaultMappingFileName,
                                explicitMappingFileName)

    def test_get_observations_remapper(self):
        config = self.setup_config()
        fieldName = 'sst'

        explicitMappingPath = '{}/maps'.format(self.test_dir)
        os.makedirs(explicitMappingPath)

        fileBase = 'map_obs_sst_1.0x1.0degree_to_0.5x0.5degree_bilinear.nc'

        defaultMappingFileName = '{}/{}'.format(self.test_dir, fileBase)
        explicitMappingFileName = '{}/{}'.format(explicitMappingPath, fileBase)

        for mappingFileName, setName in [(defaultMappingFileName, False),
                                         (explicitMappingFileName, True)]:

            remapper = self.setup_obs_remapper(config, fieldName)

            assert (os.path.abspath(mappingFileName) ==
                    os.path.abspath(remapper.mappingFileName))
            assert os.path.exists(mappingFileName)

            assert isinstance(remapper.sourceDescriptor,
                              LatLonGridDescriptor)
            assert isinstance(remapper.destinationDescriptor,
                              LatLonGridDescriptor)

            if not setName:
                # copy the mapping file so it exists in the 'maps' dir
                shutil.copyfile(defaultMappingFileName,
                                explicitMappingFileName)

    def test_compute_climatology(self):
        config = self.setup_config()
        calendar = 'gregorian_noleap'
        ds = self.open_test_ds(config, calendar)

        assert('month' not in ds.coords.keys())
        assert('daysInMonth' not in ds.coords.keys())

        # test add_months_and_days_in_month
        ds = add_years_months_days_in_month(ds, calendar)

        self.assertArrayEqual(ds.month.values, [1, 2, 3])
        self.assertArrayEqual(numpy.round(ds.daysInMonth.values), [31, 28, 31])

        # test compute_climatology on a data set
        monthNames = 'JFM'
        monthValues = constants.monthDictionary[monthNames]
        dsClimatology = compute_climatology(ds, monthValues, calendar)

        assert('Time' not in dsClimatology.dims.keys())

        self.assertEqual(list(dsClimatology.data_vars.keys()), ['mld'])

        climFileName = '{}/refSeasonalClim.nc'.format(self.datadir)
        refClimatology = xarray.open_dataset(climFileName)
        self.assertArrayApproxEqual(dsClimatology.mld.values,
                                    refClimatology.mld.values)

        # test compute_climatology on a data array
        mldClimatology = compute_climatology(ds.mld, monthValues, calendar)

        assert('Time' not in mldClimatology.dims)

        self.assertArrayApproxEqual(dsClimatology.mld.values,
                                    mldClimatology.values)

        # for good measure...
        self.assertArrayApproxEqual(mldClimatology.values,
                                    refClimatology.mld.values)

    def test_compute_monthly_climatology(self):
        config = self.setup_config()
        calendar = 'gregorian_noleap'
        ds = self.open_test_ds(config, calendar)

        monthlyClimatology = compute_monthly_climatology(ds, calendar)

        assert(len(monthlyClimatology.month) == 3)

        self.assertEqual(list(monthlyClimatology.data_vars.keys()), ['mld'])

        climFileName = '{}/refMonthlyClim.nc'.format(self.datadir)
        refClimatology = xarray.open_dataset(climFileName)

        self.assertArrayApproxEqual(monthlyClimatology.mld.values,
                                    refClimatology.mld.values)

        self.assertArrayApproxEqual(monthlyClimatology.month.values,
                                    refClimatology.month.values)


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
