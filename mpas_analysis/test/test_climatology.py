"""
Unit test infrastructure for climatologies.

Xylar Asay-Davis
04/03/2017
"""

import pytest
import tempfile
import shutil
import os
import numpy
import xarray

from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.shared.generalized_reader.generalized_reader \
    import open_multifile_dataset
from mpas_analysis.configuration.MpasAnalysisConfigParser \
    import MpasAnalysisConfigParser
from mpas_analysis.shared.climatology import climatology
from mpas_analysis.shared.constants import constants

@pytest.mark.usefixtures("loaddatadir")
class TestClimatology(TestCase):

    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def setup_config(self, autocloseFileLimitFraction=0.5,
                     maxChunkSize=10000):
        config = MpasAnalysisConfigParser()
        config.add_section('input')
        config.set('input', 'autocloseFileLimitFraction',
                   str(autocloseFileLimitFraction))
        config.set('input', 'maxChunkSize', str(maxChunkSize))
        config.set('input', 'mpasMeshName', 'QU240')

        config.add_section('output')
        config.set('output', 'baseDirectory', self.test_dir)
        config.set('output', 'mappingSubdirectory', '.')
        config.set('output', 'mpasClimatologySubdirectory', 'clim/mpas')
        config.set('output', 'mpasRegriddedClimSubdirectory',
                   'clim/mpas/regrid')

        config.add_section('climatology')
        config.set('climatology', 'startYear', '2')
        config.set('climatology', 'endYear', '2')
        config.set('climatology', 'comparisonLatResolution', '0.5')
        config.set('climatology', 'comparisonLonResolution', '0.5')

        config.set('climatology', 'overwriteMapping', 'False')
        config.set('climatology', 'overwriteMpasClimatology', 'False')
        config.set('climatology', 'mpasInterpolationMethod', 'bilinear')

        config.add_section('oceanObservations')
        config.set('oceanObservations', 'interpolationMethod', 'bilinear')
        config.set('oceanObservations', 'climatologySubdirectory', 'clim/obs')
        config.set('oceanObservations', 'regriddedClimSubdirectory',
                   'clim/obs/regrid')

        return config

    def test_write_mpas_mapping_file(self):
        config = self.setup_config()
        mpasMeshFileName = '{}/mpasMesh.nc'.format(self.datadir)
        climatology.write_mpas_mapping_file(config, mpasMeshFileName)

        mappingFileName = '{}/map_QU240_to_0.5x0.5degree_' \
                          'bilinear.nc'.format(self.test_dir)
        assert os.path.exists(mappingFileName)

        mappingFileName = '{}/mapping.nc'.format(self.test_dir)
        config.set('climatology', 'mpasMappingFile', mappingFileName)
        climatology.write_mpas_mapping_file(config, mpasMeshFileName)
        assert os.path.exists(mappingFileName)

    def test_write_observations_mapping_file(self):
        config = self.setup_config()
        gridFileName = '{}/obsGrid.nc'.format(self.datadir)
        componentName = 'ocean'
        fieldName = 'sst'
        climatology.write_observations_mapping_file(config,
                                                    componentName,
                                                    fieldName,
                                                    gridFileName,
                                                    latVarName='lat',
                                                    lonVarName='lon')

        mappingFileName = '{}/map_obs_{}_1.0x1.0degree_to_0.5x0.5degree_' \
                          'bilinear.nc'.format(self.test_dir, fieldName)
        assert os.path.exists(mappingFileName)

        mappingFileName = '{}/mapping.nc'.format(self.test_dir)
        config.set('oceanObservations', 'sstClimatologyMappingFile',
                   mappingFileName)
        climatology.write_observations_mapping_file(config,
                                                    componentName,
                                                    fieldName,
                                                    gridFileName,
                                                    latVarName='lat',
                                                    lonVarName='lon')
        assert os.path.exists(mappingFileName)

    def test_get_mpas_climatology_file_names(self):
        config = self.setup_config()
        fieldName = 'sst'
        monthNames = 'JFM'
        (climatologyFileName, regriddedFileName) = \
            climatology.get_mpas_climatology_file_names(config, fieldName,
                                                        monthNames)
        expectedClimatologyFileName = '{}/clim/mpas/sst_QU240_JFM_' \
                                      'years0002-0002.nc'.format(self.test_dir)
        self.assertEqual(climatologyFileName, expectedClimatologyFileName)

        expectedRegriddedFileName = '{}/clim/mpas/regrid/sst_QU240_to_' \
                                    '0.5x0.5degree_JFM_' \
                                    'years0002-0002.nc'.format(self.test_dir)
        self.assertEqual(regriddedFileName, expectedRegriddedFileName)

    def test_get_observation_climatology_file_names(self):
        config = self.setup_config()
        fieldName = 'sst'
        monthNames = 'JFM'
        gridFileName = '{}/obsGrid.nc'.format(self.datadir)
        componentName = 'ocean'
        (climatologyFileName, regriddedFileName) = \
            climatology.get_observation_climatology_file_names(
                config, fieldName, monthNames, componentName, gridFileName,
                latVarName='lat', lonVarName='lon')
        expectedClimatologyFileName = '{}/clim/obs/sst_1.0x1.0degree_' \
                                      'JFM.nc'.format(self.test_dir)
        self.assertEqual(climatologyFileName, expectedClimatologyFileName)

        expectedRegriddedFileName = '{}/clim/obs/regrid/sst_1.0x1.0degree_' \
                                    'to_0.5x0.5degree_' \
                                    'JFM.nc'.format(self.test_dir)
        self.assertEqual(regriddedFileName, expectedRegriddedFileName)

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

    def test_compute_climatology(self):
        config = self.setup_config()
        calendar = 'gregorian_noleap'
        ds = self.open_test_ds(config, calendar)

        assert('month' not in ds.coords.keys())
        assert('daysInMonth' not in ds.coords.keys())

        # test add_months_and_days_in_month
        ds = climatology.add_months_and_days_in_month(ds, calendar)

        self.assertArrayEqual(ds.month.values, [1, 2, 3])
        self.assertArrayEqual(numpy.round(ds.daysInMonth.values), [31, 28, 31])

        # test compute_climatology on a data set
        monthNames = 'JFM'
        monthValues = constants.monthDictionary[monthNames]
        dsClimatology = climatology.compute_climatology(ds, monthValues,
                                                        calendar)

        assert('Time' not in dsClimatology.dims.keys())

        self.assertEqual(dsClimatology.data_vars.keys(), ['mld'])

        climFileName = '{}/refSeasonalClim.nc'.format(self.datadir)
        refClimatology = xarray.open_dataset(climFileName)
        self.assertArrayApproxEqual(dsClimatology.mld.values,
                                    refClimatology.mld.values)

        # test compute_climatology on a data array
        mldClimatology = climatology.compute_climatology(ds.mld, monthValues,
                                                         calendar)

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

        monthlyClimatology = climatology.compute_monthly_climatology(ds,
                                                                     calendar)

        assert(len(monthlyClimatology.month) == 3)

        self.assertEqual(monthlyClimatology.data_vars.keys(), ['mld'])

        climFileName = '{}/refMonthlyClim.nc'.format(self.datadir)
        refClimatology = xarray.open_dataset(climFileName)

        self.assertArrayApproxEqual(monthlyClimatology.mld.values,
                                    refClimatology.mld.values)

        self.assertArrayApproxEqual(monthlyClimatology.month.values,
                                    refClimatology.month.values)

    def test_update_start_end_year(self):
        config = self.setup_config()
        calendar = 'gregorian_noleap'
        ds = self.open_test_ds(config, calendar)

        changed, startYear, endYear = \
            climatology.update_start_end_year(ds, config, calendar)

        assert(not changed)
        assert(startYear == 2)
        assert(endYear == 2)

        config.set('climatology', 'endYear', '50')
        ds = self.open_test_ds(config, calendar)

        with self.assertWarns('climatology start and/or end year different from requested'):
            changed, startYear, endYear = \
                climatology.update_start_end_year(ds, config, calendar)

        assert(changed)
        assert(startYear == 2)
        assert(endYear == 2)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
