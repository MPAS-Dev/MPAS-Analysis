"""
Unit test infrastructure for climatologies.

Xylar Asay-Davis
04/11/2017
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
from mpas_analysis.shared.climatology import \
    get_lat_lon_comparison_descriptor, get_remapper, \
    get_mpas_climatology_file_names, get_observation_climatology_file_names, \
    add_years_months_days_in_month, compute_climatology, \
    compute_monthly_climatology, update_start_end_year, cache_climatologies
from mpas_analysis.shared.grid import MpasMeshDescriptor, LatLonGridDescriptor
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
        config.set('output', 'mpasRemappedClimSubdirectory',
                   'clim/mpas/remap')

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
            get_lat_lon_comparison_descriptor(config)

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
            get_lat_lon_comparison_descriptor(config)

        obsDescriptor = LatLonGridDescriptor()
        obsDescriptor.read(fileName=gridFileName, latVarName='lat',
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

        for mappingFileName, setName in [(defaultMappingFileName,  False),
                                         (explicitMappingFileName, True)]:
            if setName:
                config.set('input', 'mappingDirectory', explicitMappingPath)

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

        for mappingFileName, setName in [(defaultMappingFileName,  False),
                                         (explicitMappingFileName, True)]:

            if setName:
                config.set('input', 'mappingDirectory', explicitMappingPath)

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

    def test_get_mpas_climatology_file_names(self):
        config = self.setup_config()
        fieldName = 'sst'
        monthNames = 'JFM'

        remapper = self.setup_mpas_remapper(config)

        (climatologyFileName, climatologyPrefix, remappedFileName) = \
            get_mpas_climatology_file_names(
                config, fieldName, monthNames,
                remapper.sourceDescriptor.meshName,
                remapper.destinationDescriptor.meshName)
        expectedClimatologyFileName = '{}/clim/mpas/sst_QU240_JFM_' \
                                      'year0002.nc'.format(self.test_dir)
        self.assertEqual(climatologyFileName, expectedClimatologyFileName)

        expectedClimatologyPrefix = '{}/clim/mpas/sst_QU240_' \
                                    'JFM'.format(self.test_dir)
        self.assertEqual(climatologyPrefix, expectedClimatologyPrefix)

        expectedRemappedFileName = '{}/clim/mpas/remap/sst_QU240_to_' \
                                   '0.5x0.5degree_JFM_' \
                                   'year0002.nc'.format(self.test_dir)
        self.assertEqual(remappedFileName, expectedRemappedFileName)

    def test_get_observation_climatology_file_names(self):
        config = self.setup_config()
        fieldName = 'sst'
        monthNames = 'JFM'
        componentName = 'ocean'

        remapper = self.setup_obs_remapper(config, fieldName)

        (climatologyFileName, remappedFileName) = \
            get_observation_climatology_file_names(
                config, fieldName, monthNames, componentName, remapper)
        expectedClimatologyFileName = '{}/clim/obs/sst_1.0x1.0degree_' \
                                      'JFM.nc'.format(self.test_dir)
        self.assertEqual(climatologyFileName, expectedClimatologyFileName)

        expectedRemappedFileName = '{}/clim/obs/remap/sst_1.0x1.0degree_' \
                                   'to_0.5x0.5degree_' \
                                   'JFM.nc'.format(self.test_dir)
        self.assertEqual(remappedFileName, expectedRemappedFileName)

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

        self.assertEqual(dsClimatology.data_vars.keys(), ['mld'])

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
            update_start_end_year(ds, config, calendar)

        assert(not changed)
        assert(startYear == 2)
        assert(endYear == 2)

        config.set('climatology', 'endYear', '50')
        ds = self.open_test_ds(config, calendar)

        with self.assertWarns('climatology start and/or end year different '
                              'from requested'):
            changed, startYear, endYear = \
                update_start_end_year(ds, config, calendar)

        assert(changed)
        assert(startYear == 2)
        assert(endYear == 2)

    def cache_climatologies_setup(self):
        config = self.setup_config()
        calendar = 'gregorian_noleap'
        ds = self.open_test_ds(config, calendar)
        fieldName = 'mld'
        climFileName = '{}/refSeasonalClim.nc'.format(self.datadir)
        refClimatology = xarray.open_dataset(climFileName)

        remapper = self.setup_mpas_remapper(config)

        return {'config': config, 'calendar': calendar, 'ds': ds,
                'fieldName': fieldName, 'climFileName': climFileName,
                'refClimatology': refClimatology, 'remapper': remapper}

    def test_jan_1yr_climo_test1(self):
        setup = self.cache_climatologies_setup()
        # test1: Just January, 1-year climatologies are cached; only one file
        #        is produced with suffix year0002; a second run of
        #        cache_climatologies doesn't modify any files
        test1 = {'monthNames': 'Jan',
                 'monthValues': [1],
                 'yearsPerCacheFile': 1,
                 'expectedSuffixes': ['year0002'],
                 'expectedModified': [False],
                 # weird value because first time step of Jan. missing in ds
                 'expectedDays': 30.958333,
                 'expectedMonths': 1,
                 'refClimatology': None}
        self.cache_climatologies_driver(test1, **setup)

    def test_jfm_1yr_climo_test2(self):
        setup = self.cache_climatologies_setup()
        # same as test1 but with JFM
        test2 = {'monthNames': 'JFM',
                 'monthValues': constants.monthDictionary['JFM'],
                 'yearsPerCacheFile': 1,
                 'expectedSuffixes': ['year0002'],
                 'expectedModified': [False],
                 # weird value because first time step of Jan. missing in ds
                 'expectedDays': 89.958333,
                 'expectedMonths': 3,
                 'refClimatology': setup['refClimatology']}
        self.cache_climatologies_driver(test2, **setup)

    def test_jan_2yr_climo_test3(self):
        setup = self.cache_climatologies_setup()
        # test3: 2-year climatologies are cached; 2 files are produced
        #        with suffix years0002-0003 (the "individual" climatology
        #        file) and year0002 (the "aggregated" climatology file);
        #        a second tries to update the "individual" cache file
        #        because it appears to be incomplete but does not attempt
        #        to update the aggregated climatology file because no
        #        additional years were processed and the file was already
        #        complete for the span of years present
        test3 = {'monthNames': 'Jan',
                 'monthValues': [1],
                 'yearsPerCacheFile': 2,
                 'expectedSuffixes': ['years0002-0003', 'year0002'],
                 'expectedModified': [True, False],
                 # weird value because first time step of Jan. missing in ds
                 'expectedDays': 30.958333,
                 'expectedMonths': 1,
                 'refClimatology': None}
        self.cache_climatologies_driver(test3, **setup)

    def test_jfm_2yr_climo_test4(self):
        setup = self.cache_climatologies_setup()
        # test4: same as test3 but with JFM
        test4 = {'monthNames': 'JFM',
                 'monthValues': constants.monthDictionary['JFM'],
                 'yearsPerCacheFile': 2,
                 'expectedSuffixes': ['years0002-0003', 'year0002'],
                 'expectedModified': [True, False],
                 # weird value because first time step of Jan. missing in ds
                 'expectedDays': 89.958333,
                 'expectedMonths': 3,
                 'refClimatology': setup['refClimatology']}
        self.cache_climatologies_driver(test4, **setup)

    def cache_climatologies_driver(self, test, config, fieldName,
                                   ds, remapper, calendar, **kwargs):
        monthNames = test['monthNames']
        monthValues = test['monthValues']
        yearsPerCacheFile = test['yearsPerCacheFile']
        expectedSuffixes = test['expectedSuffixes']
        expectedModified = test['expectedModified']
        expectedDays = test['expectedDays']
        expectedMonths = test['expectedMonths']
        refClimatology = test['refClimatology']

        (climatologyFileName, climatologyPrefix) = \
            get_mpas_climatology_file_names(
                config, fieldName, monthNames,
                remapper.sourceDescriptor.meshName)

        config.set('climatology', 'yearsPerCacheFile',
                   str(yearsPerCacheFile))
        # once without cache files
        dsClimatology = cache_climatologies(
            ds, monthValues, config, climatologyPrefix, calendar,
            printProgress=True)

        if refClimatology is not None:
            self.assertArrayApproxEqual(dsClimatology.mld.values,
                                        refClimatology.mld.values)

        self.assertEqual(dsClimatology.attrs['totalMonths'],
                         expectedMonths)
        self.assertApproxEqual(dsClimatology.attrs['totalDays'],
                               expectedDays)
        dsClimatology.close()

        fingerprints = []
        for suffix in expectedSuffixes:
            expectedClimatologyFileName = '{}/clim/mpas/mld_QU240_' \
                                          '{}_{}.nc'.format(
                                              self.test_dir, monthNames,
                                              suffix)
            assert os.path.exists(expectedClimatologyFileName)

            dsClimatology = xarray.open_dataset(expectedClimatologyFileName)
            fingerprints.append(dsClimatology.fingerprintClimo)

        # try it again with cache files saved
        dsClimatology = cache_climatologies(
            ds, monthValues, config,  climatologyPrefix, calendar,
            printProgress=True)

        if refClimatology is not None:
            self.assertArrayApproxEqual(dsClimatology.mld.values,
                                        refClimatology.mld.values)

        self.assertEqual(dsClimatology.attrs['totalMonths'],
                         expectedMonths)
        self.assertApproxEqual(dsClimatology.attrs['totalDays'],
                               expectedDays)
        dsClimatology.close()

        for index, suffix in enumerate(expectedSuffixes):
            expectedClimatologyFileName = '{}/clim/mpas/mld_QU240_' \
                                          '{}_{}.nc'.format(
                                              self.test_dir, monthNames,
                                              suffix)

            dsClimatology = xarray.open_dataset(expectedClimatologyFileName)
            fingerprintCheck = dsClimatology.fingerprintClimo

            # Check whether the given file was modified, and whether
            # this was the expected result
            fileWasModified = fingerprints[index] != fingerprintCheck
            assert fileWasModified == expectedModified[index]

            # remove the cache file for the next try
            os.remove(expectedClimatologyFileName)


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
