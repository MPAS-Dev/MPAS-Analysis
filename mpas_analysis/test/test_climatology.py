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
from functools import partial

from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.shared.generalized_reader \
    import open_multifile_dataset
from mpas_analysis.configuration.MpasAnalysisConfigParser \
    import MpasAnalysisConfigParser
from mpas_analysis.shared.climatology import Climatology, \
    MpasClimatology, ObservationClimatology
from mpas_analysis.shared.grid import MpasMeshDescriptor, LatLonGridDescriptor
from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.timekeeping.utility import \
    add_years_months_days_in_month

from mpas_analysis.shared.analysis_task import AnalysisTask


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
        config.read('config.default')
        config.set('input', 'baseDirectory', str(self.datadir))
        config.set('input', 'autocloseFileLimitFraction',
                   str(autocloseFileLimitFraction))
        config.set('input', 'maxChunkSize', str(maxChunkSize))
        config.set('input', 'mpasMeshName', 'QU240')

        config.set('output', 'baseDirectory', self.test_dir)
        config.set('output', 'mappingSubdirectory', '.')

        config.set('climatology', 'startYear', '2')
        config.set('climatology', 'endYear', '2')

        return config

    def setup_task(self, config):
        task = AnalysisTask(config=config,
                            taskName='genericClimatology',
                            componentName='ocean',
                            tags=['climatology', 'horizontalMap'])
        task.setup_and_check()
        return task

    def setup_mpas_climatology(self, config, task):

        fieldName = 'mld'
        monthNames = 'JFM'
        streamName = 'timeSeriesStats'

        mpasMeshFileName = '{}/mpasMesh.nc'.format(self.datadir)

        climatology = MpasClimatology(task=task,
                                      fieldName=fieldName,
                                      monthNames=monthNames,
                                      streamName=streamName,
                                      meshFileName=mpasMeshFileName,
                                      comparisonGrid='latlon',
                                      mappingFileSection='climatology',
                                      mappingFileOption='mpasMappingFile',
                                      mappingFilePrefix='map',
                                      method=config.get(
                                          'climatology',
                                          'mpasInterpolationMethod'))

        return climatology

    def setup_obs_climatology(self, config, task):
        gridFileName = '{}/obsGrid.nc'.format(self.datadir)
        fieldName = 'sst'
        monthNames = 'JFM'

        obsDescriptor = LatLonGridDescriptor.read(fileName=gridFileName,
                                                  latVarName='lat',
                                                  lonVarName='lon')

        climatology = \
            ObservationClimatology(
                task=task,
                fieldName=fieldName,
                monthNames=monthNames,
                obsGridDescriptor=obsDescriptor,
                comparisonGrid='latlon',
                mappingFileSection='oceanObservations',
                mappingFileOption='sstClimatologyMappingFile',
                mappingFilePrefix='map_obs_{}'.format(fieldName),
                method=config.get('oceanObservations',
                                  'interpolationMethod'))

        return climatology

    def open_ds_part(self, task, inputFileNames, startDate, endDate):
        variableMap = {'mld': ['timeMonthly_avg_tThreshMLD'],
                       'Time': [['xtime_startMonthly', 'xtime_endMonthly']]}

        variableList = ['mld']

        ds = open_multifile_dataset(
            fileNames=inputFileNames,
            calendar=task.calendar,
            config=task.config,
            timeVariableName='Time',
            variableList=variableList,
            variableMap=variableMap,
            startDate=startDate,
            endDate=endDate)

        return ds

    def open_test_ds(self, task):
        fileNames = ['{}/timeSeries.0002-{:02d}-01.nc'.format(self.datadir,
                                                              month)
                     for month in [1, 2, 3]]
        ds = self.open_ds_part(task, fileNames, None, None)
        assert(len(ds.Time) == 3)
        return ds

    def test_mpas_remapping(self):
        config = self.setup_config()
        task = self.setup_task(config)

        defaultMappingFileName = '{}/map_QU240_to_0.5x0.5degree_' \
                                 'bilinear.nc'.format(self.test_dir)

        explicitMappingFileName = '{}/mapping.nc'.format(self.test_dir)

        for mappingFileName, setName in [(defaultMappingFileName,  False),
                                         (explicitMappingFileName, True)]:
            if setName:
                config.set('climatology', 'mpasMappingFile', mappingFileName)

            climatology = self.setup_mpas_climatology(config, task)

            remapper = climatology.remapper

            assert (os.path.abspath(mappingFileName) ==
                    os.path.abspath(remapper.mappingFileName))
            assert os.path.exists(mappingFileName)

            assert isinstance(remapper.sourceDescriptor,
                              MpasMeshDescriptor)
            assert isinstance(remapper.destinationDescriptor,
                              LatLonGridDescriptor)

    def test_observation_remapping(self):
        config = self.setup_config()
        task = self.setup_task(config)

        defaultMappingFileName = '{}/map_obs_sst_1.0x1.0degree_to_' \
                                 '0.5x0.5degree_bilinear.nc'.format(
                                     self.test_dir)

        explicitMappingFileName = '{}/mapping.nc'.format(self.test_dir)

        for mappingFileName, setName in [(defaultMappingFileName,  False),
                                         (explicitMappingFileName, True)]:

            if setName:
                config.set('oceanObservations', 'sstClimatologyMappingFile',
                           mappingFileName)

            climatology = self.setup_obs_climatology(config, task)

            remapper = climatology.remapper

            assert (os.path.abspath(mappingFileName) ==
                    os.path.abspath(remapper.mappingFileName))
            assert os.path.exists(mappingFileName)

            assert isinstance(remapper.sourceDescriptor,
                              LatLonGridDescriptor)
            assert isinstance(remapper.destinationDescriptor,
                              LatLonGridDescriptor)

    def test_mpas_climatology_file_names(self):
        config = self.setup_config()
        task = self.setup_task(config)

        climatology = self.setup_mpas_climatology(config, task)

        expectedClimatologyFileName = '{}/clim/mpas/mld_QU240_JFM_' \
                                      'year0002.nc'.format(self.test_dir)
        self.assertEqual(climatology.climatologyFileName,
                         expectedClimatologyFileName)

        expectedClimatologyPrefix = '{}/clim/mpas/mld_QU240_' \
                                    'JFM'.format(self.test_dir)
        self.assertEqual(climatology.climatologyPrefix,
                         expectedClimatologyPrefix)

        expectedRemappedFileName = '{}/clim/mpas/remapped/mld_QU240_to_' \
                                   '0.5x0.5degree_JFM_' \
                                   'year0002.nc'.format(self.test_dir)
        self.assertEqual(climatology.remappedFileName,
                         expectedRemappedFileName)

    def test_observation_climatology_file_names(self):
        config = self.setup_config()
        task = self.setup_task(config)

        climatology = self.setup_obs_climatology(config, task)

        expectedClimatologyFileName = '{}/clim/obs/sst_1.0x1.0degree_' \
                                      'JFM.nc'.format(self.test_dir)
        self.assertEqual(climatology.climatologyFileName,
                         expectedClimatologyFileName)

        expectedRemappedFileName = '{}/clim/obs/remapped/sst_1.0x1.0degree_' \
                                   'to_0.5x0.5degree_' \
                                   'JFM.nc'.format(self.test_dir)
        self.assertEqual(climatology.remappedFileName,
                         expectedRemappedFileName)

    def test_climatology_compute(self):
        config = self.setup_config()
        task = self.setup_task(config)
        ds = self.open_test_ds(task)

        assert('month' not in ds.coords.keys())
        assert('daysInMonth' not in ds.coords.keys())

        # test add_months_and_days_in_month
        ds = add_years_months_days_in_month(ds, task.calendar)

        self.assertArrayEqual(ds.month.values, [1, 2, 3])
        self.assertArrayEqual(numpy.round(ds.daysInMonth.values), [31, 28, 31])

        climatology = self.setup_mpas_climatology(config, task)
        dsClimatology = climatology.compute(ds, maskVaries=False)

        assert('Time' not in dsClimatology.dims.keys())

        self.assertEqual(dsClimatology.data_vars.keys(), ['mld'])

        climFileName = '{}/refSeasonalClim.nc'.format(self.datadir)
        refClimatology = xarray.open_dataset(climFileName)
        self.assertArrayApproxEqual(dsClimatology.mld.values,
                                    refClimatology.mld.values)

        # test compute_climatology on a data array
        mldClimatology = climatology.compute(ds.mld, maskVaries=False)

        assert('Time' not in mldClimatology.dims)

        self.assertArrayApproxEqual(dsClimatology.mld.values,
                                    mldClimatology.values)

        # for good measure...
        self.assertArrayApproxEqual(mldClimatology.values,
                                    refClimatology.mld.values)

    def test_compute_monthly_climatology(self):
        config = self.setup_config()
        task = self.setup_task(config)
        ds = self.open_test_ds(task)

        climatology = Climatology(task)
        dsMonthly = climatology.compute_monthly(ds, maskVaries=False)

        assert(len(dsMonthly.month) == 3)

        self.assertEqual(dsMonthly.data_vars.keys(), ['mld'])

        climFileName = '{}/refMonthlyClim.nc'.format(self.datadir)
        refClimatology = xarray.open_dataset(climFileName)

        self.assertArrayApproxEqual(dsMonthly.mld.values,
                                    refClimatology.mld.values)

        self.assertArrayApproxEqual(dsMonthly.month.values,
                                    refClimatology.month.values)

    def test_jan_1yr_climo_test1(self):
        task, refClimatology = self.cache_climatologies_setup()
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
        self.cache_climatologies_driver(test1, task)

    def test_jfm_1yr_climo_test2(self):
        task, refClimatology = self.cache_climatologies_setup()
        # same as test1 but with JFM
        test2 = {'monthNames': 'JFM',
                 'monthValues': constants.monthDictionary['JFM'],
                 'yearsPerCacheFile': 1,
                 'expectedSuffixes': ['year0002'],
                 'expectedModified': [False],
                 # weird value because first time step of Jan. missing in ds
                 'expectedDays': 89.958333,
                 'expectedMonths': 3,
                 'refClimatology': refClimatology}
        self.cache_climatologies_driver(test2, task)

    def test_jan_2yr_climo_test3(self):
        task, refClimatology = self.cache_climatologies_setup()
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
        self.cache_climatologies_driver(test3, task)

    def test_jfm_2yr_climo_test4(self):
        task, refClimatology = self.cache_climatologies_setup()
        # test4: same as test3 but with JFM
        test4 = {'monthNames': 'JFM',
                 'monthValues': constants.monthDictionary['JFM'],
                 'yearsPerCacheFile': 2,
                 'expectedSuffixes': ['years0002-0003', 'year0002'],
                 'expectedModified': [True, False],
                 # weird value because first time step of Jan. missing in ds
                 'expectedDays': 89.958333,
                 'expectedMonths': 3,
                 'refClimatology': refClimatology}
        self.cache_climatologies_driver(test4, task)

    def cache_climatologies_setup(self):
        config = self.setup_config()
        task = self.setup_task(config)
        climFileName = '{}/refSeasonalClim.nc'.format(self.datadir)
        refClimatology = xarray.open_dataset(climFileName)
        return task, refClimatology

    def cache_climatologies_driver(self, test, task):
        monthNames = test['monthNames']
        yearsPerCacheFile = test['yearsPerCacheFile']
        expectedSuffixes = test['expectedSuffixes']
        expectedModified = test['expectedModified']
        expectedDays = test['expectedDays']
        expectedMonths = test['expectedMonths']
        refClimatology = test['refClimatology']

        task.config.set('climatology', 'yearsPerCacheFile',
                        str(yearsPerCacheFile))

        fieldName = 'mld'
        streamName = 'timeSeriesStats'

        mpasMeshFileName = '{}/mpasMesh.nc'.format(self.datadir)

        climatology = MpasClimatology(task=task,
                                      fieldName=fieldName,
                                      monthNames=monthNames,
                                      streamName=streamName,
                                      meshFileName=mpasMeshFileName,
                                      comparisonGrid='latlon',
                                      mappingFileSection='climatology',
                                      mappingFileOption='mpasMappingFile',
                                      mappingFilePrefix='map',
                                      method=task.config.get(
                                          'climatology',
                                          'mpasInterpolationMethod'))

        # once without cache files
        openDataSetFunc = partial(self.open_ds_part, task)
        dsClimatology = climatology.cache(openDataSetFunc=openDataSetFunc,
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
        dsClimatology = climatology.cache(openDataSetFunc=openDataSetFunc,
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
