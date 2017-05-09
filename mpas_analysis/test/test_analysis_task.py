"""
Unit tests for utility functions in run_analysis

Xylar Asay-Davis
"""

import pytest
import shutil
import tempfile
import os

from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.shared.analysis_task import AnalysisTask
from mpas_analysis.configuration.MpasAnalysisConfigParser \
    import MpasAnalysisConfigParser


@pytest.mark.usefixtures("loaddatadir")
class TestAnalysisTask(TestCase):
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def setup_config(self):
        config = MpasAnalysisConfigParser()
        config.read('config.default')
        config.set('input', 'baseDirectory', str(self.datadir))
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

    def test_setup_and_check(self):
        config = self.setup_config()
        task = self.setup_task(config)
        # make sure everything was set up as expected
        self.assertEqual(task.calendar, 'gregorian_noleap')
        self.assertEqual(task.runDirectory, '{}/.'.format(self.datadir))
        self.assertEqual(task.historyDirectory, '{}/.'.format(self.datadir))
        self.assertEqual(task.plotsDirectory, '{}/plots'.format(self.test_dir))
        assert(os.path.exists(task.plotsDirectory))
        assert(task.namelistMap is not None)
        assert(task.streamMap is not None)
        assert(task.variableMap is not None)
        assert(config.has_option('climatology', 'startDate'))
        assert(config.has_option('climatology', 'endDate'))

    def test_check_generate(self):

        def doTest(generate, expectedResults):
            config = MpasAnalysisConfigParser()
            config.add_section('output')
            config.set('output', 'generate', generate)
            for taskName in expectedResults:
                genericTask = AnalysisTask(config=config,
                                           taskName=taskName,
                                           componentName=cores[taskName],
                                           tags=tags[taskName])
                expectedResult = expectedResults[taskName]
                result = genericTask.check_generate()
                self.assertEqual(result, expectedResult)

        # Comments from config.template about how generate works:
        #
        # a list of analyses to generate.  Valid names are:
        #   'timeSeriesOHC', 'timeSeriesSST', 'regriddedSST',
        #   'regriddedSSS', 'regriddedMLD', 'timeSeriesSeaIceAreaVol',
        #   'regriddedSeaIceConcThick'
        # the following shortcuts exist:
        #   'all' -- all analyses will be run
        #   'all_timeSeries' -- all time-series analyses will be run
        #   'all_regriddedHorizontal' -- all analyses involving regridded
        #                                horizontal fields will be run
        #   'all_ocean' -- all ocean analyses will be run
        #   'all_seaIce' -- all sea-ice analyses will be run
        #   'no_timeSeriesOHC' -- skip 'timeSeriesOHC' (and similarly with the
        #                             other analyses).
        #   'no_ocean', 'no_timeSeries', etc. -- in analogy to 'all_*', skip
        #                                        the given category of analysis
        # an equivalent syntax can be used on the command line to override this
        # option:
        #    ./run_analysis.py config.analysis --generate \
        #         all,no_ocean,all_timeSeries

        cores = {'timeSeriesOHC': 'ocean',
                 'timeSeriesSST': 'ocean',
                 'indexNino34': 'ocean',
                 'meridionalHeatTransport': 'ocean',
                 'streamfunctionMOC': 'ocean',
                 'climatologyMapSST': 'ocean',
                 'climatologyMapMLD': 'ocean',
                 'climatologyMapSSS': 'ocean',
                 'timeSeriesSeaIceAreaVol': 'seaIce',
                 'climatologyMapSeaIceConcThick': 'seaIce'}

        tags = {'timeSeriesOHC': ['timeSeries', 'ohc'],
                'timeSeriesSST': ['timeSeries', 'sst'],
                'indexNino34': ['index', 'nino'],
                'meridionalHeatTransport': ['climatology', 'mht'],
                'streamfunctionMOC': ['climatology', 'timeSeries',
                                      'streamfunction', 'moc'],
                'climatologyMapSST': ['climatology', 'horizontalMap', 'sst'],
                'climatologyMapMLD': ['climatology', 'horizontalMap', 'mld'],
                'climatologyMapSSS': ['climatology', 'horizontalMap', 'sss'],
                'timeSeriesSeaIceAreaVol': ['timeSeries'],
                'climatologyMapSeaIceConcThick': ['climatology',
                                                  'horizontalMap']}

        # test 'all'
        expectedResults = {}
        for taskName in cores:
            expectedResults[taskName] = True
        doTest("['all']", expectedResults)

        # test 'all_<category>' and ['all', 'no_<category>']
        allTags = []
        for taskName in tags:
            allTags.extend(tags[taskName])

        for tag in set(allTags):
            expectedResults = {}
            for taskName in tags:
                expectedResults[taskName] = \
                    (tag in tags[taskName])
            doTest("['all_{}']".format(tag), expectedResults)

            expectedResults = {}
            for taskName in tags:
                expectedResults[taskName] = \
                    (tag not in tags[taskName])
            doTest("['all', 'no_{}']".format(tag), expectedResults)

        # test 'all_<core>' and ['all', 'no_<core>']
        for core in set(cores.values()):
            expectedResults = {}
            for taskName in cores:
                expectedResults[taskName] = \
                    (cores[taskName] == core)
            doTest("['all_{}']".format(core), expectedResults)

            expectedResults = {}
            for taskName in cores:
                expectedResults[taskName] = \
                    (cores[taskName] != core)
            doTest("['all','no_{}']".format(core), expectedResults)

        # test each analysis individually
        for taskName in cores:
            expectedResults = {}
            for otherAnalysis in cores:
                expectedResults[otherAnalysis] = \
                    (taskName == otherAnalysis)
            doTest("['{}']".format(taskName), expectedResults)

        # test a non-existent analysis
        expectedResults = {}
        for taskName in cores:
            expectedResults[taskName] = False
        doTest("['fakeAnalysis']", expectedResults)

        # test ['all', 'no_ocean', 'all_timeSeries']
        expectedResults = {}
        for taskName in cores:
            expectedResults[taskName] = True
        for taskName in cores:
            if cores[taskName] == 'ocean':
                expectedResults[taskName] = False
        for taskName in tags:
            if 'timeSeries' in tags[taskName]:
                expectedResults[taskName] = True
        doTest("['all', 'no_ocean', 'all_timeSeries']", expectedResults)

        # test ['all', 'no_timeSeriesOHC']
        expectedResults = {}
        for taskName in cores:
            expectedResults[taskName] = True
        expectedResults['timeSeriesOHC'] = False
        doTest("['all', 'no_timeSeriesOHC']", expectedResults)

    def test_update_start_end_date(self):
        config = self.setup_config()
        task = self.setup_task(config)

        inputFileNames = \
            task.get_input_file_names(streamName='timeSeriesStats',
                                      startAndEndDateSection='climatology')

        timeCache = task.cache_multifile_dataset_times(
            inputFileNames, streamName='timeSeriesStats',
            timeVariableName='Time')

        # make sure the times have been cached
        assert(os.path.exists(
            '{}/timecache/ocean_timeSeriesStats_times.pickle'.format(
                self.test_dir)))

        for index, fileName in enumerate(timeCache.keys()):
            assert(fileName == os.path.abspath(inputFileNames[index]))
            assert(timeCache[fileName]['years'][0] == 2)
            assert(timeCache[fileName]['months'][0] == index+1)

        changed = task.update_start_end_date(section='climatology',
                                             streamName='timeSeriesStats')

        assert(not changed)
        assert(config.getint('climatology', 'startYear') == 2)
        assert(config.getint('climatology', 'endYear') == 2)
        assert(config.get('climatology', 'startDate') == '0002-01-01_00:00:00')
        assert(config.get('climatology', 'endDate') == '0002-12-31_23:59:59')

        config.set('climatology', 'endYear', '50')

        with self.assertWarns('climatology start and/or end year different '
                              'from requested'):
            changed = task.update_start_end_date(section='climatology',
                                                 streamName='timeSeriesStats')

        assert(changed)
        assert(config.getint('climatology', 'startYear') == 2)
        assert(config.getint('climatology', 'endYear') == 2)
        assert(config.get('climatology', 'startDate') == '0002-01-01_00:00:00')
        assert(config.get('climatology', 'endDate') == '0002-12-31_23:59:59')

    def test_get_input_file_names(self):
        config = self.setup_config()
        task = self.setup_task(config)

        inputFileNames = \
            task.get_input_file_names(streamName='timeSeriesStats',
                                      startDate='0002-01-01_00:00:00',
                                      endDate='0002-03-01_00:00:00')

        for index, month in enumerate([1, 2]):
            expectedFileName = \
                '{}/./timeSeries.0002-{:02d}-01.nc'.format(str(self.datadir),
                                                           month)
            assert(inputFileNames[index] == expectedFileName)

        inputFileNames = \
            task.get_input_file_names(streamName='timeSeriesStats',
                                      startAndEndDateSection='climatology')

        for index, month in enumerate([1, 2, 3]):
            expectedFileName = \
                '{}/./timeSeries.0002-{:02d}-01.nc'.format(str(self.datadir),
                                                           month)
            assert(inputFileNames[index] == expectedFileName)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
