"""
Unit tests for utility functions in run_analysis

Xylar Asay-Davis
"""

import pytest
from mpas_analysis.test import TestCase
from mpas_analysis.shared.analysis_task import AnalysisTask
from mpas_analysis.configuration.MpasAnalysisConfigParser \
    import MpasAnalysisConfigParser


class TestAnalysisTask(TestCase):

    def test_checkGenerate(self):

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
        #   'timeSeriesOHC', 'timeSeriesSST', 'climatologyMapSST',
        #   'climatologyMapSSS', 'climatologyMapMLD', 'timeSeriesSeaIceAreaVol',
        #   'climatologyMapSeaIceConcThick'
        # the following shortcuts exist:
        #   'all' -- all analyses will be run
        #   'all_timeSeries' -- all time-series analyses will be run
        #   'all_horizontalMap' -- all analyses involving remapped
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


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
