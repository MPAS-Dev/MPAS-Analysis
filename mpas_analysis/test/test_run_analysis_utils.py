"""
Unit tests for utility functions in run_analysis

Xylar Asay-Davis
02/03/2017
"""

import pytest
from mpas_analysis.test import TestCase
from run_analysis import checkGenerate
from mpas_analysis.configuration.MpasAnalysisConfigParser \
    import MpasAnalysisConfigParser


class TestRunAnalysisUtils(TestCase):

    def test_checkGenerate(self):

        def doTest(generate, expectedResults):
            config = MpasAnalysisConfigParser()
            config.add_section('output')
            config.set('output', 'generate', generate)
            for analysisName in expectedResults:
                expectedResult = expectedResults[analysisName]
                result = checkGenerate(
                    config, analysisName=analysisName,
                    mpasCore=cores[analysisName],
                    analysisCategory=categories[analysisName])
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
                 'timeSeriesNino34': 'ocean',
                 'timeSeriesMHT': 'ocean',
                 'timeSeriesMOC': 'ocean',
                 'regriddedSST': 'ocean',
                 'regriddedMLD': 'ocean',
                 'regriddedSSS': 'ocean',
                 'timeSeriesSeaIceAreaVol': 'seaIce',
                 'regriddedSeaIceConcThick': 'seaIce'}

        categories = {'timeSeriesOHC': 'timeSeries',
                      'timeSeriesSST': 'timeSeries',
                      'timeSeriesNino34': 'timeSeries',
                      'timeSeriesMHT': 'timeSeries',
                      'timeSeriesMOC': 'timeSeries',
                      'regriddedSST': 'regriddedHorizontal',
                      'regriddedMLD': 'regriddedHorizontal',
                      'regriddedSSS': 'regriddedHorizontal',
                      'timeSeriesSeaIceAreaVol': 'timeSeries',
                      'regriddedSeaIceConcThick': 'regriddedHorizontal'}

        # test 'all'
        expectedResults = {}
        for analysisName in cores:
            expectedResults[analysisName] = True
        doTest("['all']", expectedResults)

        # test 'all_<category>' and ['all', 'no_<category>']
        for category in set(categories.values()):
            expectedResults = {}
            for analysisName in categories:
                expectedResults[analysisName] = \
                    (categories[analysisName] == category)
            doTest("['all_{}']".format(category), expectedResults)

            expectedResults = {}
            for analysisName in categories:
                expectedResults[analysisName] = \
                    (categories[analysisName] != category)
            doTest("['all', 'no_{}']".format(category), expectedResults)

        # test 'all_<core>' and ['all', 'no_<core>']
        for core in set(cores.values()):
            expectedResults = {}
            for analysisName in cores:
                expectedResults[analysisName] = \
                    (cores[analysisName] == core)
            doTest("['all_{}']".format(core), expectedResults)

            expectedResults = {}
            for analysisName in cores:
                expectedResults[analysisName] = \
                    (cores[analysisName] != core)
            doTest("['all','no_{}']".format(core), expectedResults)

        # test each analysis individually
        for analysisName in cores:
            expectedResults = {}
            for otherAnalysis in cores:
                expectedResults[otherAnalysis] = \
                    (analysisName == otherAnalysis)
            doTest("['{}']".format(analysisName), expectedResults)

        # test a non-existent analysis
        expectedResults = {}
        for analysisName in cores:
            expectedResults[analysisName] = False
        doTest("['fakeAnalysis']", expectedResults)

        # test ['all', 'no_ocean', 'all_timeSeries']
        expectedResults = {}
        for analysisName in cores:
            expectedResults[analysisName] = True
        for analysisName in cores:
            if cores[analysisName] == 'ocean':
                expectedResults[analysisName] = False
        for analysisName in categories:
            if categories[analysisName] == 'timeSeries':
                expectedResults[analysisName] = True
        doTest("['all', 'no_ocean', 'all_timeSeries']", expectedResults)

        # test ['all', 'no_timeSeriesOHC']
        expectedResults = {}
        for analysisName in cores:
            expectedResults[analysisName] = True
        expectedResults['timeSeriesOHC'] = False
        doTest("['all', 'no_timeSeriesOHC']", expectedResults)


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
