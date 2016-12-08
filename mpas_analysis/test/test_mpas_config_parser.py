"""
Unit test infrastructure for MpasAnalysisConfigParser

Xylar Asay-Davis
12/05/2016
"""

import pytest
from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.configuration.MpasAnalysisConfigParser \
    import MpasAnalysisConfigParser


@pytest.mark.usefixtures("loaddatadir")
class TestMPASAnalysisConfigParser(TestCase):
    def setup_config(self):
        configPath = self.datadir.join('config.analysis')
        self.config = MpasAnalysisConfigParser()
        self.config.read(str(configPath))

    def check_container(self, container, container_type, item_type):
        assert isinstance(container, container_type)
        for item in container:
            assert isinstance(item, item_type)

    def test_read_config(self):
        self.setup_config()

        colorMapName = self.config.get('sst_modelvsobs', 'cmapDiff')
        self.assertEqual(colorMapName, 'coolwarm')

        self.assertEqual(self.config.getint('Test', 'testInt'), 15)
        self.assertEqual(self.config.getfloat('Test', 'testFloat'), 18.)
        self.assertEqual(self.config.getboolean('Test', 'testBool'), True)

        testList = self.config.getExpression('sst_modelvsobs',
                                             'cmapIndicesModelObs')
        self.check_container(testList, list, int)
        self.assertEqual(testList, [0, 40, 80, 110, 140, 170, 200, 230, 255])

        testList = self.config.getExpression('sst_modelvsobs',
                                             'cmapIndicesModelObs',
                                             elementType=float)
        self.check_container(testList, list, float)
        self.assertEqual(testList, [0., 40., 80., 110., 140., 170., 200.,
                                    230., 255.])

        testList = self.config.getExpression('sst_modelvsobs',
                                             'comparisonTimes')
        self.check_container(testList, list, str)
        self.assertEqual(testList, ['JFM', 'JAS', 'ANN'])

        testList = self.config.getExpression('Test', 'testList')
        self.check_container(testList, list, float)
        self.assertEqual(testList, [0.5, 0.1, 0.5])

        testTuple = self.config.getExpression('Test', 'testTuple')
        assert isinstance(testTuple, tuple)
        self.assertEqual(testTuple, (5, 0.1, 'string'))

        testDict = self.config.getExpression('Test', 'testDict')
        assert isinstance(testDict, dict)
        self.assertEqual(testDict, {'key1': 'string',
                                    'key2': -12,
                                    'key3': False})


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
