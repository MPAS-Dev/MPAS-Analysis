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
Unit test infrastructure for MpasAnalysisConfigParser

Xylar Asay-Davis, Phillip J. Wolfram
01/31/2017
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import pytest
import six
from six.moves import configparser
from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.configuration import MpasAnalysisConfigParser
from mpas_analysis.test import requires_numpy


@pytest.mark.usefixtures("loaddatadir")
class TestMPASAnalysisConfigParser(TestCase):
    def setup_config(self):
        configPath = self.datadir.join('analysis.cfg')
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
        self.assertEqual(self.config.getExpression('Test', 'testInt'), 15)

        self.assertEqual(self.config.getfloat('Test', 'testFloat'), 18.0)
        self.assertEqual(self.config.getExpression('Test', 'testFloat'), 18.0)

        self.assertEqual(self.config.getfloat('Test', 'testFloat2'), 3.)
        self.assertEqual(self.config.getExpression('Test', 'testFloat2'), 3.)

        self.assertEqual(self.config.getboolean('Test', 'testBool'), True)
        self.assertEqual(self.config.getExpression('Test', 'testBool'), True)

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

        with six.assertRaisesRegex(
                self,
                configparser.NoOptionError,
                "No option 'doesntexist' in section: 'Test'"):
            self.config.getExpression(str('Test'), str('doesntexist'))

    @requires_numpy
    def test_read_config_numpy(self):
        self.setup_config()

        # tests numpy evaluation capability
        import numpy as np
        for testname in ['testNumpyarange' + str(ii) for ii in np.arange(3)]:
            self.assertArrayEqual(self.config.getExpression('TestNumpy',
                                                            testname,
                                                            usenumpyfunc=True),
                                  np.arange(0, 1, 10))
        for testname in ['testNumpylinspace' + str(ii) for ii in np.arange(3)]:
            self.assertArrayEqual(self.config.getExpression('TestNumpy',
                                                            testname,
                                                            usenumpyfunc=True),
                                  np.linspace(0, 1, 10))
        for testNumpy in ['testNumpypi' + str(ii) for ii in np.arange(3)] + \
                ['testNumpyPi']:
            self.assertEqual(self.config.getExpression('TestNumpy', testNumpy,
                                                       usenumpyfunc=True),
                             np.pi)
        with six.assertRaisesRegex(
                self,
                AssertionError,
                "'__' is not allowed in .* for `usenumpyfunc=True`"):
            self.config.getExpression('TestNumpy', 'testBadStr',
                                      usenumpyfunc=True),

    def test_get_with_default(self):
        self.setup_config()

        def check_get_with_default(name, value, dtype):
            # test an options that doesn't exist using getWithDefault
            var = self.config.getWithDefault('sst_modelvsobs', name, value)
            assert isinstance(var, dtype)
            self.assertEqual(var, value)

        # test several types with getWithDefault
        check_get_with_default(name='aBool', value=True, dtype=bool)
        check_get_with_default(name='anInt', value=1, dtype=six.integer_types)
        check_get_with_default(name='aFloat', value=1.0, dtype=float)
        check_get_with_default(name='aList', value=[1, 2, 3], dtype=list)
        check_get_with_default(name='aTuple', value=(1, 2, 3), dtype=tuple)
        check_get_with_default(name='aDict', value={'blah': 1}, dtype=dict)
        check_get_with_default(name='aStr', value='blah',
                               dtype=six.string_types)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
