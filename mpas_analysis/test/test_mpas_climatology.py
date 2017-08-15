"""
Unit test infrastructure for climatologies.

Xylar Asay-Davis
04/11/2017
"""

import pytest
import tempfile
import shutil
import os
from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.configuration.MpasAnalysisConfigParser \
    import MpasAnalysisConfigParser
from mpas_analysis.shared.climatology import MpasClimatology, \
    get_comparison_descriptor
from mpas_analysis.shared.constants import constants


@pytest.mark.usefixtures("loaddatadir")
class TestMpasClimatology(TestCase):

    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def setup_config(self):
        config = MpasAnalysisConfigParser()
        config.read('{}/config.GQU240'.format(self.datadir))
        config.set('input', 'baseDirectory', str(self.datadir))
        config.set('output', 'baseDirectory', str(self.test_dir))
        return config

    def setup_climatology(self, config, seasons, comparisonGridNames=None):
        variableList = ['timeMonthly_avg_ssh', 'timeMonthly_avg_tThreshMLD']
        climatologyTask = \
            MpasClimatology(config=config,
                            variableList=variableList,
                            taskSuffix='SSH_MLD',
                            componentName='ocean',
                            comparisonGridNames=comparisonGridNames,
                            seasons=seasons,
                            tags=['climatology', 'ssh', 'mld'])

        climatologyTask.setup_and_check()

        climatologyTask.run()

        return climatologyTask

    def test_seasons_none(self):
        config = self.setup_config()
        climatologyTask = \
            self.setup_climatology(config=config,
                                   seasons=['none'])

        for month in constants.abrevMonthNames:
            fileName = climatologyTask.get_ncclimo_file_name(
                    season=month, stage='unmasked')
            assert os.path.exists(fileName)

    def test_seasons(self):
        config = self.setup_config()
        seasons = ['JFM', 'JAS', 'FM', 'ON', 'ANN']
        climatologyTask = \
            self.setup_climatology(config=config,
                                   seasons=seasons)

        for season in seasons:
            fileName = climatologyTask.get_ncclimo_file_name(
                    season=season, stage='unmasked')
            assert os.path.exists(fileName)

    def test_remap(self):
        config = self.setup_config()
        seasons = ['JFM', 'JAS', 'FM', 'ON', 'ANN']
        comparisonGridNames = ['latlon', 'antarctic']
        climatologyTask = \
            self.setup_climatology(config=config,
                                   seasons=seasons,
                                   comparisonGridNames=comparisonGridNames)

        for season in seasons:
            fileName = climatologyTask.get_ncclimo_file_name(
                    season=season, stage='masked')
            assert os.path.exists(fileName)

            for comparisonGridName in comparisonGridNames:
                comparisonDescriptor = \
                    get_comparison_descriptor(config, comparisonGridName)
                fileName = climatologyTask.get_ncclimo_file_name(
                        season=season, stage='remapped',
                        comparisonDescriptor=comparisonDescriptor)
            assert os.path.exists(fileName)


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
