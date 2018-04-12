# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#
"""
Unit test infrastructure for MpasClimatologyTask.

Xylar Asay-Davis
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import pytest
import tempfile
import shutil
import os

from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.configuration import MpasAnalysisConfigParser
from mpas_analysis.shared.climatology import MpasClimatologyTask, \
    RemapMpasClimatologySubtask
from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories


@pytest.mark.usefixtures("loaddatadir")
class TestMpasClimatologyTask(TestCase):
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def setup_config(self):
        configPath = self.datadir.join('config.QU240')
        config = MpasAnalysisConfigParser()
        config.read(str(configPath))
        config.set('input', 'baseDirectory', str(self.datadir))
        config.set('output', 'baseDirectory', str(self.test_dir))
        return config

    def setup_task(self):
        config = self.setup_config()
        mpasClimatologyTask = MpasClimatologyTask(config=config,
                                                  componentName='ocean')

        mpasClimatologyTask.setup_and_check()
        return mpasClimatologyTask

    def setup_subtask(self, mpasClimatologyTask):
        parentTask = AnalysisTask(
                config=mpasClimatologyTask.config, taskName='fake',
                componentName=mpasClimatologyTask.componentName,
                tags=['climatology'])
        climatologyName = 'ssh'
        variableList = ['timeMonthly_avg_ssh']
        seasons = [mpasClimatologyTask.seasons[0]]

        remapSubtask = RemapMpasClimatologySubtask(
                mpasClimatologyTask, parentTask, climatologyName,
                variableList, seasons, comparisonGridNames=['latlon'])

        remapSubtask.setup_and_check()
        return remapSubtask

    def add_variables(self, mpasClimatologyTask):
        variableList = ['timeMonthly_avg_ssh', 'timeMonthly_avg_tThreshMLD']
        seasons = ['JFM', 'JJA', 'ANN']
        mpasClimatologyTask.add_variables(variableList=variableList,
                                          seasons=seasons)

        return variableList, seasons

    def test_add_variables(self):
        mpasClimatologyTask = self.setup_task()
        variableList, seasons = self.add_variables(mpasClimatologyTask)

        assert(variableList == mpasClimatologyTask.variableList)
        assert(seasons == mpasClimatologyTask.seasons)

        # add a variable and season already in the list
        mpasClimatologyTask.add_variables(variableList=[variableList[0]],
                                          seasons=[seasons[-1]])

        # make sure the lists still match (extra redundant varible and season
        # weren't added)
        assert(variableList == mpasClimatologyTask.variableList)
        assert(seasons == mpasClimatologyTask.seasons)

    def test_get_file_name(self):
        mpasClimatologyTask = self.setup_task()
        variableList, seasons = self.add_variables(mpasClimatologyTask)

        fileName = mpasClimatologyTask.get_file_name(season='JFM')
        assert(fileName == '{}/clim/mpas/unmasked_oQU240/'
               'mpaso_JFM_000201_000203_climo.nc'.format(str(self.test_dir)))

    def test_run_analysis(self):
        mpasClimatologyTask = self.setup_task()
        self.add_variables(mpasClimatologyTask)

        config = mpasClimatologyTask.config
        logsDirectory = build_config_full_path(config, 'output',
                                               'logsSubdirectory')
        make_directories(logsDirectory)
        make_directories('{}/configs/'.format(logsDirectory))

        mpasClimatologyTask.run(writeLogFile=False)

        for season in mpasClimatologyTask.seasons:
            fileName = mpasClimatologyTask.get_file_name(season=season)
            assert(os.path.exists(fileName))

    def test_update_climatology_bounds_and_create_symlinks(self):
        mpasClimatologyTask = self.setup_task()
        config = mpasClimatologyTask.config

        # first make sure the start and end years stay unchanged when we use
        # the start and end years already in the config file
        startYear = 2
        endYear = 2
        startDate = '{:04d}-01-01_00:00:00'.format(startYear)
        endDate = '{:04d}-12-31_23:59:59'.format(endYear)

        mpasClimatologyTask._update_climatology_bounds()
        mpasClimatologyTask._create_symlinks()

        assert(mpasClimatologyTask.startYear == startYear)
        assert(mpasClimatologyTask.endYear == endYear)
        assert(mpasClimatologyTask.startDate == startDate)
        assert(mpasClimatologyTask.endDate == endDate)

        # Now, set the the start and end years out of range and make sure they
        # get changed back to the values that are in range
        startYear = 1
        endYear = 5
        startDate = '{:04d}-01-01_00:00:00'.format(startYear)
        endDate = '{:04d}-12-31_23:59:59'.format(endYear)

        config.set('climatology', 'startYear', str(startYear))
        config.set('climatology', 'endYear', str(endYear))
        config.set('climatology', 'startDate', startDate)
        config.set('climatology', 'endDate', endDate)

        mpasClimatologyTask._update_climatology_bounds()
        mpasClimatologyTask._create_symlinks()

        startYear = 2
        endYear = 2
        startDate = '{:04d}-01-01_00:00:00'.format(startYear)
        endDate = '{:04d}-12-31_23:59:59'.format(endYear)

        assert(mpasClimatologyTask.startYear == startYear)
        assert(mpasClimatologyTask.endYear == endYear)
        assert(mpasClimatologyTask.startDate == startDate)
        assert(mpasClimatologyTask.endDate == endDate)

    def test_subtask_run_analysis(self):
        mpasClimatologyTask = self.setup_task()
        self.add_variables(mpasClimatologyTask)
        remapSubtask = self.setup_subtask(mpasClimatologyTask)

        config = mpasClimatologyTask.config
        logsDirectory = build_config_full_path(config, 'output',
                                               'logsSubdirectory')
        make_directories(logsDirectory)
        make_directories('{}/configs/'.format(logsDirectory))

        mpasClimatologyTask.run(writeLogFile=False)
        remapSubtask.run(writeLogFile=False)

        for season in remapSubtask.seasons:
            fileName = remapSubtask.get_masked_file_name(season=season)
            assert(os.path.exists(fileName))

            fileName = remapSubtask.get_remapped_file_name(
                    season=season, comparisonGridName='latlon')
            assert(os.path.exists(fileName))

    def test_subtask_get_file_name(self):
        mpasClimatologyTask = self.setup_task()
        variableList, seasons = self.add_variables(mpasClimatologyTask)
        remapSubtask = self.setup_subtask(mpasClimatologyTask)

        fileName = remapSubtask.get_masked_file_name(season='JFM')
        assert(fileName == '{}/clim/mpas/masked/ssh_oQU240/'
               'mpaso_JFM_000201_000203_climo.nc'.format(str(self.test_dir)))

        fileName = remapSubtask.get_remapped_file_name(
                season='JFM', comparisonGridName='latlon')
        assert(fileName == '{}/clim/mpas/remapped/ssh_oQU240_to_0.5x0.5degree/'
               'mpaso_JFM_000201_000203_climo.nc'.format(str(self.test_dir)))
