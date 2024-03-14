# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/main/LICENSE
"""
Unit test infrastructure for MpasClimatologyTask.

Xylar Asay-Davis
"""

import pytest
import tempfile
import shutil
import os

from mpas_tools.config import MpasConfigParser

from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.shared.climatology import MpasClimatologyTask, \
    RemapMpasClimatologySubtask
from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.analysis_task import \
    update_time_bounds_from_file_names
from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories
from mpas_analysis.shared.constants import constants


@pytest.mark.usefixtures("loaddatadir")
class TestMpasClimatologyTask(TestCase):
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def setup_config(self):
        configPath = self.datadir.join('QU240.cfg')
        config = MpasConfigParser()
        config.add_from_file(str(configPath))
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
        seasons = list(mpasClimatologyTask.variableList.keys())

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

    def verify_variables_for_season(self, mpasClimatologyTask, variableList,
                                    season):
        assert(season in mpasClimatologyTask.variableList.keys())
        monthValues = constants.monthDictionary[season]
        monthNames = [constants.abrevMonthNames[month-1] for month in
                      monthValues]
        assert(variableList == mpasClimatologyTask.variableList[season])
        for monthName in monthNames:
            assert(monthName in mpasClimatologyTask.variableList.keys())
            assert(variableList ==
                   mpasClimatologyTask.variableList[monthName])

    def test_add_variables(self):
        mpasClimatologyTask = self.setup_task()
        variableList, seasons = self.add_variables(mpasClimatologyTask)

        for season in seasons:
            self.verify_variables_for_season(mpasClimatologyTask, variableList,
                                             season)

        # add a variable and season already in the list
        mpasClimatologyTask.add_variables(variableList=[variableList[0]],
                                          seasons=[seasons[-1]])

        for season in seasons:
            self.verify_variables_for_season(mpasClimatologyTask, variableList,
                                             season)

    def test_get_file_name(self):
        mpasClimatologyTask = self.setup_task()
        variableList, seasons = self.add_variables(mpasClimatologyTask)

        fileName = mpasClimatologyTask.get_file_name(season='JFM')
        assert(fileName == '{}/clim/mpas/avg/unmasked_oQU240/'
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

        for season in mpasClimatologyTask.variableList:
            fileName = mpasClimatologyTask.get_file_name(season=season)
            assert(os.path.exists(fileName))

    def test_update_climatology_bounds_and_create_symlinks(self):
        mpasClimatologyTask = self.setup_task()
        config = mpasClimatologyTask.config

        # make sure the start and end years stay unchanged when we use
        # the start and end years already in the config file
        startYear = 2
        endYear = 2
        startDate = '{:04d}-01-01_00:00:00'.format(startYear)
        endDate = '{:04d}-12-31_23:59:59'.format(endYear)

        update_time_bounds_from_file_names(config, 'climatology', 'ocean')
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

        with self.assertRaisesRegex(ValueError,
                                    'climatology start and/or end year '
                                    'different from requested'):
            update_time_bounds_from_file_names(config, 'climatology', 'ocean')

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
        assert(fileName == '{}/clim/mpas/avg/masked/ssh_oQU240/'
               'mpaso_JFM_000201_000203_climo.nc'.format(str(self.test_dir)))

        fileName = remapSubtask.get_remapped_file_name(
            season='JFM', comparisonGridName='latlon')
        assert(fileName == '{}/clim/mpas/avg/remapped/'
               'ssh_oQU240_to_0.5x0.5degree/mpaso_JFM_000201_000203_climo.nc'
               ''.format(str(self.test_dir)))
