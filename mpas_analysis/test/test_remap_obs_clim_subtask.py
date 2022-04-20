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
Unit test infrastructure for MpasClimatologyTask.

Xylar Asay-Davis
"""

import pytest
import tempfile
import shutil
import os
import xarray
from pyremap import LatLonGridDescriptor

from mpas_tools.config import MpasConfigParser

from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.shared.climatology import RemapObservedClimatologySubtask
from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories


class RemapObservedMLDClimatology(RemapObservedClimatologySubtask):
    """
    A subtask for reading and remapping MLD observations
    """

    # Authors
    # -------
    # Xylar Asay-Davis

    def get_observation_descriptor(self, fileName):
        """
        get a MeshDescriptor for the observation grid

        Parameters
        ----------
        fileName : str
            observation file name describing the source grid

        Returns
        -------
        obsDescriptor : ``MeshDescriptor``
            The descriptor for the observation grid
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # create a descriptor of the observation grid using the lat/lon
        # coordinates
        obsDescriptor = LatLonGridDescriptor.read(fileName=fileName,
                                                  latVarName='lat',
                                                  lonVarName='lon')
        return obsDescriptor

    def build_observational_dataset(self, fileName):
        """
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Parameters
        ----------
        fileName : str
            observation file name

        Returns
        -------
        dsObs : ``xarray.Dataset``
            The observational dataset
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # Load MLD observational data
        dsObs = xarray.open_dataset(fileName)

        dsObs.coords['month'] = dsObs['month']
        dsObs.coords['year'] = dsObs['year']

        return dsObs


@pytest.mark.usefixtures("loaddatadir")
class TestRemapObsClimSubtask(TestCase):
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def setup_config(self):
        configPath = self.datadir.join('remap_obs.cfg')
        config = MpasConfigParser()
        config.add_from_file(str(configPath))
        config.set('input', 'baseDirectory', str(self.datadir))
        config.set('diagnostics', 'base_path', str(self.datadir))
        config.set('oceanObservations', 'obsSubdirectory', '.')
        config.set('output', 'baseDirectory', str(self.test_dir))
        return config

    def setup_subtask(self):
        config = self.setup_config()
        parent = AnalysisTask(config, 'taskName', 'ocean')
        obsFileName = '{}/mld_4.0x4.0degree.nc'.format(str(self.datadir))
        remapObsTask = RemapObservedMLDClimatology(
            parentTask=parent,
            seasons=['ANN'], fileName=obsFileName,
            outFilePrefix='mld',
            comparisonGridNames=['latlon', 'antarctic'])

        remapObsTask.setup_and_check()
        return remapObsTask

    def test_subtask_run_analysis(self):
        remapSubtask = self.setup_subtask()

        config = remapSubtask.config
        logsDirectory = build_config_full_path(config, 'output',
                                               'logsSubdirectory')
        make_directories(logsDirectory)
        make_directories('{}/configs/'.format(logsDirectory))

        remapSubtask.run(writeLogFile=False)

        for comparisonGridName in remapSubtask.comparisonGridNames:
            for season in remapSubtask.seasons:
                for stage in ['original', 'climatology', 'remapped']:
                    fileName = remapSubtask.get_file_name(
                        season=season, stage=stage,
                        comparisonGridName=comparisonGridName)
                    assert (os.path.exists(fileName))

    def test_subtask_get_file_name(self):
        remapSubtask = self.setup_subtask()

        for comparisonGridName in ['latlon', 'antarctic', None]:
            fileName = remapSubtask.get_file_name(
                stage='original', comparisonGridName=comparisonGridName)
            assert (fileName == '{}/clim/obs/mld_4.0x4.0degree.nc'.format(
                str(self.test_dir)))

        stage = 'climatology'
        for comparisonGridName in ['latlon', 'antarctic', None]:
            fileName = remapSubtask.get_file_name(stage=stage,
                                                  season='JFM',
                                                  comparisonGridName='latlon')
            assert (fileName == '{}/clim/obs/mld_4.0x4.0degree_JFM.nc'.format(
                str(self.test_dir)))

        stage = 'remapped'
        fileName = remapSubtask.get_file_name(stage=stage, season='JFM',
                                              comparisonGridName='latlon')
        assert (fileName == '{}/clim/obs/remapped/mld_4.0x4.0degree_to_'
                            '0.5x0.5degree_JFM.nc'.format(str(self.test_dir)))

        fileName = remapSubtask.get_file_name(stage=stage, season='JFM',
                                              comparisonGridName='antarctic')
        assert (fileName == '{}/clim/obs/remapped/mld_4.0x4.0degree_to_'
                            '6000.0x6000.0km_10.0km_Antarctic_stereo_JFM.nc'.format(
            str(self.test_dir)))
