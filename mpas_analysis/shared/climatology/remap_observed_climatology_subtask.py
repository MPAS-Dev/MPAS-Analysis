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
import os
import os.path
import xarray as xr

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories
from mpas_analysis.shared.io import write_netcdf_with_fill

from mpas_analysis.shared.climatology.climatology import get_remapper, \
    remap_and_write_climatology, compute_climatology

from mpas_analysis.shared.climatology.comparison_descriptors import \
    get_comparison_descriptor


class RemapObservedClimatologySubtask(AnalysisTask):
    """
    An analysis task for comparison of 2D model fields against observations.

    Attributes
    ----------
    seasons : list of str
       A list of seasons (keys in ``constants.monthDictionary``) over
       which the climatology should be computed.

    fileName : str
        The name of the observation file

    outFilePrefix : str
        The prefix in front of output files and mapping files, typically the
        name of the field being remapped

    comparisonGridNames : list of str
        The name(s) of the comparison grid to use for remapping.
    """

    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, seasons, fileName, outFilePrefix,
                 comparisonGridNames=['latlon'],
                 subtaskName='remapObservations'):

        """
        Construct one analysis subtask for each plot (i.e. each season and
        comparison grid) and a subtask for computing climatologies.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent (main) task for this subtask

        seasons : list of str
           A list of seasons (keys in ``constants.monthDictionary``) over
           which the climatology should be computed.

        fileName : str
            The name of the observation file

        outFilePrefix : str
            The prefix in front of output files and mapping files, typically
            the name of the field being remapped

        comparisonGridNames : list of str
            optional
            The name(s) of the comparison grid to use for remapping.

        subtaskName : str, optional
            The name of the subtask
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.seasons = seasons
        self.fileName = fileName
        self.outFilePrefix = outFilePrefix
        self.comparisonGridNames = comparisonGridNames

        config = parentTask.config
        taskName = parentTask.taskName
        tags = parentTask.tags
        componentName = parentTask.componentName

        # call the constructor from the base class (AnalysisTask)
        super(RemapObservedClimatologySubtask, self).__init__(
            config=config, taskName=taskName, subtaskName=subtaskName,
            componentName=componentName, tags=tags)

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(RemapObservedClimatologySubtask, self).setup_and_check()

        # we set up the remappers here because ESFM_RegridWeightGen seems to
        # have trouble if it runs in another process (or in several at once)
        self._setup_remappers(self.fileName)

        # build the observational data set and write it out to a file, to
        # be read back in during the run_task() phase
        obsFileName = self.get_file_name(stage='original')
        if not os.path.exists(obsFileName):
            ds = self.build_observational_dataset(self.fileName)
            write_netcdf_with_fill(ds, obsFileName)

    def run_task(self):
        """
        Performs remapping of obsrevations to the comparsion grid
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config

        obsFileName = self.get_file_name(stage='original')
        if not os.path.isfile(obsFileName):
            raise OSError('Obs file {} not found.'.format(
                obsFileName))

        for comparisonGridName in self.comparisonGridNames:
            for season in self.seasons:

                remappedFileName = self.get_file_name(
                    stage='remapped',
                    season=season,
                    comparisonGridName=comparisonGridName)

                if not os.path.exists(remappedFileName):

                    ds = xr.open_dataset(obsFileName)

                    climatologyFileName = self.get_file_name(
                        stage='climatology',
                        season=season,
                        comparisonGridName=comparisonGridName)
                    if 'month' in ds.coords and 'year' in ds.coords:
                        # this data set is not yet a climatology, so compute
                        # the climatology
                        monthValues = constants.monthDictionary[season]
                        seasonalClimatology = compute_climatology(
                            ds, monthValues, maskVaries=True)
                    else:
                        # We don't have month or year arrays to compute a
                        # climatology so assume this already is one
                        seasonalClimatology = ds

                    write_netcdf_with_fill(seasonalClimatology,
                                           climatologyFileName)

                    remapper = self.remappers[comparisonGridName]

                    if remapper.map_filename is None:
                        # no need to remap because the observations are on the
                        # comparison grid already
                        os.symlink(climatologyFileName, remappedFileName)
                    else:
                        remap_and_write_climatology(
                            config, seasonalClimatology,
                            climatologyFileName,
                            remappedFileName, remapper,
                            logger=self.logger)

    def get_observation_descriptor(self, fileName):
        """
        get a MeshDescriptor for the observation grid.  A subclass derived from
        this class must override this method to create the appropriate
        descriptor

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

        return None

    def build_observational_dataset(self, fileName):
        """
        read in the data sets for observations, and possibly rename some
        variables and dimensions.  A subclass derived from this class must
        override this method to create the appropriate data set

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

        return None

    def get_file_name(self, stage, season=None, comparisonGridName=None):

        """
        Given config options, the name of a field and a string identifying the
        months in a seasonal climatology, returns the full path for MPAS
        climatology files before and after remapping.

        Parameters
        ----------
        stage : {'original', 'climatology', 'remapped'}
            The stage of the masking and remapping process

        season : str, optional
            One of the seasons in ``constants.monthDictionary``

        comparisonGridName : str, optional
            The name of the comparison grid to use for remapping.

        Returns
        -------
        fileName : str
            The path to the climatology file for the specified season.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config
        obsSection = '{}Observations'.format(self.componentName)
        if comparisonGridName is None:
            # just needed for getting the obs. grid name, so doesn't matter
            # which comparison grid
            remapper = self.remappers[self.comparisonGridNames[0]]
        else:
            remapper = self.remappers[comparisonGridName]

        obsGridName = remapper.src_descriptor.mesh_name

        outFilePrefix = self.outFilePrefix

        if stage in ['original', 'climatology']:
            climatologyDirectory = build_config_full_path(
                config=config, section='output',
                relativePathOption='climatologySubdirectory',
                relativePathSection=obsSection)

            make_directories(climatologyDirectory)

            if stage == 'original':
                fileName = '{}/{}_{}.nc'.format(
                    climatologyDirectory, outFilePrefix, obsGridName)
            else:
                fileName = '{}/{}_{}_{}.nc'.format(
                    climatologyDirectory, outFilePrefix, obsGridName, season)

        elif stage == 'remapped':
            remappedDirectory = build_config_full_path(
                config=config, section='output',
                relativePathOption='remappedClimSubdirectory',
                relativePathSection=obsSection)

            make_directories(remappedDirectory)

            comparisonGridName = remapper.dst_descriptor.mesh_name
            fileName = '{}/{}_{}_to_{}_{}.nc'.format(
                remappedDirectory, outFilePrefix, obsGridName,
                comparisonGridName, season)

        else:
            raise ValueError('Unknown stage {}'.format(stage))

        return fileName

    def _setup_remappers(self, fileName):
        """
        Set up the remappers for remapping from observations to the comparison
        grids.

        Parameters
        ----------
        fileName : str
            The name of the observation file used to determine the source grid
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config

        sectionName = '{}Observations'.format(self.componentName)

        obsDescriptor = self.get_observation_descriptor(fileName)

        outFilePrefix = self.outFilePrefix
        self.remappers = {}
        for comparisonGridName in self.comparisonGridNames:
            comparisonDescriptor = get_comparison_descriptor(
                config, comparison_grid_name=comparisonGridName)

            self.remappers[comparisonGridName] = get_remapper(
                config=config,
                sourceDescriptor=obsDescriptor,
                comparisonDescriptor=comparisonDescriptor,
                mappingFilePrefix='map_obs_{}'.format(outFilePrefix),
                method=config.get(sectionName,
                                  'interpolationMethod'),
                logger=self.logger)
