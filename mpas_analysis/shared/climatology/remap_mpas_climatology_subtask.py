
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr
import numpy
import os

from ..analysis_task import AnalysisTask

from ..constants import constants

from ..io.utility import build_config_full_path, make_directories
from ..io import write_netcdf

from .climatology import get_remapper
from .comparison_descriptors import get_comparison_descriptor

from ..grid import MpasMeshDescriptor

from ..mpas_xarray import mpas_xarray


class RemapMpasClimatologySubtask(AnalysisTask):  # {{{
    '''
    An analysis tasks for computing climatologies from output from the
    ``timeSeriesStatsMonthly`` analysis member.

    Attributes
    ----------

    climatologyName :  str
        A name that describes the climatology (e.g. a short version of
        the important field(s) in the climatology) used to name the
        subdirectories for each stage of the climatology

    variableList : list of str
        A list of variable names in ``timeSeriesStatsMonthly`` to be
        included in the climatologies

    iselValues : dict
        A dictionary of dimensions and indices (or ``None``) used to extract
        a slice of the MPAS field.

    seasons : list of str
        A list of seasons (keys in ``shared.constants.monthDictionary``)
        over which the climatology should be computed or ['none'] if only
        monthly climatologies are needed.

    comparisonGridNames : list of {'latlon', 'antarctic'}
        The name(s) of the comparison grid to use for remapping.

    restartFileName : str
        If ``comparisonGridName`` is not ``None``, the name of a restart
        file from which the MPAS mesh can be read.

    Authors
    -------
    Xylar Asay-Davis
    '''

    def __init__(self, mpasClimatologyTask, parentTask, climatologyName,
                 variableList, seasons, comparisonGridNames=['latlon'],
                 iselValues=None):
        # {{{
        '''
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped

        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        climatologyName : str
            A name that describes the climatology (e.g. a short version of
            the important field(s) in the climatology) used to name the
            subdirectories for each stage of the climatology

        variableList : list of str
            A list of variable names in ``timeSeriesStatsMonthly`` to be
            included in the climatologies

        seasons : list of str, optional
            A list of seasons (keys in ``shared.constants.monthDictionary``)
            to be computed or ['none'] (not ``None``) if only monthly
            climatologies are needed.

        comparisonGridNames : list of {'latlon', 'antarctic'}, optional
            The name(s) of the comparison grid to use for remapping.

        iselValues : dict, optional
            A dictionary of dimensions and indices (or ``None``) used to
            extract a slice of the MPAS field(s).

        Authors
        -------
        Xylar Asay-Davis
        '''
        tags = ['climatology']

        # call the constructor from the base class (AnalysisTask)
        super(RemapMpasClimatologySubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            subtaskName='remapMpasClimatology',
            componentName=parentTask.componentName,
            tags=tags)

        self.variableList = variableList
        self.seasons = seasons
        self.comparisonGridNames = comparisonGridNames
        self.iselValues = iselValues
        self.climatologyName = climatologyName
        self.mpasClimatologyTask = mpasClimatologyTask

        self.run_after(mpasClimatologyTask)

        parentTask.add_subtask(self)

        # this is a stopgap until MPAS implements the _FillValue attribute
        # correctly
        self._fillValue = -9.99999979021476795361e+33

        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        IOError :
            If a restart file is not available from which to read mesh
            information or if no history files are available from which to
            compute the climatology in the desired time range.

        Authors
        -------
        Xylar Asay-Davis
        '''

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(RemapMpasClimatologySubtask, self).setup_and_check()

        try:
            self.restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS restart file found: need at least one '
                          'restart file to perform remapping of '
                          'climatologies.')

        # we set up the remapper here because ESFM_RegridWeightGen seems to
        # have trouble if it runs in another process (or in several at once)
        self._setup_remappers()

        # don't add the variables and seasons to mpasClimatologyTask until
        # we're sure this subtask is supposed to run
        self.mpasClimatologyTask.add_variables(self.variableList, self.seasons)

        self._setup_file_names()

        # make the mapping directory, because doing so within each process
        # seems to be giving ESMF_RegridWeightGen some trouble
        mappingSubdirectory = build_config_full_path(self.config, 'output',
                                                     'mappingSubdirectory')
        make_directories(mappingSubdirectory)

        # }}}

    def run_task(self):  # {{{
        '''
        Compute the requested climatologies

        Authors
        -------
        Xylar Asay-Davis
        '''

        self.logger.info('\nRemapping climatology {}'.format(
            self.climatologyName))

        dsMask = xr.open_dataset(self.mpasClimatologyTask.inputFiles[0])
        dsMask = mpas_xarray.subset_variables(dsMask, self.variableList)
        iselValues = {'Time': 0}
        if self.iselValues is not None:
            iselValues.update(self.iselValues)
        # select only Time=0 and possibly only the desired vertical
        # slice
        dsMask = dsMask.isel(**iselValues)

        for season in self.seasons:
            self._mask_climatologies(season, dsMask)

        for comparisonGridName in self.comparisonGridNames:

            for season in self.seasons:

                maskedClimatologyFileName = self.get_file_name(
                        season, 'masked', comparisonGridName)

                remappedFileName = self.get_file_name(
                        season, 'remapped', comparisonGridName)

                if not os.path.exists(remappedFileName):
                    self._remap(inFileName=maskedClimatologyFileName,
                                outFileName=remappedFileName,
                                remapper=self.remappers[comparisonGridName],
                                comparisonGridName=comparisonGridName)
        # }}}

    def get_file_name(self, season, stage, comparisonGridName=None):  # {{{
        """
        Given config options, the name of a field and a string identifying the
        months in a seasonal climatology, returns the full path for MPAS
        climatology files before and after remapping.

        Parameters
        ----------
        season : str
            One of the seasons in ``constants.monthDictionary``

        stage : {'masked', 'remapped'}
            The stage of the masking and remapping process

        comparisonGridName : {'latlon', 'antarctic'}, optional
            The name of the comparison grid to use for remapping.

        Returns
        -------
        fileName : str
            The path to the climatology file for the specified season.

        Authors
        -------
        Xylar Asay-Davis
        """

        if stage == 'remapped':
            key = (season, stage, comparisonGridName)
        else:
            key = (season, stage)
        fileName = self._outputFiles[key]

        return fileName  # }}}

    def _setup_remappers(self):  # {{{
        """
        Set up the remappers for remapping from the MPAS to the comparison
        grids.

        Authors
        -------
        Xylar Asay-Davis
        """
        config = self.config

        # make reamppers
        mappingFilePrefix = 'map'
        self.remappers = {}
        for comparisonGridName in self.comparisonGridNames:

            comparisonDescriptor = get_comparison_descriptor(
                    config, comparisonGridName)
            self.comparisonGridName = comparisonDescriptor.meshName
            mpasDescriptor = MpasMeshDescriptor(
                self.restartFileName, meshName=config.get('input',
                                                          'mpasMeshName'))
            self.mpasMeshName = mpasDescriptor.meshName

            self.remappers[comparisonGridName] = get_remapper(
                config=config, sourceDescriptor=mpasDescriptor,
                comparisonDescriptor=comparisonDescriptor,
                mappingFilePrefix=mappingFilePrefix,
                method=config.get('climatology', 'mpasInterpolationMethod'),
                logger=self.logger)

        # }}}

    def customize_masked_climatology(self, climatology):  # {{{
        """
        Override this function to customize the climatology during the masking
        phase (before remapping)

        Parameters
        ----------
        climatology : ``xarray.Dataset```
            The MPAS climatology data set that has had a mask added but has
            not yet been remapped

        Returns
        -------
        climatology : ``xarray.Dataset```
            The same data set with any custom fields added or modifications
            made

        Authors
        -------
        Xylar Asay-Davis
        """

        return climatology  # }}}

    def customize_remapped_climatology(self, climatology):  # {{{
        """
        Override this function to customize the climatology after remapping

        Parameters
        ----------
        climatology : ``xarray.Dataset```
            The MPAS climatology data set that has been remapped

        Returns
        -------
        climatology : ``xarray.Dataset```
            The same data set with any custom fields added or modifications
            made

        Authors
        -------
        Xylar Asay-Davis
        """

        return climatology  # }}}

    def _setup_file_names(self):  # {{{
        """
        Create a dictionary of file names and directories for this climatology

        Authors
        -------
        Xylar Asay-Davis
        """

        config = self.config
        climatologyBaseDirectory = build_config_full_path(
            config, 'output', 'mpasClimatologySubdirectory')

        mpasMeshName = config.get('input', 'mpasMeshName')

        comparisonFullMeshNames = {}
        for comparisonGridName in self.comparisonGridNames:
            comparisonDescriptor = get_comparison_descriptor(
                    config, comparisonGridName)
            comparisonFullMeshNames[comparisonGridName] = \
                comparisonDescriptor.meshName

        keys = []
        for season in self.seasons:
            stage = 'masked'
            keys.append((season, stage))
            stage = 'remapped'
            for comparisonGridName in self.comparisonGridNames:
                keys.append((season, stage, comparisonGridName))

        self._outputDirs = {}
        self._outputFiles = {}

        for key in keys:
            season = key[0]
            stage = key[1]
            if stage == 'remapped':
                comparisonGridName = key[2]

            stageDirectory = '{}/{}'.format(climatologyBaseDirectory, stage)

            if stage == 'masked':
                directory = '{}/{}_{}'.format(
                        stageDirectory, self.climatologyName,
                        mpasMeshName)
            elif stage == 'remapped':
                directory = '{}/{}_{}_to_{}'.format(
                        stageDirectory,
                        self.climatologyName,
                        mpasMeshName,
                        comparisonFullMeshNames[comparisonGridName])

            make_directories(directory)

            monthValues = sorted(constants.monthDictionary[season])
            startMonth = monthValues[0]
            endMonth = monthValues[-1]

            suffix = '{:04d}{:02d}_{:04d}{:02d}_climo'.format(
                    self.mpasClimatologyTask.startYear, startMonth,
                    self.mpasClimatologyTask.endYear, endMonth)

            if season in constants.abrevMonthNames:
                season = '{:02d}'.format(monthValues[0])
            fileName = '{}/{}_{}_{}.nc'.format(
                    directory, self.mpasClimatologyTask.ncclimoModel, season,
                    suffix)

            self._outputDirs[key] = directory
            self._outputFiles[key] = fileName
        # }}}

    def _mask_climatologies(self, season, dsMask):  # {{{
        '''
        For each season, creates a masked version of the climatology

        Parameters
        ----------
        season : str
            The name of the season to be masked

        dsMask : ``xarray.Dataset`` object
            A data set (from the first input file) that can be used to
            determine the mask in MPAS output files.

        Author
        ------
        Xylar Asay-Davis
        '''

        climatologyFileName = self.mpasClimatologyTask.get_file_name(season)

        maskedClimatologyFileName = self.get_file_name(season, 'masked')

        if not os.path.exists(maskedClimatologyFileName):
            # slice and mask the data set
            climatology = xr.open_dataset(climatologyFileName)
            climatology = mpas_xarray.subset_variables(climatology,
                                                       self.variableList)
            iselValues = {'Time': 0}
            if self.iselValues is not None:
                iselValues.update(self.iselValues)
            # select only Time=0 and possibly only the desired vertical
            # slice
            climatology = climatology.isel(**iselValues)

            # add valid mask as a variable, useful for remapping later
            climatology['validMask'] = \
                xr.DataArray(numpy.ones(climatology.dims['nCells']),
                             dims=['nCells'])
            # mask the data set
            for variableName in self.variableList:
                climatology[variableName] = \
                    climatology[variableName].where(
                        dsMask[variableName] != self._fillValue)

            # customize (if this function has been overridden)
            climatology = self.customize_masked_climatology(climatology)

            write_netcdf(climatology, maskedClimatologyFileName)
        # }}}

    def _remap(self, inFileName, outFileName, remapper, comparisonGridName):
        # {{{
        """
        Performs remapping either using ``ncremap`` or the native python code,
        depending on the requested setting and the comparison grid

        Parameters
        ----------
        inFileName : str
            The name of the input file to be remapped.

        outFileName : str
            The name of the output file to which the remapped data set should
            be written.

        remapper : ``Remapper`` object
            A remapper that can be used to remap files or data sets to a
            comparison grid.

        comparisonGridNames : {'latlon', 'antarctic'}
            The name of the comparison grid to use for remapping.

        Authors
        -------
        Xylar Asay-Davis
        """
        if remapper.mappingFileName is None:
            # no remapping is needed
            return

        useNcremap = self.config.getboolean('climatology', 'useNcremap')

        if comparisonGridName != 'latlon':
            # ncremap doesn't support grids other than lat/lon
            useNcremap = False

        renormalizationThreshold = self.config.getfloat(
            'climatology', 'renormalizationThreshold')

        if useNcremap:
            remapper.remap_file(inFileName=inFileName,
                                outFileName=outFileName,
                                overwrite=True,
                                renormalize=renormalizationThreshold,
                                logger=self.logger)

            remappedClimatology = xr.open_dataset(outFileName)
            remappedClimatology.load()
            remappedClimatology.close()
        else:

            climatologyDataSet = xr.open_dataset(inFileName)

            remappedClimatology = remapper.remap(climatologyDataSet,
                                                 renormalizationThreshold)

        # customize (if this function has been overridden)
        remappedClimatology = self.customize_remapped_climatology(
                remappedClimatology)

        write_netcdf(remappedClimatology, outFileName)
        # }}}

    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
