# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray
import os
import subprocess
from distutils.spawn import find_executable
import dask

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.climatology.climatology import \
    get_unmasked_mpas_climatology_directory, \
    get_unmasked_mpas_climatology_file_name, \
    get_climatology_op_directory

from mpas_analysis.shared.io.utility import make_directories, \
    get_files_year_month

from mpas_analysis.shared.io import write_netcdf

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.mpas_xarray.mpas_xarray import subset_variables


class MpasClimatologyTask(AnalysisTask):  # {{{
    '''
    An analysis tasks for computing climatologies from output from the
    ``timeSeriesStatsMonthly*`` analysis members.

    Attributes
    ----------

    variableList : dict of lists
        A dictionary with seasons as keys and a list of variable names in
        the stream to be included in the climatologies for each season in the
        values.

    allVariables : list of str
        A list of all available variable names in the stream used to raise an
        exception when an unavailable variable is requested

    inputFiles : list of str
        A list of input files used to compute the climatologies.

    ncclimoModel : {'mpaso', 'mpascice'}
        The name of the component expected by ``ncclimo``

    startDate, endDate : str
        The start and end dates of the climatology as strings

    startYear, endYear : int
        The start and end years of the climatology

    seasonSubtasks : dict
        If using xarray to compute climatologies, a dictionary of subtasks, one
        for each possible season

    op : {'avg', 'min', 'max'}
         operator for monthly stats

    streamName : str
        The name of the stream to read from, one of
        ``timeSeriesStatsMonthlyOutput``,
        ``timeSeriesStatsMonthlyMinOutput``,
        ``timeSeriesStatsMonthlyMaxOutput``
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, componentName, taskName=None,
                 op='avg'):  # {{{
        '''
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        config : ``MpasAnalysisConfigParser``
            Contains configuration options

        componentName : {'ocean', 'seaIce'}
            The name of the component (same as the folder where the task
            resides)

        op : {'avg', 'min', 'max'}, optioinal
             operator for monthly stats

        taskName : str, optional
            the name of the task, defaults to mpasClimatology<ComponentName>
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.variableList = {}

        self.op = op
        if op == 'avg':
            self.streamName = 'timeSeriesStatsMonthlyOutput'
        elif op == 'min':
            self.streamName = 'timeSeriesStatsMonthlyMinOutput'
        elif op == 'max':
            self.streamName = 'timeSeriesStatsMonthlyMaxOutput'
        else:
            raise ValueError('Unexpected monthly stats operator {}'.format(op))

        tags = ['climatology', op]

        if componentName == 'ocean':
            self.ncclimoModel = 'mpaso'
        elif componentName == 'seaIce':
            self.ncclimoModel = 'mpascice'
        else:
            raise ValueError('component {} is not supported by ncclimo.\n'
                             'Check with Charlie Zender and Xylar Asay-Davis\n'
                             'about getting it added'.format(componentName))

        if taskName is None:
            suffix = componentName[0].upper() + componentName[1:] + \
                op[0].upper() + op[1:]
            taskName = 'mpasClimatology{}'.format(suffix)

        self.allVariables = None
        self.useNcclimo = config.getboolean('climatology', 'useNcclimo')

        # call the constructor from the base class (AnalysisTask)
        super(MpasClimatologyTask, self).__init__(
            config=config,
            taskName=taskName,
            componentName=componentName,
            tags=tags)

        ncclimoParallelMode = config.get('execute', 'ncclimoParallelMode')
        if self.useNcclimo and ncclimoParallelMode in ['bck', 'mpi']:
            self.subprocessCount = 12

        self.seasonSubtasks = {}

        parallelTaskCount = config.getint('execute', 'parallelTaskCount')
        if not self.useNcclimo:
            self.subprocessCount = 1

            # setup one subtask for each possible season that could be added
            for season in constants.monthDictionary:
                if season in constants.abrevMonthNames:
                    subprocessCount = max(parallelTaskCount // 12, 1)
                else:
                    # something to experiment with
                    subprocessCount = 2
                self.seasonSubtasks[season] = MpasClimatologySeasonSubtask(
                    self, season, subprocessCount=subprocessCount)
                self.add_subtask(self.seasonSubtasks[season])

            # make sure each season runs after the months that make up that
            # season
            for season in constants.monthDictionary:
                if season in constants.abrevMonthNames:
                    continue
                monthValues = constants.monthDictionary[season]
                monthNames = [constants.abrevMonthNames[month-1] for month in
                              monthValues]
                for monthName in monthNames:
                    self.seasonSubtasks[season].run_after(
                            self.seasonSubtasks[monthName])
        # }}}

    def add_variables(self, variableList, seasons=None):  # {{{
        '''
        Add one or more variables and optionally one or more seasons for which
        to compute climatologies.

        Parameters
        ----------
        variableList : list of str
            A list of variable names in the stream to be included in the
            climatologies

        seasons : list of str, optional
            A list of seasons (keys in ``shared.constants.monthDictionary``)
            to be computed or ``None`` if only monthly
            climatologies are needed.

        Raises
        ------
        ValueError
            if this funciton is called before this task has been set up (so
            the list of available variables has not yet been set) or if one
            or more of the requested variables is not available in the stream.
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        if self.allVariables is None:
            raise ValueError('add_variables() can only be called after '
                             'setup_and_check() in MpasClimatologyTask.\n'
                             'Presumably tasks were added in the wrong order '
                             'or add_variables() is being called in the wrong '
                             'place.')

        if seasons is None:
            seasons = list(constants.abrevMonthNames)

        for variable in variableList:
            if variable not in self.allVariables:
                raise ValueError(
                    '{} is not available in {} output:\n{}'.format(
                        variable, self.streamName, self.allVariables))

            for season in seasons:
                if season not in self.variableList:
                    self.variableList[season] = []
                if variable not in self.variableList[season]:
                    self.variableList[season].append(variable)

        # add variables to individual months as well, since those will
        # be computed first
        for season in seasons:
            if season not in constants.abrevMonthNames:
                monthValues = constants.monthDictionary[season]
                monthNames = [constants.abrevMonthNames[month-1] for month in
                              monthValues]
                self.add_variables(variableList, seasons=monthNames)

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
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(MpasClimatologyTask, self).setup_and_check()

        if self.op == 'avg':
            self.check_analysis_enabled(
                analysisOptionName='config_am_timeseriesstatsmonthly_enable',
                raiseException=True)
        elif self.op == 'min':
            self.check_analysis_enabled(
                analysisOptionName='config_AM_timeSeriesStatsMonthlyMin_enable',
                raiseException=True)
        elif self.op == 'max':
            self.check_analysis_enabled(
                analysisOptionName='config_AM_timeSeriesStatsMonthlyMax_enable',
                raiseException=True)

        self.startYear, self.endYear = self.get_start_and_end()

        self.startDate = '{:04d}-01-01_00:00:00'.format(self.startYear)
        self.endDate = '{:04d}-12-31_23:59:59'.format(self.endYear)

        # get a list of timeSeriesSta output files from the streams file,
        # reading only those that are between the start and end dates
        self.inputFiles = self.historyStreams.readpath(
            self.streamName, startDate=self.startDate, endDate=self.endDate,
            calendar=self.calendar)

        if len(self.inputFiles) == 0:
            raise IOError('No files were found in stream {} between {} and '
                          '{}.'.format(self.streamName, self.startDate,
                                       self.endDate))

        self.symlinkDirectory = self._create_symlinks()

        with xarray.open_dataset(self.inputFiles[0]) as ds:
            self.allVariables = list(ds.data_vars.keys())

        # }}}

    def run_task(self):  # {{{
        '''
        Compute the requested climatologies
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        if len(self.variableList.keys()) == 0:
            # nothing to do
            return

        if not self.useNcclimo:
            # subtasks will take care of it, so nothing to do
            return

        self.logger.info('\nComputing MPAS climatologies from files:\n'
                         '    {} through\n    {}'.format(
                             os.path.basename(self.inputFiles[0]),
                             os.path.basename(self.inputFiles[-1])))

        seasonsToCheck = list(constants.abrevMonthNames)

        for season in self.variableList:
            if season not in seasonsToCheck:
                seasonsToCheck.append(season)

        allExist = True
        for season in seasonsToCheck:

            climatologyFileName = self.get_file_name(season)
            climatologyDirectory = get_unmasked_mpas_climatology_directory(
                self.config, self.op)

            if not os.path.exists(climatologyFileName):
                allExist = False
                break

        if allExist:
            for season in seasonsToCheck:
                if season not in self.variableList:
                    continue
                # make sure all the necessary variables are also present
                with xarray.open_dataset(self.get_file_name(season)) as ds:
                    for variableName in self.variableList[season]:
                        if variableName not in ds.variables:
                            allExist = False
                            break

        if not allExist:
            self._compute_climatologies_with_ncclimo(
                inDirectory=self.symlinkDirectory,
                outDirectory=climatologyDirectory)

        # }}}

    def get_start_and_end(self):  # {{{
        """
        Get the start and end years and dates for the climatology.  This
        function is provided to allow a custom method for setting the start
        and end years of the climatology.  By default, they are read from the
        climatology section of the config file

        Returns
        -------
        startYear, endYear : int
           The start and end years of the climatology

        """
        # Authors
        # -------
        # Xylar Asay-Davis

        startYear = self.config.getint('climatology', 'startYear')
        endYear = self.config.getint('climatology', 'endYear')

        return startYear, endYear

        # }}}

    def get_file_name(self, season):  # {{{
        """

        Returns the full path for MPAS climatology file produced by ncclimo.

        Parameters
        ----------
        season : str
            One of the seasons in ``constants.monthDictionary``

        Returns
        -------
        fileName : str
            The path to the climatology file for the specified season.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        return get_unmasked_mpas_climatology_file_name(self.config, season,
                                                       self.componentName,
                                                       self.op)

        # }}}

    def _create_symlinks(self):  # {{{
        """
        Create symlinks to monthly mean files so they have the expected file
        naming convention for ncclimo.

        Returns
        -------
        symlinkDirectory : str
            The path to the symlinks created for each timeSeriesStatsMonthly
            input file
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config

        fileNames = sorted(self.inputFiles)
        years, months = get_files_year_month(fileNames,
                                             self.historyStreams,
                                             self.streamName)

        climatologyOpDirectory = get_climatology_op_directory(config, self.op)

        symlinkDirectory = '{}/source_symlinks'.format(
            climatologyOpDirectory)

        make_directories(symlinkDirectory)

        for inFileName, year, month in zip(fileNames, years, months):
            outFileName = '{}/{}.hist.am.timeSeriesStatsMonthly.{:04d}-' \
                '{:02d}-01.nc'.format(symlinkDirectory, self.ncclimoModel,
                                      year, month)

            if not os.path.exists(outFileName):
                os.symlink(inFileName, outFileName)

        return symlinkDirectory

        # }}}

    def _compute_climatologies_with_ncclimo(self, inDirectory, outDirectory,
                                            remapper=None,
                                            remappedDirectory=None):  # {{{
        '''
        Uses ncclimo to compute monthly, seasonal and/or annual climatologies.

        Parameters
        ----------
        inDirectory : str
            The run directory containing timeSeriesStatsMonthly output

        outDirectory : str
            The output directory where climatologies will be written

        remapper : ``shared.intrpolation.Remapper`` object, optional
            If present, a remapper that defines the source and desitnation
            grids for remapping the climatologies.

        remappedDirectory : str, optional
            If present, the path where remapped climatologies should be
            written. By default, remapped files are stored in the same
            directory as the climatologies on the source grid.  Has no effect
            if ``remapper`` is ``None``.

        Raises
        ------
        OSError
            If ``ncclimo`` is not in the system path.
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        if find_executable('ncclimo') is None:
            raise OSError('ncclimo not found. Make sure the latest nco '
                          'package is installed: \n'
                          'conda install nco\n'
                          'Note: this presumes use of the conda-forge '
                          'channel.')

        parallelMode = self.config.get('execute', 'ncclimoParallelMode')

        seasons = [season for season in self.variableList
                   if season not in constants.abrevMonthNames]

        variableList = []
        for season in self.variableList:
            variableList.extend(self.variableList[season])

        # include each variable only once
        variableList = sorted(list(set(variableList)))

        if len(seasons) == 0:
            seasons = ['none']

        args = ['ncclimo',
                '-4',
                '--clm_md=mth',
                '-a', 'sdd',
                '-m', self.ncclimoModel,
                '-p', parallelMode,
                '-v', ','.join(variableList),
                '--seasons={}'.format(','.join(seasons)),
                '-s', '{:04d}'.format(self.startYear),
                '-e', '{:04d}'.format(self.endYear),
                '-i', inDirectory,
                '-o', outDirectory]

        if remapper is not None:
            args.extend(['-r', remapper.mappingFileName])
            if remappedDirectory is not None:
                args.extend(['-O', remappedDirectory])

        self.logger.info('running: {}'.format(' '.join(args)))
        for handler in self.logger.handlers:
            handler.flush()

        # set an environment variable to make sure we're not using czender's
        # local version of NCO instead of one we have intentionally loaded
        env = os.environ.copy()
        env['NCO_PATH_OVERRIDE'] = 'No'

        process = subprocess.Popen(args, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, env=env)
        stdout, stderr = process.communicate()

        if stdout:
            stdout = stdout.decode('utf-8')
            for line in stdout.split('\n'):
                self.logger.info(line)
        if stderr:
            stderr = stderr.decode('utf-8')
            for line in stderr.split('\n'):
                self.logger.error(line)

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode,
                                                ' '.join(args))

        # }}}
    # }}}


class MpasClimatologySeasonSubtask(AnalysisTask):  # {{{
    '''
    An analysis subtasks for computing climatologies from output from the
    ``timeSeriesStatsMonthly`` analysis member for a single month or season.

    Attributes
    ----------

    season : str
        The season of the climatology

    parentTask : ``MpasClimatologyTask``
        The task that this subtask belongs to.
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, season, subtaskName=None,
                 subprocessCount=1):  # {{{
        '''
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        parentTask : ``MpasClimatologyTask``
            The task that this subtask belongs to.

        season : str
            A keys in ``shared.constants.monthDictionary``

        subtaskName : str, optional
            the name of the subtask, defaults to season

        subprocessCount : int, optional
            The number of subprocesses (dask threads) this task is allowed to
            use
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.season = season

        if subtaskName is None:
            subtaskName = season

        self.subprocessCount = subprocessCount
        self.parentTask = parentTask

        # call the constructor from the base class (AnalysisTask)
        super(MpasClimatologySeasonSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName=subtaskName)

        # }}}

    def run_task(self):  # {{{
        '''
        Compute the requested climatologies
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        season = self.season
        parentTask = self.parentTask
        if season not in parentTask.variableList:
            # nothing to do
            return

        variableList = parentTask.variableList[season]

        if len(variableList) == 0:
            # nothing to do
            return

        self.logger.info('\nComputing MPAS climatology from files:\n'
                         '    {} through\n    {}'.format(
                             os.path.basename(parentTask.inputFiles[0]),
                             os.path.basename(parentTask.inputFiles[-1])))

        climatologyFileName = parentTask.get_file_name(season)
        climatologyDirectory = get_unmasked_mpas_climatology_directory(
            self.config)

        allExist = False
        if os.path.exists(climatologyFileName):
            allExist = True
            # make sure all the necessary variables are also present
            with xarray.open_dataset(climatologyFileName) as ds:
                for variableName in variableList:
                    if variableName not in ds.variables:
                        allExist = False
                        break

        if not allExist:
            with dask.config.set(scheduler='threads'):
                self._compute_climatologies_with_xarray(
                    inDirectory=parentTask.symlinkDirectory,
                    outDirectory=climatologyDirectory)

        # }}}

    def _compute_climatologies_with_xarray(self, inDirectory, outDirectory):
        # {{{
        '''
        Uses xarray to compute seasonal and/or annual climatologies.

        Parameters
        ----------
        inDirectory : str
            The run directory containing timeSeriesStatsMonthly output

        outDirectory : str
            The output directory where climatologies will be written
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        def _preprocess(ds):
            # drop unused variables during preprocessing because only the
            # variables we want are guaranteed to be in all the files
            return subset_variables(ds, variableList)

        season = self.season
        parentTask = self.parentTask
        variableList = parentTask.variableList[season]

        chunkSize = self.config.getint('input', 'maxChunkSize')

        if season in constants.abrevMonthNames:
            # this is an individual month, so create a climatology from
            # timeSeriesStatsMonthlyOutput

            fileNames = sorted(parentTask.inputFiles)
            years, months = get_files_year_month(
                fileNames,  self.historyStreams,
                'timeSeriesStatsMonthlyOutput')

            with xarray.open_mfdataset(parentTask.inputFiles,
                                       concat_dim='Time',
                                       chunks={'nCells': chunkSize},
                                       decode_cf=False, decode_times=False,
                                       preprocess=_preprocess) as ds:

                ds.coords['year'] = ('Time', years)
                ds.coords['month'] = ('Time', months)
                month = constants.abrevMonthNames.index(season) + 1
                climatologyFileName = parentTask.get_file_name(season)
                self.logger.info('computing climatology {}'.format(
                    os.path.basename(climatologyFileName)))

                ds = ds.where(ds.month == month, drop=True)
                ds = ds.mean(dim='Time')
                ds.compute(num_workers=self.subprocessCount)
                write_netcdf(ds, climatologyFileName)
        else:
            outFileName = parentTask.get_file_name(season=season)
            self.logger.info('computing climatology {}'.format(
                os.path.basename(outFileName)))
            fileNames = []
            weights = []
            for month in constants.monthDictionary[season]:
                monthName = constants.abrevMonthNames[month-1]
                fileNames.append(parentTask.get_file_name(season=monthName))
                weights.append(constants.daysInMonth[month-1])

            with xarray.open_mfdataset(fileNames, concat_dim='weight',
                                       chunks={'nCells': chunkSize},
                                       decode_cf=False, decode_times=False,
                                       preprocess=_preprocess) as ds:
                ds.coords['weight'] = ('weight', weights)
                ds = ((ds.weight*ds).sum(dim='weight') /
                      ds.weight.sum(dim='weight'))
                ds.compute(num_workers=self.subprocessCount)
                write_netcdf(ds, outFileName)

        # }}}
    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
