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
    get_unmasked_mpas_climatology_file_name

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories, get_files_year_month

from mpas_analysis.shared.io import write_netcdf

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.mpas_xarray.mpas_xarray import subset_variables


class MpasClimatologyTask(AnalysisTask):  # {{{
    '''
    An analysis tasks for computing climatologies from output from the
    ``timeSeriesStatsMonthly`` analysis member.

    Attributes
    ----------

    variableList : dict of lists
        A dictionary with seasons as keys and a list of variable names in
        ``timeSeriesStatsMonthly`` to be included in the climatologies for each
        season in the values.

    allVariables : list of str
        A list of all available variable names in ``timeSeriesStatsMonthly``
        used to raise an exception when an unavailable variable is requested

    inputFiles : list of str
        A list of input files used to compute the climatologies.

    ncclimoModel : {'mpaso', 'mpascice'}
        The name of the component expected by ``ncclimo``

    startDate, endDate : str
        The start and end dates of the climatology as strings

    startYear, endYear : int
        The start and end years of the climatology
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, componentName, taskName=None):  # {{{
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

        taskName : str, optional
            the name of the task, defaults to mpasClimatology<ComponentName>
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.variableList = {}

        tags = ['climatology']

        if componentName == 'ocean':
            self.ncclimoModel = 'mpaso'
        elif componentName == 'seaIce':
            self.ncclimoModel = 'mpascice'
        else:
            raise ValueError('component {} is not supported by ncclimo.\n'
                             'Check with Charlie Zender and Xylar Asay-Davis\n'
                             'about getting it added'.format(componentName))

        if taskName is None:
            suffix = componentName[0].upper() + componentName[1:]
            taskName = 'mpasClimatology{}'.format(suffix)

        self.allVariables = None
        self.useNcclimo = config.getboolean('climatology', 'useNcclimo')

        # call the constructor from the base class (AnalysisTask)
        super(MpasClimatologyTask, self).__init__(
            config=config,
            taskName=taskName,
            componentName=componentName,
            tags=tags)

        if self.useNcclimo and \
                (config.get('execute', 'ncclimoParallelMode') == 'bck'):
            self.subprocessCount = 12

        parallelTaskCount = config.getint('execute', 'parallelTaskCount')
        if not self.useNcclimo:
            self.subprocessCount = parallelTaskCount

        # }}}

    def add_variables(self, variableList, seasons=None):  # {{{
        '''
        Add one or more variables and optionally one or more seasons for which
        to compute climatologies.

        Parameters
        ----------
        variableList : list of str
            A list of variable names in ``timeSeriesStatsMonthly`` to be
            included in the climatologies

        seasons : list of str, optional
            A list of seasons (keys in ``shared.constants.monthDictionary``)
            to be computed or ``None`` if only monthly
            climatologies are needed.

        Raises
        ------
        ValueError
            if this funciton is called before this task has been set up (so
            the list of available variables has not yet been set) or if one
            or more of the requested variables is not available in the
            ``timeSeriesStatsMonthly`` output.
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
                    '{} is not available in timeSeriesStatsMonthly '
                    'output:\n{}'.format(variable, self.allVariables))

            for season in seasons:
                if season not in self.variableList:
                    self.variableList[season] = []
                if variable not in self.variableList[season]:
                    self.variableList[season].append(variable)
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

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        self.startYear, self.endYear = self.get_start_and_end()

        self.startDate = '{:04d}-01-01_00:00:00'.format(self.startYear)
        self.endDate = '{:04d}-12-31_23:59:59'.format(self.endYear)

        # get a list of timeSeriesSta output files from the streams file,
        # reading only those that are between the start and end dates
        streamName = 'timeSeriesStatsMonthlyOutput'
        self.inputFiles = self.historyStreams.readpath(
            streamName, startDate=self.startDate, endDate=self.endDate,
            calendar=self.calendar)

        if len(self.inputFiles) == 0:
            raise IOError('No files were found in stream {} between {} and '
                          '{}.'.format(streamName, self.startDate,
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

        self.logger.info('\nComputing MPAS climatologies from files:\n'
                         '    {} through\n    {}'.format(
                             os.path.basename(self.inputFiles[0]),
                             os.path.basename(self.inputFiles[-1])))

        if self.useNcclimo:
            seasonsToCheck = list(constants.abrevMonthNames)
        else:
            seasonsToCheck = []

        for season in self.variableList:
            if season not in seasonsToCheck:
                seasonsToCheck.append(season)

        allExist = True
        for season in seasonsToCheck:

            climatologyFileName = self.get_file_name(season)
            climatologyDirectory = get_unmasked_mpas_climatology_directory(
                self.config)

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
            if self.useNcclimo:
                self._compute_climatologies_with_ncclimo(
                    inDirectory=self.symlinkDirectory,
                    outDirectory=climatologyDirectory)
            else:
                self._compute_climatologies_with_xarray(
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
                                                       self.componentName)

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
                                             'timeSeriesStatsMonthlyOutput')

        climatologyBaseDirectory = build_config_full_path(
            config, 'output', 'mpasClimatologySubdirectory')

        symlinkDirectory = '{}/source_symlinks'.format(
            climatologyBaseDirectory)

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

        fileNames = sorted(self.inputFiles)
        years, months = get_files_year_month(fileNames,
                                             self.historyStreams,
                                             'timeSeriesStatsMonthlyOutput')

        chunkSize = self.config.getint('input', 'maxChunkSize')

        with xarray.open_mfdataset(self.inputFiles, concat_dim='Time',
                                   chunks={'nCells': chunkSize},
                                   decode_cf=False, decode_times=False) as ds:

            ds.coords['year'] = ('Time', years)
            ds.coords['month'] = ('Time', months)
            with dask.config.set(scheduler='threads'):
                for month in range(1, 13):
                    monthName = constants.abrevMonthNames[month-1]
                    climatologyFileName = self.get_file_name(monthName)
                    self.logger.info('computing climatology {}'.format(
                        os.path.basename(climatologyFileName)))
                    variableList = []
                    for season in self.variableList:
                        monthValues = constants.monthDictionary[season]
                        if month in monthValues:
                            variableList.extend(self.variableList[season])
                    # just the unique variables
                    variableList = sorted(list(set(variableList)))

                    if len(variableList) == 0:
                        continue

                    dsMonth = subset_variables(ds, variableList)
                    dsMonth = dsMonth.where(dsMonth.month == month, drop=True)
                    dsMonth = dsMonth.mean(dim='Time')
                    dsMonth.compute(num_workers=self.subprocessCount)
                    write_netcdf(dsMonth, climatologyFileName)

        for season in self.variableList:
            outFileName = self.get_file_name(season=season)
            self.logger.info('computing climatology {}'.format(
                os.path.basename(outFileName)))
            fileNames = []
            weights = []
            for month in constants.monthDictionary[season]:
                monthName = constants.abrevMonthNames[month-1]
                fileNames.append(self.get_file_name(season=monthName))
                weights.append(constants.daysInMonth[month-1])
            with xarray.open_mfdataset(fileNames, concat_dim='weight',
                                       chunks={'nCells': chunkSize},
                                       decode_cf=False, decode_times=False) \
                    as ds:
                ds = subset_variables(ds, self.variableList[season])
                ds.coords['weight'] = ('weight', weights)
                ds = ((ds.weight*ds).sum(dim='weight') /
                      ds.weight.sum(dim='weight'))
                ds.compute(num_workers=self.subprocessCount)
                write_netcdf(ds, outFileName)

        # }}}
    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
