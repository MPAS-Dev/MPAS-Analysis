import xarray
import os
import warnings
import subprocess
from distutils.spawn import find_executable

from ..analysis_task import AnalysisTask

from ..constants import constants

from ..io.utility import build_config_full_path, make_directories


class MpasClimatologyTask(AnalysisTask):  # {{{
    '''
    An analysis tasks for computing climatologies from output from the
    ``timeSeriesStatsMonthly`` analysis member.

    Attributes
    ----------

    variableList : list of str
        A list of variable names in ``timeSeriesStatsMonthly`` to be
        included in the climatologies

    seasons : list of str
        A list of seasons (keys in ``shared.constants.monthDictionary``)
        over which the climatology should be computed or ['none'] if only
        monthly climatologies are needed.

    inputFiles : list of str
        A list of input files used to compute the climatologies.

    ncclimoModel : {'mpaso', 'mpascice'}
        The name of the component expected by ``ncclimo``

    startDate, endDate : str
        The start and end dates of the climatology as strings

    startYear, endYear : int
        The start and end years of the climatology

    Authors
    -------
    Xylar Asay-Davis
    '''

    def __init__(self, config, componentName):  # {{{
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

        Authors
        -------
        Xylar Asay-Davis
        '''
        self.variableList = []
        self.seasons = []

        tags = ['climatology']

        if componentName == 'ocean':
            self.ncclimoModel = 'mpaso'
        elif componentName == 'seaIce':
            self.ncclimoModel = 'mpascice'
        else:
            raise ValueError('component {} is not supported by ncclimo.\n'
                             'Check with Charlie Zender and Xylar Asay-Davis\n'
                             'about getting it added'.format(componentName))

        suffix = componentName[0].upper() + componentName[1:]
        taskName = 'mpasClimatology{}'.format(suffix)

        # call the constructor from the base class (AnalysisTask)
        super(MpasClimatologyTask, self).__init__(
            config=config,
            taskName=taskName,
            componentName=componentName,
            tags=tags)

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
            to be computed or ['none'] (not ``None``) if only monthly
            climatologies are needed.

        Authors
        -------
        Xylar Asay-Davis
        '''

        for variable in variableList:
            if variable not in self.variableList:
                self.variableList.append(variable)

        if seasons is not None:
            for season in seasons:
                if season not in self.seasons:
                    self.seasons.append(season)

        self._setup_file_names()

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
        super(MpasClimatologyTask, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        # get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        startDate = self.config.get('climatology', 'startDate')
        endDate = self.config.get('climatology', 'endDate')
        streamName = 'timeSeriesStatsMonthlyOutput'
        self.inputFiles = self.historyStreams.readpath(
                streamName, startDate=startDate, endDate=endDate,
                calendar=self.calendar)

        if len(self.inputFiles) == 0:
            raise IOError('No files were found in stream {} between {} and '
                          '{}.'.format(streamName, startDate, endDate))

        self._update_climatology_bounds_from_file_names()

        # }}}

    def run_task(self):  # {{{
        '''
        Compute the requested climatologies

        Authors
        -------
        Xylar Asay-Davis
        '''

        if len(self.variableList) == 0:
            # nothing to do
            return

        self.logger.info('\nComputing MPAS climatologies from files:\n'
                         '    {} through\n    {}'.format(
                                 os.path.basename(self.inputFiles[0]),
                                 os.path.basename(self.inputFiles[-1])))

        if self.seasons[0] is 'none':
            seasonsToCheck = ['{:02d}'.format(month) for month in range(1, 13)]
        else:
            seasonsToCheck = self.seasons

        allExist = True
        for season in seasonsToCheck:

            climatologyFileName, climatologyDirectory = \
                self.get_file_name(season, returnDir=True)

            if not os.path.exists(climatologyFileName):
                allExist = False
                break

        if allExist:
            # make sure all the necessary variables are also present
            ds = xarray.open_dataset(self.get_file_name(seasonsToCheck[0],
                                                        returnDir=False))

            for variableName in self.variableList:
                if variableName not in ds.variables:
                    allExist = False
                    break

        if not allExist:
            self._compute_climatologies_with_ncclimo(
                    inDirectory=self.historyDirectory,
                    outDirectory=climatologyDirectory)

        # }}}

    def get_file_name(self, season, returnDir=False):  # {{{
        """
        Given config options, the name of a field and a string identifying the
        months in a seasonal climatology, returns the full path for MPAS
        climatology files before and after remapping.

        Parameters
        ----------
        season : str
            One of the seasons in ``constants.monthDictionary``

        mpasMeshName : str
            The name of the MPAS mesh

        returnDir : bool, optional
            Return the directory as well

        Returns
        -------
        fileName : str
            The path to the climatology file for the specified season.

        Authors
        -------
        Xylar Asay-Davis
        """

        fileName = self._outputFiles[season]

        if returnDir:
            directory = self._outputDirs[season]
            return fileName, directory
        else:
            return fileName  # }}}

    def _update_climatology_bounds_from_file_names(self):  # {{{
        """
        Update the start and end years and dates for climatologies based on the
        years actually available in the list of files.

        Authors
        -------
        Xylar Asay-Davis
        """

        config = self.config

        requestedStartYear = config.getint('climatology', 'startYear')
        requestedEndYear = config.getint('climatology', 'endYear')

        dates = sorted([fileName[-13:-6] for fileName in self.inputFiles])
        years = [int(date[0:4]) for date in dates]
        months = [int(date[5:7]) for date in dates]

        # search for the start of the first full year
        firstIndex = 0
        while(firstIndex < len(years) and months[firstIndex] != 1):
            firstIndex += 1
        startYear = years[firstIndex]

        # search for the end of the last full year
        lastIndex = len(years)-1
        while(lastIndex >= 0 and months[lastIndex] != 12):
            lastIndex -= 1
        endYear = years[lastIndex]

        if startYear != requestedStartYear or endYear != requestedEndYear:
            message = "climatology start and/or end year different from " \
                      "requested\n" \
                      "requestd: {:04d}-{:04d}\n" \
                      "actual:   {:04d}-{:04d}\n".format(requestedStartYear,
                                                         requestedEndYear,
                                                         startYear,
                                                         endYear)
            warnings.warn(message)
            config.set('climatology', 'startYear', str(startYear))
            config.set('climatology', 'endYear', str(endYear))

            startDate = '{:04d}-01-01_00:00:00'.format(startYear)
            config.set('climatology', 'startDate', startDate)
            endDate = '{:04d}-12-31_23:59:59'.format(endYear)
            config.set('climatology', 'endDate', endDate)

        else:
            startDate = config.get('climatology', 'startDate')
            endDate = config.get('climatology', 'endDate')

        self.startDate = startDate
        self.endDate = endDate
        self.startYear = startYear
        self.endYear = endYear

        # }}}

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

        self._outputDirs = {}
        self._outputFiles = {}

        directory = '{}/unmasked_{}'.format(climatologyBaseDirectory,
                                            mpasMeshName)

        make_directories(directory)
        for season in self.seasons:
            monthValues = sorted(constants.monthDictionary[season])
            startMonth = monthValues[0]
            endMonth = monthValues[-1]

            suffix = '{:04d}{:02d}_{:04d}{:02d}_climo'.format(
                    self.startYear, startMonth, self.endYear, endMonth)

            if season in constants.abrevMonthNames:
                season = '{:02d}'.format(monthValues[0])
            fileName = '{}/{}_{}_{}.nc'.format(directory, self.ncclimoModel,
                                               season, suffix)

            self._outputDirs[season] = directory
            self._outputFiles[season] = fileName

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

        Author
        ------
        Xylar Asay-Davis
        '''

        if find_executable('ncclimo') is None:
            raise OSError('ncclimo not found. Make sure the latest nco '
                          'package is installed: \n'
                          'conda install nco\n'
                          'Note: this presumes use of the conda-forge '
                          'channel.')

        parallelMode = self.config.get('execute', 'ncclimoParallelMode')

        args = ['ncclimo',
                '--clm_md=mth',
                '-a', 'sdd',
                '-m', self.ncclimoModel,
                '-p', parallelMode,
                '-v', ','.join(self.variableList),
                '--seasons={}'.format(','.join(self.seasons)),
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
            self.logger.info(stdout)
        if stderr:
            for line in stderr.split('\n'):
                self.logger.error(line)

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode,
                                                ' '.join(args))

        # }}}
    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
