'''
Defines the base class for analysis tasks.

Authors
-------
Xylar Asay-Davis
'''

import warnings
import numpy

from .timekeeping.utility import days_to_datetime, string_to_days_since_date

from .io import NameList, StreamsFile, build_config_full_path, \
    make_directories

from .variable_namelist_stream_maps.ocean_maps import oceanNamelistMap, \
    oceanStreamMap, oceanVariableMap

from .variable_namelist_stream_maps.sea_ice_maps import seaIceNamelistMap, \
    seaIceStreamMap, seaIceVariableMap


class AnalysisTask(object):  # {{{
    '''
    The base class for analysis tasks.

    Attributes
    ----------
    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    taskName : str
        The name of the task, typically the same as the class name except
        starting with lowercase (e.g. 'myTask' for class 'MyTask')

    componentName :  {'ocean', 'seaIce'}
        The name of the component (same as the folder where the task
        resides)

    tags :  list of str
        Tags used to describe the task (e.g. 'timeSeries', 'climatology',
        horizontalMap', 'index', 'transect').  These are used to determine
        which tasks are generated (e.g. 'all_transect' or 'no_climatology'
        in the 'generate' flags)

    prerequisiteTasks :  list of str, optional
        Names of tasks that must complete before this task can run.
        Typically, this will include one or more tasks of the form
        ``cache<Component><StreamName>Times``, e.g.
        ``cacheOceanTimeSeriesStatsTimes``

    runDirectory : str
        the base input directory for namelists, streams files and restart files

    historyDirectory : str
        the base input directory for history files

    plotsDirectory : str
        the directory for writing plots (which is also created if it doesn't
        exist)

    namelist : ``NameList`` object
        the namelist reader

    runStreams : ``StreamsFile`` object
        the streams file reader for streams in the run directory (e.g. restart
        files)

    historyStreams : ``StreamsFile`` object
        the streams file reader for streams in the history directory (most
        streams other than restart files)

    calendar : {'gregorian', 'gregorian_noleap'}
        the name of the calendar

    namelistMap : dict
        A map between names of namelist options used by MPAS-Analysis and
        those in various MPAS versions

    streamMap  : dict
        a map between names of streams used by MPAS-Analysis and those in
        various MPAS versions

    variableMap : dict
        a map between names of variables within streams used by MPAS-Analysis
        and those in various MPAS versions

    Authors
    -------
    Xylar Asay-Davis

    '''
    @classmethod
    def create_tasks(cls, config):  # {{{
        """
        Construct a single analysis task.  This function can be replaced in
        child classes so that a list of 2 or more tasks is returned, allowing
        a task to also provide its own prerequisite tasks.

        Parameters
        ----------
        config : MpasAnalysisConfigParser object
            Contains configuration options

        Authors
        -------
        Xylar Asay-Davis
        """

        tasks = [cls(config=config)]
        return tasks  # }}}

    def __init__(self, config, taskName, componentName, tags=[],
                 prerequisiteTasks=None):  # {{{
        '''
        Construct the analysis task.

        Individual tasks (children classes of this base class) should first
        call this method to perform basic initialization, then, define the
        `taskName`, `componentName` and list of `tags` for the task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        taskName :  str
            The name of the task, typically the same as the class name except
            starting with lowercase (e.g. 'myTask' for class 'MyTask')

        componentName :  {'ocean', 'seaIce'}
            The name of the component (same as the folder where the task
            resides)

        tags :  list of str, optional
            Tags used to describe the task (e.g. 'timeSeries', 'climatology',
            horizontalMap', 'index', 'transect').  These are used to determine
            which tasks are generated (e.g. 'all_transect' or 'no_climatology'
            in the 'generate' flags)

        prerequisiteTasks :  list of str, optional
            Names of tasks that must complete before this task can run.
            Typically, this will include one or more tasks of the form
            ``cache<Component><StreamName>Times``, e.g.
            ``cacheOceanTimeSeriesStatsTimes``

        Authors
        -------
        Xylar Asay-Davis
        '''
        self.config = config
        self.taskName = taskName
        self.componentName = componentName
        self.tags = tags
        self.prerequisiteTasks = prerequisiteTasks  # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis (e.g. reading namelists and
        streams files).

        Individual tasks (children classes of this base class) should first
        call this method to perform basic setup, then, check whether the
        configuration is correct for a given analysis and perform additional,
        analysis-specific setup.  For example, this function could check if
        necessary observations and other data files are found, then, determine
        the list of files to be read when the analysis is run.

        If the task includes ``climatology``, ``timeSeries`` or ``index`` tags,
        ``startDate`` and ``endDate`` config options are computed from
        ``startYear`` and ``endYear``config options.

        After this call, the following attributes are set.

        Attributes
        ----------
        runDirectory : str
            the base input directory for namelists, streams files and restart
            files

        historyDirectory : str
            the base input directory for history files

        plotsDirectory : str
            the directory for writing plots (which is also created if it
            doesn't exist)

        namelist : ``NameList`` object
            the namelist reader

        runStreams : ``StreamsFile`` object
            the streams file reader for streams in the run directory (e.g.
            restart files)

        historyStreams : ``StreamsFile`` object
            the streams file reader for streams in the history directory (most
            streams other than restart files)

        calendar : {'gregorian', 'gregorian_noleap'}
            the name of the calendar

        namelistMap : dict
            A map between names of namelist options used by MPAS-Analysis and
            those in various MPAS versions

        streamMap  : dict
            a map between names of streams used by MPAS-Analysis and those in
            various MPAS versions

        variableMap : dict
            a map between names of variables within streams used by
            MPAS-Analysis and those in various MPAS versions

        Authors
        -------
        Xylar Asay-Davis
        '''

        # read parameters from config file
        # the run directory contains the restart files
        self.runDirectory = build_config_full_path(self.config, 'input',
                                                   'runSubdirectory')
        # if the history directory exists, use it; if not, fall back on
        # runDirectory
        self.historyDirectory = build_config_full_path(
            self.config, 'input',
            '{}HistorySubdirectory'.format(self.componentName),
            defaultPath=self.runDirectory)

        self.plotsDirectory = build_config_full_path(self.config, 'output',
                                                     'plotsSubdirectory')
        namelistFileName = build_config_full_path(
            self.config,  'input',
            '{}NamelistFileName'.format(self.componentName))
        self.namelist = NameList(namelistFileName)

        streamsFileName = build_config_full_path(
            self.config, 'input',
            '{}StreamsFileName'.format(self.componentName))
        self.runStreams = StreamsFile(streamsFileName,
                                      streamsdir=self.runDirectory)
        self.historyStreams = StreamsFile(streamsFileName,
                                          streamsdir=self.historyDirectory)

        self.calendar = self.namelist.get('config_calendar_type')

        make_directories(self.plotsDirectory)

        if self.componentName == 'ocean':
            self.namelistMap = oceanNamelistMap
            self.streamMap = oceanStreamMap
            self.variableMap = oceanVariableMap
        elif self.componentName == 'seaIce':
            self.namelistMap = seaIceNamelistMap
            self.streamMap = seaIceStreamMap
            self.variableMap = seaIceVariableMap
        else:
            self.namelistMap = None
            self.streamMap = None
            self.variableMap = None

        # set the start and end dates for each type of analysis
        for tag in ['climatology', 'timeSeries', 'index']:
            if tag in self.tags:
                self.set_start_end_date(section=tag)

        # }}}

    def run(self):  # {{{
        '''
        Runs the analysis task.

        Individual tasks (children classes of this base class) should first
        call this method to perform any common steps in an analysis task,
        then, perform the steps required to run the analysis task.

        Authors
        -------
        Xylar Asay-Davis
        '''
        return  # }}}

    def check_generate(self):
        # {{{
        '''
        Determines if this analysis should be generated, based on the
        `generate` config option and `taskName`, `componentName` and
        `tags`.

        Individual tasks do not need to create their own versions of this
        function.

        Returns
        -------
        generate : bool
            Whether or not this task should be run.

        Raises
        ------
        ValueError : If one of `self.taskName`, `self.componentName`
            or `self.tags` has not been set.

        Authors
        -------
        Xylar Asay-Davis
        '''

        for memberName in ['taskName', 'componentName', 'tags']:
            if not hasattr(self, memberName):
                raise ValueError('Analysis tasks must define self.{} in their '
                                 '__init__ method.'.format(memberName))

        if (not isinstance(self.tags, list) and
                self.tags is not None):
            raise ValueError('Analysis tasks\'s member self.tags '
                             'must be NOne or a list of strings.')

        config = self.config
        generateList = config.getExpression('output', 'generate')
        generate = False
        for element in generateList:
            if '_' in element:
                (prefix, suffix) = element.split('_', 1)
            else:
                prefix = element
                suffix = None

            allSuffixes = [self.componentName]
            if self.tags is not None:
                allSuffixes = allSuffixes + self.tags
            noSuffixes = [self.taskName] + allSuffixes
            if prefix == 'all':
                if (suffix in allSuffixes) or (suffix is None):
                    generate = True
            elif prefix == 'no':
                if suffix in noSuffixes:
                    generate = False
            elif element == self.taskName:
                generate = True

        return generate  # }}}

    def check_analysis_enabled(self, analysisOptionName, default=False,
                               raiseException=True):
        '''
        Check to make sure a given analysis is turned on, issuing a warning or
        raising an exception if not.

        Parameters
        ----------
        analysisOptionName : str
            The name of a boolean namelist option indicating whether the given
            analysis member is enabled

        default : bool, optional
            If no analysis option with the given name can be found, indicates
            whether the given analysis is assumed to be enabled by default.

        raiseException : bool, optional
            Whether

        Returns
        -------
        enabled : bool
            Whether the given analysis is enabled

        Raises
        ------
        RuntimeError
            If the given analysis option is not found and ``default`` is not
            ``True`` or if the analysis option is found and is ``False``.  The
            exception is only raised if ``raiseException = True``.

        Authors
        -------
        Xylar Asay-Davis
        '''

        try:
            if self.namelistMap is None:
                optionName = analysisOptionName
            else:
                optionName = self.namelist.find_option(
                    self.namelistMap[analysisOptionName])
            enabled = self.namelist.getbool(optionName)
        except ValueError:
            enabled = default
            if default:
                message = 'WARNING: namelist option {} not found.\n' \
                          'This likely indicates that the simulation you ' \
                          'are analyzing was run with an\n' \
                          'older version of MPAS-O that did not support ' \
                          'this flag.  Assuming enabled.'.format(
                              analysisOptionName)
                warnings.warn(message)

        if not enabled and raiseException:
            raise RuntimeError('*** MPAS-Analysis relies on {} = .true.\n'
                               '*** Make sure to enable this analysis '
                               'member.'.format(analysisOptionName))

        return enabled

    def set_start_end_date(self, section):  # {{{
        '''
        Set the start and end dates in the ``config`` correspond to the start
        and end years in a given category of analysis

        Parameters
        ----------
        section : str
            The name of a section in the config file containing ``startYear``
            and ``endYear`` options. ``section`` is typically one of
            ``climatology``, ``timeSeries`` or ``index``

        Authors
        -------
        Xylar Asay-Davis
        '''

        if not self.config.has_option(section, 'startDate'):
            startDate = '{:04d}-01-01_00:00:00'.format(
                self.config.getint(section, 'startYear'))
            self.config.set(section, 'startDate', startDate)
        if not self.config.has_option(section, 'endDate'):
            endDate = '{:04d}-12-31_23:59:59'.format(
                self.config.getint(section, 'endYear'))
            self.config.set(section, 'endDate', endDate)  # }}}

    def update_start_end_date(self, section, streamName):  # {{{
        '''
        Update the start and end dates (and years) based on the times found
        in the given stream.  Cache the times if they are not already cached.

        Parameters
        ----------
        section : str
            The name of a section in the config file containing ``startYear``
            and ``endYear`` options. ``section`` is typically one of
            ``climatology``, ``timeSeries`` or ``index``

        streamName : str
            The name of a stream from which to read (and cache) the times

        Returns
        -------
        changed : bool
            Whether the start and end dates were updated.

        Authors
        -------
        Xylar Asay-Davis
        '''

        startDate = self.config.get(section, 'startDate')
        endDate = self.config.get(section, 'endDate')
        startDate = string_to_days_since_date(dateString=startDate,
                                              calendar=self.calendar)
        endDate = string_to_days_since_date(dateString=endDate,
                                            calendar=self.calendar)

        inFileNames = self.get_input_file_names(
            streamName, startAndEndDateSection=section)

        fullTimeCache = self.cache_multifile_dataset_times(
            inFileNames, streamName, timeVariableName='Time')

        # find only those cached times between starDate and endDate
        times = []
        for fileName in fullTimeCache:
            localTimes = fullTimeCache[fileName]['times']
            mask = numpy.logical_and(localTimes >= startDate,
                                     localTimes < endDate)
            if numpy.count_nonzero(mask) == 0:
                continue

            times.extend(list(localTimes[mask]))

        requestedStartYear = self.config.getint('climatology', 'startYear')
        requestedEndYear = self.config.getint('climatology', 'endYear')

        startYear = days_to_datetime(numpy.amin(times),
                                     calendar=self.calendar).year
        endYear = days_to_datetime(numpy.amax(times),
                                   calendar=self.calendar).year
        changed = False
        if startYear != requestedStartYear or endYear != requestedEndYear:
            message = "{} start and/or end year different from " \
                      "requested\n" \
                      "requested: {:04d}-{:04d}\n" \
                      "actual:    {:04d}-{:04d}\n".format(section,
                                                          requestedStartYear,
                                                          requestedEndYear,
                                                          startYear,
                                                          endYear)
            warnings.warn(message)
            self.config.set(section, 'startYear', str(startYear))
            self.config.set(section, 'endYear', str(endYear))

            startDate = '{:04d}-01-01_00:00:00'.format(startYear)
            self.config.set(section, 'startDate', startDate)
            endDate = '{:04d}-12-31_23:59:59'.format(endYear)
            self.config.set(section, 'endDate', endDate)  # }}}

            changed = True

        return changed  # }}}

    def get_input_file_names(self, streamName,
                             startDate=None, endDate=None,
                             startAndEndDateSection=None):  # {{{
        '''
        Get a list of input files corresponding to the given stream and
        optionally bounded by the start and end dates found in the given
        section of the config file.

        Parameters
        ----------
        streamName : str
            The name of a stream to check.  If ``self.streamMap`` is defined,
            the streamName will be mapped to the corresponding name in the
            streams file

        startDate, endDate : float, optional
            start and end date to use in determining which files to include in
            the list

        startAndEndDateSection : str, optional
            If ``startDate`` and ``endDate`` arguments are not supplied, the
            name of a section in the config file containing ``startDate`` and
            ``endDate`` options to use instead. ``startAndEndDateSection`` is
            typically one of ``climatology``, ``timeSeries`` or ``index``.

        Raises
        ------
        RuntimeError
            If no files are found in the desired date range.

        Authors
        -------
        Xylar Asay-Davis
        '''

        if startDate is None and endDate is None and \
                startAndEndDateSection is not None:
            startDate = self.config.get(startAndEndDateSection, 'startDate')
            endDate = self.config.get(startAndEndDateSection, 'endDate')

        if self.streamMap is not None:
            streamName = \
                self.historyStreams.find_stream(self.streamMap[streamName])
        inputFileNames = self.historyStreams.readpath(streamName,
                                                      startDate=startDate,
                                                      endDate=endDate,
                                                      calendar=self.calendar)

        if len(inputFileNames) == 0:
            raise RuntimeError('No input files could be found in stream {} '
                               'between {} and {}'.format(streamName,
                                                          startDate, endDate))
        return inputFileNames  # }}}

    def add_prerequisite_tasks(self, prerequisiteTasks):  # {{{
        '''
        Add one or more task names to the list of prerequisite tasks

        Parameters
        ----------
        prerequisiteTasks : list of str
            The names of the tasks to add as prerequisites

        Authors
        -------
        Xylar Asay-Davis
        '''

        if self.prerequisiteTasks is None:
            self.prerequisiteTasks = prerequisiteTasks
        else:
            for task in prerequisiteTasks:
                if task not in self.prerequisiteTasks:
                    self.prerequisiteTasks.apend(task)

        # }}}
# }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
