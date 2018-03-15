# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#
'''
Defines the base class for analysis tasks.

Authors
-------
Xylar Asay-Davis

'''

from __future__ import absolute_import, division, print_function, \
    unicode_literals

from multiprocessing import Process, Value
import time
import traceback
import logging
import sys

from mpas_analysis.shared.io import NameList, StreamsFile
from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories


class AnalysisTask(Process):  # {{{
    '''
    The base class for analysis tasks.

    Attributes
    ----------
    config : ``MpasAnalysisConfigParser``
        Contains configuration options

    taskName : str
        The name of the task, typically the same as the class name except
        starting with lowercase (e.g. 'myTask' for class 'MyTask')

    componentName : {'ocean', 'seaIce'}
        The name of the component (same as the folder where the task
        resides)

    tags : list of str
        Tags used to describe the task (e.g. 'timeSeries', 'climatology',
        horizontalMap', 'index', 'transect').  These are used to determine
        which tasks are generated (e.g. 'all_transect' or 'no_climatology'
        in the 'generate' flags)

    runDirectory : str
        The base input directory for namelists, streams files and restart files

    historyDirectory : str
        The base input directory for history files

    plotsDirectory : str
        The directory for writing plots (which is also created if it doesn't
        exist)

    namelist : ``shared.io.NameList``
        the namelist reader

    runStreams : ``shared.io.StreamsFile``
        the streams file reader for streams in the run directory (e.g. restart
        files)

    historyStreams : ``shared.io.StreamsFile``
        the streams file reader for streams in the history directory (most
        streams other than restart files)

    calendar : {'gregorian', 'gregoraian_noleap'}
        The calendar used in the MPAS run

    runAfterTasks : list of ``AnalysisTasks``
        tasks that must be complete before this task can run

    subtasks : ``OrderedDict`` of ``AnalysisTasks``
        Subtasks of this task, with subtask names as keys

    xmlFileNames : list of strings
        The XML file associated with each plot produced by this analysis, empty
        if no plots were produced

    logger : ``logging.Logger``
        A logger for output during the run phase of an analysis task

    Authors
    -------
    Xylar Asay-Davis

    '''

    # flags for run status
    UNSET = 0
    READY = 1
    BLOCKED = 2
    RUNNING = 3
    SUCCESS = 4
    FAIL = 5

    def __init__(self, config, taskName, componentName, tags=[],
                 subtaskName=None):  # {{{
        '''
        Construct the analysis task.

        Individual tasks (children classes of this base class) should first
        call this method to perform basic initialization, then, define the
        ``taskName``, ``componentName`` and list of ``tags`` for the task.

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

        subtaskName : str, optional
            If this is a subtask of ``taskName``, the name of the subtask

        Authors
        -------
        Xylar Asay-Davis
        '''
        if subtaskName is None:
            self.fullTaskName = taskName
            self.printTaskName = taskName
        else:
            self.fullTaskName = '{}_{}'.format(taskName, subtaskName)
            self.printTaskName = '{}: {}'.format(taskName, subtaskName)

        # call the constructor from the base class (Process)
        super(AnalysisTask, self).__init__(name=self.fullTaskName)

        self.config = config
        self.taskName = taskName
        self.subtaskName = subtaskName
        self.componentName = componentName
        self.tags = tags
        self.subtasks = []
        self.logger = None
        self.runAfterTasks = []
        self.xmlFileNames = []

        # non-public attributes related to multiprocessing and logging
        self.daemon = True
        self._setupStatus = None
        self._runStatus = Value('i', AnalysisTask.UNSET)
        self._stackTrace = None
        self._logFileName = None
        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis (e.g. reading namelists and
        streams files).

        After this call, the following attributes are set (see documentation
        for the class):
        runDirectory, historyDirectory, plotsDirectory, namelist, runStreams,
        historyStreams, calendar

        Individual tasks (children classes of this base class) should first
        call this method to perform basic setup, then, check whether the
        configuration is correct for a given analysis and perform additional,
        analysis-specific setup.  For example, this function could check if
        necessary observations and other data files are found, then, determine
        the list of files to be read when the analysis is run.

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

        # set the start and end dates for each type of analysis
        for tag in ['climatology', 'timeSeries', 'index']:
            if tag in self.tags:
                self.set_start_end_date(section=tag)

        # redirect output to a log file
        logsDirectory = build_config_full_path(self.config, 'output',
                                               'logsSubdirectory')

        self._logFileName = '{}/{}.log'.format(logsDirectory,
                                               self.fullTaskName)

        # }}}

    def run_task(self):  # {{{
        '''
        Run the analysis.  Each task should override this function to do the
        work of computing and/or plotting analysis

        Authors
        -------
        Xylar Asay-Davis
        '''
        return  # }}}

    def run_after(self, task):  # {{{
        '''
        Only run this task after the given task has completed.  This allows a
        task to be constructed of multiple subtasks, some of which may block
        later tasks, while allowing some subtasks to run in parallel.  It also
        allows for tasks to depend on other tasks (e.g. for computing
        climatologies or extracting time series for many variables at once).

        Parameters
        ----------
        task : ``AnalysisTask``
            The task that should finish before this one begins

        Authors
        -------
        Xylar Asay-Davis
        '''

        self.runAfterTasks.append(task)
        # }}}

    def add_subtask(self, subtask):  # {{{
        '''
        Add a subtask to this tasks.  This task always runs after the subtask
        has finished.  However, this task gets set up *before* the subtask,
        so the setup of the subtask can depend on fields defined during the
        setup of this task (the parent).

        Parameters
        ----------
        subtask : ``AnalysisTask``
            The subtask to run as part of this task

        Authors
        -------
        Xylar Asay-Davis
        '''

        if subtask not in self.subtasks:
            self.subtasks.append(subtask)
        # }}}

    def run(self, writeLogFile=True):  # {{{
        '''
        Sets up logging and then runs the analysis task.

        Parameters
        ----------
        writeLogFile : bool, optional
            If ``True``, output to stderr and stdout get written to a log file.
            Otherwise, the internal logger ``self.logger`` points to stdout
            and no log file is created.  The intention is for logging to take
            place in parallel mode but not in serial mode.

        Authors
        -------
        Xylar Asay-Davis
        '''
        # redirect output to a log file
        if writeLogFile:
            self.logger = logging.getLogger(self.fullTaskName)
            handler = logging.FileHandler(self._logFileName)
        else:
            self.logger = logging.getLogger()
            handler = logging.StreamHandler(sys.stdout)

        formatter = AnalysisFormatter()
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        self.logger.setLevel(logging.INFO)

        if writeLogFile:
            oldStdout = sys.stdout
            oldStderr = sys.stderr
            sys.stdout = StreamToLogger(self.logger, logging.INFO)
            sys.stderr = StreamToLogger(self.logger, logging.ERROR)

        startTime = time.time()
        try:
            self.run_task()
            self._runStatus.value = AnalysisTask.SUCCESS
        except (Exception, BaseException) as e:
            if isinstance(e, KeyboardInterrupt):
                raise e
            self._stackTrace = traceback.format_exc()
            self.logger.error("analysis task {} failed during run \n"
                              "{}".format(self.fullTaskName, self._stackTrace))
            self._runStatus.value = AnalysisTask.FAIL

        runDuration = time.time() - startTime
        m, s = divmod(runDuration, 60)
        h, m = divmod(int(m), 60)
        self.logger.info('Execution time: {}:{:02d}:{:05.2f}'.format(h, m, s))

        if writeLogFile:
            # restore stdout and stderr
            sys.stdout = oldStdout
            sys.stderr = oldStderr

        # remove the handlers from the logger (probably only necessary if
        # writeLogFile==False)
        self.logger.handlers = []

        # }}}

    def check_generate(self):
        # {{{
        '''
        Determines if this analysis should be generated, based on the
        ``generate`` config option and ``taskName``, ``componentName`` and
        ``tags``.

        Individual tasks do not need to create their own versions of this
        function.

        Returns
        -------
        generate : bool
            Whether or not this task should be run.

        Raises
        ------
        ValueError : If one of ``self.taskName``, ``self.componentName``
            or ``self.tags`` has not been set.

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
                             'must be None or a list of strings.')

        config = self.config
        generateList = config.getExpression('output', 'generate')
        if len(generateList) > 0 and generateList[0][0:5] == 'only_':
            # add 'all' if the first item in the list has the 'only' prefix.
            # Otherwise, we would not run any tasks
            generateList = ['all'] + generateList
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
            if prefix == 'only':
                if suffix not in allSuffixes:
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
            optionName = analysisOptionName
            enabled = self.namelist.getbool(optionName)
        except ValueError:
            enabled = default
            if default:
                print('Warning: namelist option {} not found.\n'
                      'This likely indicates that the simulation you '
                      'are analyzing was run with an\n'
                      'older version of MPAS-O that did not support '
                      'this flag.  Assuming enabled.'.format(
                              analysisOptionName))

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

# }}}


class AnalysisFormatter(logging.Formatter):  # {{{
    """
    A custom formatter for logging

    Modified from:
    https://stackoverflow.com/a/8349076/7728169

    Authors
    -------
    Xylar Asay-Davis
    """

    # printing error messages without a prefix because they are sometimes
    # errors and sometimes only warnings sent to stderr
    dbg_fmt = "DEBUG: %(module)s: %(lineno)d: %(msg)s"
    info_fmt = "%(msg)s"
    err_fmt = info_fmt

    def __init__(self, fmt=info_fmt):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._fmt = AnalysisFormatter.dbg_fmt

        elif record.levelno == logging.INFO:
            self._fmt = AnalysisFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._fmt = AnalysisFormatter.err_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._fmt = format_orig

        return result
# }}}


class StreamToLogger(object):  # {{{
    """
    Modified based on code by:
    https://www.electricmonk.nl/log/2011/08/14/redirect-stdout-and-stderr-to-a-logger-in-python/

    Copyright (C) 2011 Ferry Boender

    License: "available under the GPL" (the author does not provide more
    details)

    Fake file-like stream object that redirects writes to a logger instance.
    """
    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, str(line.rstrip()))

    def flush(self):
        pass

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
