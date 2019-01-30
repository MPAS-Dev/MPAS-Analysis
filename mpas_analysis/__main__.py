#!/usr/bin/env python
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2018 Los Alamos National Security, LLC. All rights reserved.
# Copyright (c) 2018 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2018 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

"""
Runs MPAS-Analysis via a configuration file (e.g. `config.analysis`)
specifying analysis options.
"""
# Authors
# -------
# Xylar Asay-Davis, Phillip J. Wolfram, Milena Veneziani

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import mpas_analysis

import argparse
import traceback
import sys
import pkg_resources
import shutil
import os
from collections import OrderedDict
import progressbar
import logging
import xarray

from mpas_analysis.shared.analysis_task import AnalysisFormatter

from mpas_analysis.configuration import MpasAnalysisConfigParser

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories

from mpas_analysis.shared.html import generate_html

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.analysis_task import \
    update_time_bounds_from_file_names

from mpas_analysis.shared.plot.plotting import _register_custom_colormaps, \
    _plot_color_gradients

from mpas_analysis import ocean
from mpas_analysis import sea_ice
from mpas_analysis.shared.climatology import MpasClimatologyTask, \
    RefYearMpasClimatologyTask
from mpas_analysis.shared.time_series import MpasTimeSeriesTask

from mpas_analysis.shared.io.download import download_files


def update_time_bounds_in_config(config):  # {{{
    """
    Updates the start and end year (and associated full date) for
    climatologies, time series and climate indices based on the files that are
    actually available.

    Parameters
    ----------
    config : ``MpasAnalysisConfigParser`` object
        contains config options

    """
    # By updating the bounds for each component, we should end up with the
    # more constrained time bounds if any component has less output than others
    for componentName in ['ocean', 'seaIce']:
        for section in ['climatology', 'timeSeries', 'index']:
            update_time_bounds_from_file_names(config, section, componentName)
    # }}}


def build_analysis_list(config, controlConfig):  # {{{
    """
    Build a list of analysis tasks. New tasks should be added here, following
    the approach used for existing analysis tasks.

    Parameters
    ----------
    config : ``MpasAnalysisConfigParser`` object
        contains config options

    controlConfig : ``MpasAnalysisConfigParser`` object
        contains config options for a control run, or ``None`` if no config
        file for a control run was specified

    Returns
    -------
    analyses : list of ``AnalysisTask`` objects
        A list of all analysis tasks
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    analyses = []

    # Ocean Analyses
    oceanClimatolgyTask = MpasClimatologyTask(config=config,
                                              componentName='ocean')
    oceanTimeSeriesTask = MpasTimeSeriesTask(config=config,
                                             componentName='ocean')
    oceanIndexTask = MpasTimeSeriesTask(config=config,
                                        componentName='ocean',
                                        section='index')

    oceanRefYearClimatolgyTask = RefYearMpasClimatologyTask(
            config=config,  componentName='ocean')

    analyses.append(oceanClimatolgyTask)
    analyses.append(oceanRefYearClimatolgyTask)

    analyses.append(ocean.ClimatologyMapMLD(config, oceanClimatolgyTask,
                                            controlConfig))
    analyses.append(ocean.ClimatologyMapSST(config, oceanClimatolgyTask,
                                            controlConfig))
    analyses.append(ocean.ClimatologyMapSSS(config, oceanClimatolgyTask,
                                            controlConfig))
    analyses.append(ocean.ClimatologyMapSSH(config, oceanClimatolgyTask,
                                            controlConfig))
    analyses.append(ocean.ClimatologyMapEKE(config, oceanClimatolgyTask,
                                            controlConfig))
    analyses.append(ocean.ClimatologyMapOHCAnomaly(
            config, oceanClimatolgyTask, oceanRefYearClimatolgyTask,
            controlConfig))

    analyses.append(ocean.ClimatologyMapSose(
            config, oceanClimatolgyTask, controlConfig))
    analyses.append(ocean.ClimatologyMapBGC(config, oceanClimatolgyTask,
                                            controlConfig))

    analyses.append(ocean.ClimatologyMapArgoTemperature(
            config, oceanClimatolgyTask, controlConfig))
    analyses.append(ocean.ClimatologyMapArgoSalinity(
            config, oceanClimatolgyTask, controlConfig))

    analyses.append(ocean.ClimatologyMapSchmidtko(
            config, oceanClimatolgyTask, controlConfig))

    analyses.append(ocean.ClimatologyMapAntarcticMelt(config,
                                                      oceanClimatolgyTask,
                                                      controlConfig))

    analyses.append(ocean.TimeSeriesAntarcticMelt(config, oceanTimeSeriesTask,
                                                  controlConfig))

    analyses.append(ocean.TimeSeriesTemperatureAnomaly(config,
                                                       oceanTimeSeriesTask))
    analyses.append(ocean.TimeSeriesSalinityAnomaly(config,
                                                    oceanTimeSeriesTask))
    analyses.append(ocean.TimeSeriesOHCAnomaly(config,
                                               oceanTimeSeriesTask,
                                               controlConfig))
    analyses.append(ocean.TimeSeriesSST(config, oceanTimeSeriesTask,
                                        controlConfig))
    analyses.append(ocean.MeridionalHeatTransport(config, oceanClimatolgyTask,
                                                  controlConfig))

    analyses.append(ocean.StreamfunctionMOC(config, oceanClimatolgyTask,
                                            controlConfig))
    analyses.append(ocean.IndexNino34(config, oceanIndexTask, controlConfig))

    analyses.append(ocean.WoceTransects(config, oceanClimatolgyTask,
                                        controlConfig))

    analyses.append(ocean.SoseTransects(config, oceanClimatolgyTask,
                                        controlConfig))

    analyses.append(ocean.GeojsonTransects(config, oceanClimatolgyTask,
                                           controlConfig))

    analyses.append(ocean.OceanRegionalProfiles(config, controlConfig))

    # Sea Ice Analyses
    seaIceClimatolgyTask = MpasClimatologyTask(config=config,
                                               componentName='seaIce')
    seaIceTimeSeriesTask = MpasTimeSeriesTask(config=config,
                                              componentName='seaIce')

    analyses.append(seaIceClimatolgyTask)
    analyses.append(sea_ice.ClimatologyMapSeaIceConc(
            config=config, mpasClimatologyTask=seaIceClimatolgyTask,
            hemisphere='NH', controlConfig=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceThick(
            config=config, mpasClimatologyTask=seaIceClimatolgyTask,
            hemisphere='NH', controlConfig=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceConc(
            config=config, mpasClimatologyTask=seaIceClimatolgyTask,
            hemisphere='SH', controlConfig=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceThick(
            config=config, mpasClimatologyTask=seaIceClimatolgyTask,
            hemisphere='SH', controlConfig=controlConfig))
    analyses.append(seaIceTimeSeriesTask)

    analyses.append(sea_ice.TimeSeriesSeaIce(config, seaIceTimeSeriesTask,
                                             controlConfig))

    # Iceberg Analyses
    analyses.append(sea_ice.ClimatologyMapIcebergConc(
            config=config, mpasClimatologyTask=seaIceClimatolgyTask,
            hemisphere='SH', controlConfig=controlConfig))

    return analyses  # }}}


def determine_analyses_to_generate(analyses, verbose):  # {{{
    """
    Build a list of analysis tasks to run based on the 'generate' config
    option (or command-line flag) and prerequisites and subtasks of each
    requested task.  Each task's ``setup_and_check`` method is called in the
    process.

    Parameters
    ----------
    analyses : list of ``AnalysisTask`` objects
        A list of all analysis tasks

    verbose : bool
        Whether to write out a full stack trace when exceptions occur during
        ``setup_and_check()`` calls for each task

    Returns
    -------
    analysesToGenerate : ``OrderedDict`` of ``AnalysisTask`` objects
        A dictionary of analysis tasks to run
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    analysesToGenerate = OrderedDict()
    # check which analysis we actually want to generate and only keep those
    for analysisTask in analyses:
        # update the dictionary with this task and perhaps its subtasks
        add_task_and_subtasks(analysisTask, analysesToGenerate, verbose)

    return analysesToGenerate  # }}}


def add_task_and_subtasks(analysisTask, analysesToGenerate, verbose,
                          callCheckGenerate=True):
    # {{{
    """
    If a task has been requested through the generate config option or
    if it is a prerequisite of a requested task, add it to the dictionary of
    tasks to generate.

    Parameters
    ----------
    analysisTask : ``AnalysisTask``
        A task to be added

    analysesToGenerate : ``OrderedDict`` of ``AnalysisTask``
        The list of analysis tasks to be generated, which this call may
        update to include this task and its subtasks

    verbose : bool
        Whether to write out a full stack trace when exceptions occur during
        ``setup_and_check()`` calls for each task

    callCheckGenerate : bool
        Whether the ``check_generate`` method should be call for this task to
        see if it has been requested.  We skip this for subtasks and
        prerequisites, since they are needed by another task regardless of
        whether the user specifically requested them.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    key = (analysisTask.taskName, analysisTask.subtaskName)
    if key in analysesToGenerate.keys():
        # The task was already added
        if analysisTask._setupStatus != 'success':
            ValueError("task {} already added but this version was not set up "
                       "successfully. Typically, this indicates two tasks "
                       "with the same full name".format(
                               analysisTask.fullTaskName))
        return

    # for each anlaysis task, check if we want to generate this task
    # and if the analysis task has a valid configuration
    taskTitle = analysisTask.printTaskName
    if callCheckGenerate and not analysisTask.check_generate():
        # we don't need to add this task -- it wasn't requested
        return

    # first, we should try to add the prerequisites of this task and its
    # subtasks (if they aren't also subtasks for this task)
    prereqs = analysisTask.runAfterTasks
    for subtask in analysisTask.subtasks:
        for prereq in subtask.runAfterTasks:
            if prereq not in analysisTask.subtasks:
                prereqs.extend(subtask.runAfterTasks)

    for prereq in prereqs:
        add_task_and_subtasks(prereq, analysesToGenerate, verbose,
                              callCheckGenerate=False)
        if prereq._setupStatus != 'success':
            # a prereq failed setup_and_check
            print("Warning: prerequisite of {} failed during check, "
                  "so this task will not be run".format(
                      taskTitle))
            analysisTask._setupStatus = 'fail'
            return

    # make sure all prereqs have been set up successfully before trying to
    # set up this task -- this task's setup may depend on setup in the prereqs
    try:
        analysisTask.setup_and_check()
    except (Exception, BaseException):
        if verbose:
            traceback.print_exc(file=sys.stdout)
        print("Warning: {} failed during check and will not be run".format(
            taskTitle))
        analysisTask._setupStatus = 'fail'
        return

    # next, we should try to add the subtasks.  This is done after the current
    # analysis task has been set up in case subtasks depend on information
    # from the parent task
    for subtask in analysisTask.subtasks:
        add_task_and_subtasks(subtask, analysesToGenerate, verbose,
                              callCheckGenerate=False)
        if subtask._setupStatus != 'success':
            # a subtask failed setup_and_check
            print("Warning: subtask of {} failed during check, "
                  "so this task will not be run".format(
                      taskTitle))
            analysisTask._setupStatus = 'fail'
            return

    analysesToGenerate[key] = analysisTask
    analysisTask._setupStatus = 'success'
    # }}}


def update_generate(config, generate):  # {{{
    """
    Update the 'generate' config option using a string from the command line.

    Parameters
    ----------
    config : ``MpasAnalysisConfigParser`` object
        contains config options

    generate : str
        a comma-separated string of generate flags: either names of analysis
        tasks or commands of the form ``all_<tag>`` or ``no_<tag>`` indicating
        that analysis with a given tag should be included or excluded).
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    # overwrite the 'generate' in config with a string that parses to
    # a list of string
    generateList = generate.split(',')
    generateString = ', '.join(["'{}'".format(element)
                                for element in generateList])
    generateString = '[{}]'.format(generateString)
    config.set('output', 'generate', generateString)  # }}}


def run_analysis(config, analyses):  # {{{
    """
    Run all the tasks, either in serial or in parallel

    Parameters
    ----------
    config : ``MpasAnalysisConfigParser`` object
        contains config options

    analyses : OrderedDict of ``AnalysisTask`` objects
        A dictionary of analysis tasks to run with (task, subtask) names as
        keys
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    # write the config file the log directory
    logsDirectory = build_config_full_path(config, 'output',
                                           'logsSubdirectory')

    mainRunName = config.get('runs', 'mainRunName')

    configFileName = '{}/config.{}'.format(logsDirectory, mainRunName)

    configFile = open(configFileName, 'w')
    config.write(configFile)
    configFile.close()

    parallelTaskCount = config.getWithDefault('execute', 'parallelTaskCount',
                                              default=1)

    isParallel = parallelTaskCount > 1 and len(analyses) > 1

    for analysisTask in analyses.values():
        if not analysisTask.runAfterTasks and not analysisTask.subtasks:
            analysisTask._runStatus.value = AnalysisTask.READY
        else:
            analysisTask._runStatus.value = AnalysisTask.BLOCKED

    tasksWithErrors = []
    runningTasks = {}

    # redirect output to a log file
    logsDirectory = build_config_full_path(config, 'output',
                                           'logsSubdirectory')

    logFileName = '{}/taskProgress.log'.format(logsDirectory)

    logger = logging.getLogger('mpas_analysis')
    handler = logging.FileHandler(logFileName)

    formatter = AnalysisFormatter()
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    totalTaskCount = len(analyses)
    widgets = ['Running tasks: ', progressbar.Percentage(), ' ',
               progressbar.Bar(), ' ', progressbar.ETA()]
    progress = progressbar.ProgressBar(widgets=widgets,
                                       maxval=totalTaskCount).start()

    # run each analysis task
    while True:
        # we still have tasks to run
        for analysisTask in analyses.values():
            if analysisTask._runStatus.value == AnalysisTask.BLOCKED:
                prereqs = analysisTask.runAfterTasks + analysisTask.subtasks
                prereqStatus = [prereq._runStatus.value for prereq in prereqs]
                if any([runStatus == AnalysisTask.FAIL for runStatus in
                        prereqStatus]):
                    # a prerequisite failed so this task cannot succeed
                    analysisTask._runStatus.value = AnalysisTask.FAIL
                if all([runStatus == AnalysisTask.SUCCESS for runStatus in
                        prereqStatus]):
                    # no unfinished prerequisites so we can run this task
                    analysisTask._runStatus.value = AnalysisTask.READY

        unfinishedCount = 0
        for analysisTask in analyses.values():
            if analysisTask._runStatus.value not in [AnalysisTask.SUCCESS,
                                                     AnalysisTask.FAIL]:
                unfinishedCount += 1

        progress.update(totalTaskCount-unfinishedCount)

        if unfinishedCount <= 0 and len(runningTasks.keys()) == 0:
            # we're done
            break

        # launch new tasks
        for key, analysisTask in analyses.items():
            if analysisTask._runStatus.value == AnalysisTask.READY:
                if isParallel:
                    logger.info('Running {}'.format(
                            analysisTask.printTaskName))
                    analysisTask._runStatus.value = AnalysisTask.RUNNING
                    analysisTask.start()
                    runningTasks[key] = analysisTask
                    if len(runningTasks.keys()) >= parallelTaskCount:
                        break
                else:
                    analysisTask.run(writeLogFile=False)

        if isParallel:
            # wait for a task to finish
            analysisTask = wait_for_task(runningTasks)
            key = (analysisTask.taskName, analysisTask.subtaskName)
            runningTasks.pop(key)

            taskTitle = analysisTask.printTaskName

            if analysisTask._runStatus.value == AnalysisTask.SUCCESS:
                logger.info("   Task {} has finished successfully.".format(
                        taskTitle))
            elif analysisTask._runStatus.value == AnalysisTask.FAIL:
                message = "ERROR in task {}.  See log file {} for " \
                          "details".format(taskTitle,
                                           analysisTask._logFileName)
                logger.error(message)
                print(message)
                tasksWithErrors.append(taskTitle)
            else:
                message = "Unexpected status from in task {}.  This may be " \
                          "a bug.".format(taskTitle)
                logger.error(message)
                print(message)
        else:
            if analysisTask._runStatus.value == AnalysisTask.FAIL:
                sys.exit(1)

    progress.finish()

    # blank line to make sure remaining output is on a new line
    print('')

    logger.handlers = []

    # raise the last exception so the process exits with an error
    errorCount = len(tasksWithErrors)
    if errorCount == 1:
        print("There were errors in task {}".format(tasksWithErrors[0]))
        sys.exit(1)
    elif errorCount > 0:
        print("There were errors in {} tasks: {}".format(
            errorCount, ', '.join(tasksWithErrors)))
        print("See log files in {} for details.".format(logsDirectory))
        print("The following commands may be helpful:")
        print("  cd {}".format(logsDirectory))
        print("  grep Error *.log")
        sys.exit(1)
    else:
        print('Log files for executed tasks can be found in {}'.format(
                logsDirectory))

    # }}}


def wait_for_task(runningTasks, timeout=0.1):  # {{{
    """
    Build a list of analysis modules based on the 'generate' config option.
    New tasks should be added here, following the approach used for existing
    analysis tasks.

    Parameters
    ----------
    runningTasks : dict of ``AnalysisTasks``
        The tasks that are currently running, with task names as keys

    Returns
    -------
    analysisTask : ``AnalysisTasks``
        A task that finished
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    # necessary to have a timeout so we can kill the whole thing
    # with a keyboard interrupt
    while True:
        for analysisTask in runningTasks.values():
            analysisTask.join(timeout=timeout)
            if not analysisTask.is_alive():
                return analysisTask  # }}}


def purge_output(config):
    outputDirectory = config.get('output', 'baseDirectory')
    if not os.path.exists(outputDirectory):
        print('Output directory {} does not exist.\n'
              'No purge necessary.'.format(outputDirectory))
    else:
        for subdirectory in ['plots', 'logs', 'mpasClimatology', 'mapping',
                             'timeSeries', 'html', 'mask', 'profiles']:
            option = '{}Subdirectory'.format(subdirectory)
            directory = build_config_full_path(
                    config=config, section='output',
                    relativePathOption=option)
            if os.path.exists(directory):
                print('Deleting contents of {}'.format(directory))
                if os.path.islink(directory):
                    os.unlink(directory)
                else:
                    shutil.rmtree(directory)

        for component in ['ocean', 'seaIce']:
            for subdirectory in ['climatology', 'remappedClim']:
                option = '{}Subdirectory'.format(subdirectory)
                section = '{}Observations'.format(component)
                directory = build_config_full_path(
                        config=config, section='output',
                        relativePathOption=option,
                        relativePathSection=section)
                if os.path.exists(directory):
                    print('Deleting contents of {}'.format(directory))
                    if os.path.islink(directory):
                        os.unlink(directory)
                    else:
                        shutil.rmtree(directory)


def symlink_main_run(config, defaultConfig):
    '''
    Create symlinks to the climatology and time-series directories for the
    main run that has aleady been computed so we don't have to recompute
    the analysis.
    '''

    def link_dir(section, option):
        destDirectory = build_config_full_path(config=config, section='output',
                                               relativePathOption=option,
                                               relativePathSection=section)
        if not os.path.exists(destDirectory):

            destBase, _ = os.path.split(destDirectory)

            make_directories(destBase)

            sourceDirectory = build_config_full_path(
                    config=mainConfig, section='output',
                    relativePathOption=option, relativePathSection=section)

            os.symlink(sourceDirectory, destDirectory)

    mainConfigFile = config.get('runs', 'mainRunConfigFile')
    if not os.path.exists(mainConfigFile):
        raise OSError('A main config file {} was specified but the '
                      'file does not exist'.format(mainConfigFile))
    mainConfigFiles = [mainConfigFile]
    if defaultConfig is not None:
        mainConfigFiles = [defaultConfig] + mainConfigFiles
    mainConfig = MpasAnalysisConfigParser()
    mainConfig.read(mainConfigFiles)

    for subdirectory in ['mpasClimatology', 'timeSeries', 'mapping', 'mask',
                         'profiles']:
        section = 'output'
        option = '{}Subdirectory'.format(subdirectory)
        link_dir(section=section, option=option)

    for component in ['ocean', 'seaIce']:
        for subdirectory in ['climatology', 'remappedClim']:
            section = '{}Observations'.format(component)
            option = '{}Subdirectory'.format(subdirectory)
            link_dir(section=section, option=option)


def main():

    """
    Entry point for the main script ``mpas_analysis``
    """

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-v', '--version',
                        action='version',
                        version='mpas_analysis {}'.format(
                                mpas_analysis.__version__),
                        help="Show version number and exit")
    parser.add_argument("--setup_only", dest="setup_only", action='store_true',
                        help="If only the setup phase, not the run or HTML "
                        "generation phases, should be executed.")
    parser.add_argument("--html_only", dest="html_only", action='store_true',
                        help="If only the setup and HTML generation phases, "
                        "not the run phase, should be executed.")
    parser.add_argument("-g", "--generate", dest="generate",
                        help="A list of analysis modules to generate "
                        "(nearly identical generate option in config file).",
                        metavar="ANALYSIS1[,ANALYSIS2,ANALYSIS3,...]")
    parser.add_argument("-l", "--list", dest="list", action='store_true',
                        help="List the available analysis tasks")
    parser.add_argument("-p", "--purge", dest="purge", action='store_true',
                        help="Purge the analysis by deleting the output"
                        "directory before running")
    parser.add_argument('configFiles', metavar='CONFIG',
                        type=str, nargs='*', help='config file')
    parser.add_argument("--plot_colormaps", dest="plot_colormaps",
                        action='store_true',
                        help="Make a plot displaying all available colormaps")
    parser.add_argument("--verbose", dest="verbose", action='store_true',
                        help="Verbose error reporting during setup-and-check "
                             "phase")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    for configFile in args.configFiles:
        if not os.path.exists(configFile):
            raise OSError('Config file {} not found.'.format(configFile))

    # add config.default to cover default not included in the config files
    # provided on the command line
    if pkg_resources.resource_exists('mpas_analysis', 'config.default'):
        defaultConfig = pkg_resources.resource_filename('mpas_analysis',
                                                        'config.default')
        configFiles = [defaultConfig] + args.configFiles
    else:
        print('WARNING: Did not find config.default.  Assuming other config '
              'file(s) contain a\n'
              'full set of configuration options.')
        defaultConfig = None
        configFiles = args.configFiles

    config = MpasAnalysisConfigParser()
    config.read(configFiles)

    if args.list:
        analyses = build_analysis_list(config, controlConfig=None)
        for analysisTask in analyses:
            print('task: {}'.format(analysisTask.taskName))
            print('    component: {}'.format(analysisTask.componentName)),
            print('    tags: {}'.format(', '.join(analysisTask.tags)))
        sys.exit(0)

    if args.plot_colormaps:
        _register_custom_colormaps()
        _plot_color_gradients()
        sys.exit(0)

    if config.has_option('runs', 'controlRunConfigFile'):
        controlConfigFile = config.get('runs', 'controlRunConfigFile')
        if not os.path.exists(controlConfigFile):
            raise OSError('A control config file {} was specified but the '
                          'file does not exist'.format(controlConfigFile))
        controlConfigFiles = [controlConfigFile]
        if defaultConfig is not None:
            controlConfigFiles = [defaultConfig] + controlConfigFiles
        controlConfig = MpasAnalysisConfigParser()
        controlConfig.read(controlConfigFiles)

        # replace the log directory so log files get written to this run's
        # log directory, not the control run's
        logsDirectory = build_config_full_path(config, 'output',
                                               'logsSubdirectory')

        controlConfig.set('output', 'logsSubdirectory', logsDirectory)

        print('Comparing to control run {} rather than observations. \n'
              'Make sure that MPAS-Analysis has been run previously with the '
              'control config file.'.format(controlConfig.get('runs',
                                                              'mainRunName')))
    else:
        controlConfig = None

    if args.purge:
        purge_output(config)

    if config.has_option('runs', 'mainRunConfigFile'):
        symlink_main_run(config, defaultConfig)

    if args.generate:
        update_generate(config, args.generate)

    if controlConfig is not None:
        # we want to use the "generate" option from the current run, not
        # the control config file
        controlConfig.set('output', 'generate', config.get('output',
                                                           'generate'))

    logsDirectory = build_config_full_path(config, 'output',
                                           'logsSubdirectory')
    make_directories(logsDirectory)

    update_time_bounds_in_config(config)

    file_cache_maxsize = config.getint('input', 'file_cache_maxsize')
    try:
        xarray.set_options(file_cache_maxsize=file_cache_maxsize)
    except ValueError:
        # xarray version doesn't support file_cache_maxsize yet...
        pass

    analyses = build_analysis_list(config, controlConfig)
    analyses = determine_analyses_to_generate(analyses, args.verbose)

    if not args.setup_only and not args.html_only:
        run_analysis(config, analyses)

    if not args.setup_only:
        generate_html(config, analyses, controlConfig)


def download_analysis_data():
    """
    Entry point for downloading the input data set from public repository for
    MPAS-Analysis to work. The input data set includes: pre-processed
    observations data, MPAS mapping files and MPAS regional mask files
    (which are used for the MOC computation), for a subset of MPAS meshes.
    """

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-o", "--outDir", dest="outDir", required=True,
                        help="Directory where MPAS-Analysis input data will"
                             "be downloaded")
    args = parser.parse_args()

    try:
        os.makedirs(args.outDir)
    except OSError:
        pass

    urlBase = 'https://web.lcrc.anl.gov/public/e3sm/diagnostics'
    analysisFileList = pkg_resources.resource_string(
            'mpas_analysis', 'obs/analysis_input_files').decode('utf-8')

    # remove any empty strings from the list
    analysisFileList = list(filter(None, analysisFileList.split('\n')))
    download_files(analysisFileList, urlBase, args.outDir)


if __name__ == "__main__":
    main()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
