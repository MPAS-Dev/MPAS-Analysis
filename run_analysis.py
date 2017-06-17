#!/usr/bin/env python

"""
Runs MPAS-Analysis via a configuration file (e.g. `config.analysis`)
specifying analysis options.

Authors
-------
Xylar Asay-Davis, Phillip J. Wolfram
"""

import os
import matplotlib as mpl
import argparse
import traceback
import sys
import warnings
import subprocess
import time
from collections import OrderedDict

from mpas_analysis.configuration.MpasAnalysisConfigParser \
    import MpasAnalysisConfigParser

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories


def build_analysis_list(config, isSubtask):  # {{{
    """
    Build a list of analysis tasks based on the 'generate' config option.

    Authors
    -------
    Xylar Asay-Davis
    """

    # choose the right rendering backend, depending on whether we're displaying
    # to the screen
    if not config.getboolean('plot', 'displayToScreen'):
        mpl.use('Agg')

    # analysis can only be imported after the right MPL renderer is selected
    from mpas_analysis import ocean
    from mpas_analysis import sea_ice
    from mpas_analysis.shared.cache_dataset_times_task \
        import CacheDatasetTimesTask

    # analyses will be a list of analysis classes
    analyses = []

    # add cacheOceanTimeSeriesStatsTimes task for caching the times in
    # MPAS-Ocean timeSeriesStatsMonthly output files.  This is a prerequisite
    # for all ocean analysis.
    analyses.append(CacheDatasetTimesTask(
            config=config,
            componentName='ocean',
            streamName='timeSeriesStats',
            startAndEndDateSections=['climatology', 'timeSeries', 'index'],
            namelistOption='config_am_timeseriesstatsmonthly_enable'))

    # add cacheSeaIceTimeSeriesStatsTimes task for caching the times in
    # MPAS-Ocean timeSeriesStatsMonthly output files.  This is a prerequisite
    # for all sea ice analysis.
    analyses.append(CacheDatasetTimesTask(
            config=config,
            componentName='seaIce',
            streamName='timeSeriesStats',
            startAndEndDateSections=['climatology', 'timeSeries'],
            namelistOption='config_am_timeseriesstatsmonthly_enable'))

    # Ocean Analyses
    analyses.append(ocean.TimeSeriesOHC(config))
    analyses.append(ocean.TimeSeriesSST(config))
    analyses.append(ocean.IndexNino34(config))
    analyses.append(ocean.MeridionalHeatTransport(config))
    analyses.append(ocean.StreamfunctionMOC(config))

    analyses.append(ocean.ClimatologyMapSST(config))
    analyses.append(ocean.ClimatologyMapMLD(config))
    analyses.append(ocean.ClimatologyMapSSS(config))

    # Sea Ice Analyses
    analyses.append(sea_ice.TimeSeriesSeaIce(config))
    analyses.append(sea_ice.ClimatologyMapSeaIce(config))

    possibleAnalyses = OrderedDict()
    for analysisTask in analyses:
        possibleAnalyses[analysisTask.taskName] = analysisTask

    # check which analysis we actually want to generate and only keep those
    analysesToGenerate = OrderedDict()
    for analysisTask in possibleAnalyses.itervalues():
        # update the dictionary with this task and perhaps its prerequisites
        analysesToAdd = add_task_and_prereqisites(analysisTask,
                                                  possibleAnalyses,
                                                  analysesToGenerate,
                                                  isPrerequisite=False,
                                                  isSubtask=isSubtask)
        analysesToGenerate.update(analysesToAdd)

    return analysesToGenerate  # }}}


def add_task_and_prereqisites(analysisTask, possibleAnalyses,
                              analysesToGenerate, isPrerequisite,
                              isSubtask):  # {{{
    """
    If a task has been requested through the generate config option or
    if it is a prerequisite of a requested task, add it to the dictionary of
    tasks to generate.

    Authors
    -------
    Xylar Asay-Davis
    """

    analysesToAdd = OrderedDict()
    # for each anlaysis task, check if we want to generate this task
    # and if the analysis task has a valid configuration
    if isPrerequisite or analysisTask.check_generate():
        add = False
        try:
            analysisTask.setup_and_check()
            add = True
        except:
            traceback.print_exc(file=sys.stdout)
            print "ERROR: analysis task {} failed during check and " \
                "will not be run".format(analysisTask.taskName)
        if add and not isSubtask:
            # first, we should try to add the prerequisites
            prereqs = analysisTask.prerequisiteTasks
            if prereqs is not None:
                for prereq in prereqs:
                    if prereq not in analysesToGenerate.keys():
                        prereqToAdd = add_task_and_prereqisites(
                                possibleAnalyses[prereq], possibleAnalyses,
                                analysesToGenerate, isPrerequisite=True)
                        if len(prereqToAdd.keys()) == 0:
                            # a prerequisite failed setup_and_check
                            print "ERROR: a prerequisite of analysis task {}" \
                                  " failed during check and will not be" \
                                  " run".format(analysisTask.taskName)
                            add = False
                            break
                        # the prerequisite (and its prerequisites) should be
                        # added
                        analysesToAdd.update(prereqToAdd)
        if add:
            analysesToAdd[analysisTask.taskName] = analysisTask

    return analysesToAdd  # }}}


def update_generate(config, generate):  # {{{
    """
    Update the 'generate' config option using a string from the command line.

    Authors
    -------
    Xylar Asay-Davis
    """

    # overwrite the 'generate' in config with a string that parses to
    # a list of string
    generateList = generate.split(',')
    generateString = ', '.join(["'{}'".format(element)
                                for element in generateList])
    generateString = '[{}]'.format(generateString)
    config.set('output', 'generate', generateString)  # }}}


def run_analysis(config, analyses, configFiles, isSubtask):  # {{{
    """
    Run this script once each for several parallel tasks.

    Author
    ------
    Xylar Asay-Davis
    """

    taskCount = config.getWithDefault('execute', 'parallelTaskCount',
                                      default=1)

    isParallel = not isSubtask and taskCount > 1 and len(analyses) > 1

    for analysisTask in analyses.itervalues():
        if analysisTask.prerequisiteTasks is None or isSubtask:
            analysisTask.status = 'ready'
        else:
            analysisTask.status = 'blocked'

    processes = {}
    logs = {}

    # run each analysis task
    lastException = None

    runningCount = 0
    while True:
        # we still have tasks to run
        for analysisTask in analyses.itervalues():
            if analysisTask.status == 'blocked':
                prereqStatus = [analyses[prereq].status for prereq in
                                analysisTask.prerequisiteTasks]
                if any([status == 'fail' for status in prereqStatus]):
                    # a prerequisite failed so this task cannot succeed
                    analysisTask.status = 'fail'
                if all([status == 'success' for status in prereqStatus]):
                    # no unfinished prerequisites so we can run this task
                    analysisTask.status = 'ready'

        unfinishedCount = 0
        for analysisTask in analyses.itervalues():
            if analysisTask.status not in ['success', 'fail']:
                unfinishedCount += 1

        if unfinishedCount <= 0:
            # we're done
            break

        # launch new tasks
        for taskName, analysisTask in analyses.items():
            if analysisTask.status == 'ready':
                if isParallel:
                    process, logFile = launch_task(taskName, config,
                                                   configFiles)
                    processes[taskName] = process
                    logs[taskName] = logFile
                    analysisTask.status = 'running'
                    runningCount += 1
                    if runningCount >= taskCount:
                        break
                else:
                    exception = run_task(config, analysisTask)
                    if exception is None:
                        analysisTask.status = 'success'
                    else:
                        lastException = exception
                        analysisTask.status = 'fail'

        if isParallel:
            # wait for a task to finish
            (taskName, process) = wait_for_task(processes)
            analysisTask = analyses[taskName]
            runningCount -= 1
            processes.pop(taskName)
            if process.returncode == 0:
                print "Task {} has finished successfully.".format(taskName)
                analysisTask.status = 'success'
            else:
                print "ERROR in task {}.  See log file {} for details".format(
                    taskName, logs[taskName].name)
                analysisTask.status = 'fail'
            logs[taskName].close()

    if not isParallel and config.getboolean('plot', 'displayToScreen'):
        import matplotlib.pyplot as plt
        plt.show()

    # raise the last exception so the process exits with an error
    if lastException is not None:
        raise lastException

    # }}}


def launch_task(taskName, config, configFiles):  # {{{
    """
    Launch a parallel tasks

    Authors
    -------
    Xylar Asay-Davis
    """
    thisFile = os.path.realpath(__file__)

    commandPrefix = config.getWithDefault('execute', 'commandPrefix',
                                          default='')
    if commandPrefix == '':
        commandPrefix = []
    else:
        commandPrefix = commandPrefix.split(' ')

    args = commandPrefix + [thisFile, '--subtask', '--generate', taskName] \
        + configFiles

    logFileName = '{}/{}.log'.format(logsDirectory, taskName)

    # write the command to the log file
    logFile = open(logFileName, 'w')
    logFile.write('Command: {}\n'.format(' '.join(args)))
    # make sure the command gets written before the rest of the log
    logFile.flush()
    print 'Running {}'.format(taskName)
    process = subprocess.Popen(args, stdout=logFile,
                               stderr=subprocess.STDOUT)

    return (process, logFile)  # }}}


def wait_for_task(processes):  # {{{
    """
    Wait for the next process to finish and check its status.  Returns both the
    task name and the process that finished.

    Authors
    -------
    Xylar Asay-Davis
    """

    # first, check if any process has already finished
    for taskName, process in processes.iteritems():  # python 2.7!
        if(not is_running(process)):
            return (taskName, process)

    # No process has already finished, so wait for the next one
    (pid, status) = os.waitpid(-1, 0)
    for taskName, process in processes.iteritems():
        if pid == process.pid:
            process.returncode = status
            # since we used waitpid, this won't happen automatically
            return (taskName, process)  # }}}


def is_running(process):  # {{{
    """
    Returns whether a given process is currently running

    Authors
    -------
    Xylar Asay-Davis
    """

    try:
        os.kill(process.pid, 0)
    except OSError:
        return False
    else:
        return True  # }}}


def run_task(config, analysisTask):  # {{{
    """
    Run a single analysis task, time the task, write out the config file
    (including any modifications specific to the task) and return the exception
    raised (if any)

    Authors
    -------
    Xylar Asay-Davis
    """

    # write out a copy of the configuration to document the run
    logsDirectory = build_config_full_path(config, 'output',
                                           'logsSubdirectory')
    exception = None
    try:
        startTime = time.clock()
        analysisTask.run()
        runDuration = time.clock() - startTime
        m, s = divmod(runDuration, 60)
        h, m = divmod(int(m), 60)
        print 'Execution time: {}:{:02d}:{:05.2f}'.format(h, m, s)
    except (Exception, BaseException) as e:
        if isinstance(e, KeyboardInterrupt):
            raise e
        traceback.print_exc(file=sys.stdout)
        print "ERROR: analysis task {} failed during run".format(
            analysisTask.taskName)
        exception = e

    configFileName = '{}/configs/config.{}'.format(logsDirectory,
                                                   analysisTask.taskName)
    configFile = open(configFileName, 'w')
    config.write(configFile)
    configFile.close()

    return exception  # }}}


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--subtask", dest="subtask", action='store_true',
                        help="If this is a subtask when running parallel "
                             "tasks")
    parser.add_argument("-g", "--generate", dest="generate",
                        help="A list of analysis modules to generate "
                             "(nearly identical generate option in config "
                             "file).",
                        metavar="ANALYSIS1[,ANALYSIS2,ANALYSIS3,...]")
    parser.add_argument('configFiles', metavar='CONFIG',
                        type=str, nargs='+', help='config file')
    args = parser.parse_args()

    # add config.default to cover default not included in the config files
    # provided on the command line
    defaultConfig = '{}/config.default'.format(
        os.path.dirname(os.path.realpath(__file__)))
    if os.path.exists(defaultConfig):
        configFiles = [defaultConfig] + args.configFiles
    else:
        warnings.warn('WARNING: Did not find config.default.  Assuming other '
                      'config file(s) contain a\n'
                      'full set of configuration options.')
        configFiles = args.configFiles

    config = MpasAnalysisConfigParser()
    config.read(configFiles)

    if args.generate:
        update_generate(config, args.generate)

    logsDirectory = build_config_full_path(config, 'output',
                                           'logsSubdirectory')
    make_directories(logsDirectory)
    make_directories('{}/configs/'.format(logsDirectory))

    analyses = build_analysis_list(config, args.subtask)

    run_analysis(config, analyses, configFiles, args.subtask)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
