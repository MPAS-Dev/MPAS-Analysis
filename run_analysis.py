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

from mpas_analysis.configuration import MpasAnalysisConfigParser

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories

from mpas_analysis.shared.html import generate_html

def update_generate(config, generate):  # {{{
    """
    Update the 'generate' config option using a string from the command line.

    Author: Xylar Asay-Davis
    """

    # overwrite the 'generate' in config with a string that parses to
    # a list of string
    generateList = generate.split(',')
    generateString = ', '.join(["'{}'".format(element)
                                for element in generateList])
    generateString = '[{}]'.format(generateString)
    config.set('output', 'generate', generateString)  # }}}


def run_parallel_tasks(config, analyses, configFiles, taskCount):
    # {{{
    """
    Run this script once each for several parallel tasks.

    Author
    ------
    Xylar Asay-Davis
    """

    taskNames = [analysisTask.taskName for analysisTask in analyses]

    taskCount = min(taskCount, len(taskNames))

    (processes, logs) = launch_tasks(taskNames[0:taskCount], config,
                                     configFiles)
    remainingTasks = taskNames[taskCount:]
    while len(processes) > 0:
        (taskName, process) = wait_for_task(processes)
        if process.returncode == 0:
            print "Task {} has finished successfully.".format(taskName)
        else:
            print "ERROR in task {}.  See log file {} for details".format(
                taskName, logs[taskName].name)
        logs[taskName].close()
        # remove the process from the process dictionary (no need to bother)
        processes.pop(taskName)

        if len(remainingTasks) > 0:
            (process, log) = launch_tasks(remainingTasks[0:1], config,
                                          configFiles)
            # merge the new process and log into these dictionaries
            processes.update(process)
            logs.update(log)
            remainingTasks = remainingTasks[1:]
    # }}}


def launch_tasks(taskNames, config, configFiles):  # {{{
    """
    Launch one or more tasks

    Author: Xylar Asay-Davis
    """
    thisFile = os.path.realpath(__file__)

    commandPrefix = config.getWithDefault('execute', 'commandPrefix',
                                          default='')
    if commandPrefix == '':
        commandPrefix = []
    else:
        commandPrefix = commandPrefix.split(' ')

    processes = {}
    logs = {}
    for taskName in taskNames:
        args = commandPrefix + \
            [thisFile, '--subtask', '--generate', taskName] + configFiles

        logFileName = '{}/{}.log'.format(logsDirectory, taskName)

        # write the command to the log file
        logFile = open(logFileName, 'w')
        logFile.write('Command: {}\n'.format(' '.join(args)))
        # make sure the command gets written before the rest of the log
        logFile.flush()
        print 'Running {}'.format(taskName)
        process = subprocess.Popen(args, stdout=logFile,
                                   stderr=subprocess.STDOUT)
        processes[taskName] = process
        logs[taskName] = logFile

    return (processes, logs)  # }}}


def wait_for_task(processes):  # {{{
    """
    Wait for the next process to finish and check its status.  Returns both the
    task name and the process that finished.

    Author: Xylar Asay-Davis
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

    Author: Xylar Asay-Davis
    """

    try:
        os.kill(process.pid, 0)
    except OSError:
        return False
    else:
        return True  # }}}


def build_analysis_list(config):  # {{{
    """
    Build a list of analysis modules based on the 'generate' config option.

    Author: Xylar Asay-Davis
    """

    # choose the right rendering backend, depending on whether we're displaying
    # to the screen
    if not config.getboolean('plot', 'displayToScreen'):
        mpl.use('Agg')

    # analysis can only be imported after the right MPL renderer is selected
    from mpas_analysis import ocean
    from mpas_analysis import sea_ice

    # analyses will be a list of analysis classes
    analyses = []

    # Ocean Analyses
    analyses.append(ocean.TimeSeriesOHC(config))
    analyses.append(ocean.TimeSeriesSST(config))
    analyses.append(ocean.IndexNino34(config))
    analyses.append(ocean.MeridionalHeatTransport(config))
    analyses.append(ocean.StreamfunctionMOC(config))

    analyses.append(ocean.ClimatologyMapMLD(config))
    analyses.append(ocean.ClimatologyMapSST(config))
    analyses.append(ocean.ClimatologyMapSSS(config))

    # Sea Ice Analyses
    analyses.append(sea_ice.TimeSeriesSeaIce(config))
    analyses.append(sea_ice.ClimatologyMapSeaIceConc(config, hemisphere='NH'))
    analyses.append(sea_ice.ClimatologyMapSeaIceConc(config, hemisphere='SH'))
    analyses.append(sea_ice.ClimatologyMapSeaIceThick(config, hemisphere='NH'))
    analyses.append(sea_ice.ClimatologyMapSeaIceThick(config, hemisphere='SH'))

    # check which analysis we actually want to generate and only keep those
    analysesToGenerate = []
    for analysisTask in analyses:
        # for each anlaysis module, check if we want to generate this task
        # and if the analysis task has a valid configuration
        if analysisTask.check_generate():
            add = False
            try:
                analysisTask.setup_and_check()
                add = True
            except:
                traceback.print_exc(file=sys.stdout)
                print "ERROR: analysis module {} failed during check and " \
                    "will not be run".format(analysisTask.taskName)
            if add:
                analysesToGenerate.append(analysisTask)

    return analysesToGenerate  # }}}


def run_analysis(config, analyses):  # {{{

    # run each analysis task
    lastException = None
    for analysisTask in analyses:
        # write out a copy of the configuration to document the run
        logsDirectory = build_config_full_path(config, 'output',
                                               'logsSubdirectory')
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
            print "ERROR: analysis module {} failed during run".format(
                analysisTask.taskName)
            lastException = e

        configFileName = '{}/configs/config.{}'.format(logsDirectory,
                                                       analysisTask.taskName)
        configFile = open(configFileName, 'w')
        config.write(configFile)
        configFile.close()

    if config.getboolean('plot', 'displayToScreen'):
        import matplotlib.pyplot as plt
        plt.show()

    # raise the last exception so the process exits with an error
    if lastException is not None:
        raise lastException

    return  # }}}


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--subtask", dest="subtask", action='store_true',
                        help="If this is a subtask when running parallel "
                             "tasks")
    parser.add_argument("-g", "--generate", dest="generate",
                        help="A list of analysis modules to generate "
                        "(nearly identical generate option in config file).",
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

    analyses = build_analysis_list(config)

    parallelTaskCount = config.getWithDefault('execute', 'parallelTaskCount',
                                              default=1)

    if parallelTaskCount <= 1 or len(analyses) == 1:
        run_analysis(config, analyses)
    else:
        run_parallel_tasks(config, analyses, configFiles, parallelTaskCount)

    if not args.subtask:
        generate_html(config, analyses)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
