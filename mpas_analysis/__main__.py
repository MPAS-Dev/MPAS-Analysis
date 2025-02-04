#!/usr/bin/env python
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

"""
Runs MPAS-Analysis via a configuration file (e.g. `analysis.cfg`)
specifying analysis options.
"""
# Authors
# -------
# Xylar Asay-Davis, Phillip J. Wolfram, Milena Veneziani

import mpas_analysis
import mpas_analysis.version

import argparse
import traceback
import sys
import shutil
import os
from collections import OrderedDict
import progressbar
import logging
import xarray
import time
import json
from importlib.metadata import Distribution
from importlib.resources import contents

from mache import discover_machine, MachineInfo

from mpas_tools.config import MpasConfigParser

from mpas_analysis.shared.analysis_task import AnalysisFormatter

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories, copyfile

from mpas_analysis.shared.html import generate_html

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.analysis_task import \
    update_time_bounds_from_file_names

from mpas_analysis.shared.plot.colormap import register_custom_colormaps, \
    _plot_color_gradients

from mpas_analysis import ocean
from mpas_analysis import sea_ice
from mpas_analysis.shared.climatology import MpasClimatologyTask, \
    RefYearMpasClimatologyTask
from mpas_analysis.shared.time_series import MpasTimeSeriesTask

from mpas_analysis.shared.regions import ComputeRegionMasks


def update_time_bounds_in_config(config):
    """
    Updates the start and end year (and associated full date) for
    climatologies, time series and climate indices based on the files that are
    actually available.

    Parameters
    ----------
    config : mpas_tools.config.MpasConfigParser
        contains config options

    """
    # By updating the bounds for each component, we should end up with the
    # more constrained time bounds if any component has less output than others
    for componentName in ['ocean', 'seaIce']:
        for section in ['climatology', 'timeSeries', 'index']:
            update_time_bounds_from_file_names(config, section, componentName)


def build_analysis_list(config, controlConfig):
    """
    Build a list of analysis tasks. New tasks should be added here, following
    the approach used for existing analysis tasks.

    Parameters
    ----------
    config : mpas_tools.config.MpasConfigParser
        contains config options

    controlConfig : mpas_tools.config.MpasConfigParser or None
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
    oceanClimatologyTasks = {}
    for op in ['avg', 'min', 'max']:
        oceanClimatologyTasks[op] = MpasClimatologyTask(config=config,
                                                        componentName='ocean',
                                                        op=op)
    oceanTimeSeriesTask = MpasTimeSeriesTask(config=config,
                                             componentName='ocean')
    oceanIndexTask = MpasTimeSeriesTask(config=config,
                                        componentName='ocean',
                                        section='index')

    oceanRefYearClimatologyTask = RefYearMpasClimatologyTask(
        config=config, componentName='ocean')

    oceanRegionMasksTask = ComputeRegionMasks(config=config,
                                              conponentName='ocean')

    for op in oceanClimatologyTasks:
        analyses.append(oceanClimatologyTasks[op])
    analyses.append(oceanRefYearClimatologyTask)

    analyses.append(ocean.ClimatologyMapMLD(config,
                                            oceanClimatologyTasks['avg'],
                                            controlConfig))

    analyses.append(ocean.ClimatologyMapMLDMinMax(config,
                                                  oceanClimatologyTasks,
                                                  controlConfig))

    analyses.append(ocean.ClimatologyMapSST(config,
                                            oceanClimatologyTasks['avg'],
                                            controlConfig))
    analyses.append(ocean.ClimatologyMapFluxes(config,
                                               oceanClimatologyTasks['avg'],
                                               controlConfig,
                                               fluxType='mass'))
    analyses.append(ocean.ClimatologyMapFluxes(config,
                                               oceanClimatologyTasks['avg'],
                                               controlConfig,
                                               fluxType='heat'))
    analyses.append(ocean.ClimatologyMapSSS(config,
                                            oceanClimatologyTasks['avg'],
                                            controlConfig))
    analyses.append(ocean.ClimatologyMapSSH(config,
                                            oceanClimatologyTasks['avg'],
                                            controlConfig))
    analyses.append(ocean.ClimatologyMapEKE(config,
                                            oceanClimatologyTasks['avg'],
                                            controlConfig))
    analyses.append(ocean.ClimatologyMapVel(config,
                                            oceanClimatologyTasks['avg'],
                                            controlConfig))
    analyses.append(ocean.ClimatologyMapBSF(config,
                                            oceanClimatologyTasks['avg'],
                                            controlConfig))
    analyses.append(ocean.ClimatologyMapOHCAnomaly(
        config, oceanClimatologyTasks['avg'], oceanRefYearClimatologyTask,
        controlConfig))

    analyses.append(ocean.ClimatologyMapSose(
        config, oceanClimatologyTasks['avg'], controlConfig))
    analyses.append(ocean.ClimatologyMapWoa(
        config, oceanClimatologyTasks['avg'], controlConfig))
    analyses.append(ocean.ClimatologyMapBGC(config,
                                            oceanClimatologyTasks['avg'],
                                            controlConfig))

    analyses.append(ocean.ClimatologyMapArgoTemperature(
        config, oceanClimatologyTasks['avg'], controlConfig))
    analyses.append(ocean.ClimatologyMapArgoSalinity(
        config, oceanClimatologyTasks['avg'], controlConfig))

    analyses.append(ocean.ClimatologyMapSchmidtko(
        config, oceanClimatologyTasks['avg'], controlConfig))

    analyses.append(ocean.ClimatologyMapAntarcticMelt(
        config, oceanClimatologyTasks['avg'], oceanRegionMasksTask,
        controlConfig))

    analyses.append(ocean.ConservationTask(
        config, controlConfig))

    analyses.append(ocean.RegionalTSDiagrams(
        config, oceanClimatologyTasks['avg'], oceanRegionMasksTask,
        controlConfig))

    analyses.append(ocean.TimeSeriesAntarcticMelt(config, oceanTimeSeriesTask,
                                                  oceanRegionMasksTask,
                                                  controlConfig))

    analyses.append(ocean.TimeSeriesOceanRegions(config, oceanRegionMasksTask,
                                                 controlConfig))

    analyses.append(ocean.TimeSeriesTemperatureAnomaly(config,
                                                       oceanTimeSeriesTask))
    analyses.append(ocean.TimeSeriesSalinityAnomaly(config,
                                                    oceanTimeSeriesTask))
    analyses.append(ocean.TimeSeriesOHCAnomaly(config,
                                               oceanTimeSeriesTask,
                                               controlConfig))
    analyses.append(ocean.TimeSeriesSSHAnomaly(config,
                                               oceanTimeSeriesTask,
                                               controlConfig))
    analyses.append(ocean.TimeSeriesSST(config, oceanTimeSeriesTask,
                                        controlConfig))
    analyses.append(ocean.TimeSeriesTransport(config, controlConfig))

    analyses.append(ocean.OceanHistogram(config, oceanClimatologyTasks['avg'],
                                         oceanRegionMasksTask,
                                         controlConfig))
    analyses.append(ocean.MeridionalHeatTransport(
        config, oceanClimatologyTasks['avg'], controlConfig))

    analyses.append(ocean.StreamfunctionMOC(config,
                                            oceanClimatologyTasks['avg'],
                                            controlConfig))
    analyses.append(ocean.IndexNino34(config, oceanIndexTask, controlConfig))

    analyses.append(ocean.WoceTransects(config, oceanClimatologyTasks['avg'],
                                        controlConfig))

    analyses.append(ocean.AntshipTransects(config,
                                           oceanClimatologyTasks['avg'],
                                           controlConfig))

    analyses.append(ocean.OsnapTransects(config, oceanClimatologyTasks['avg'],
                                         controlConfig))

    analyses.append(ocean.SoseTransects(config, oceanClimatologyTasks['avg'],
                                        controlConfig))

    analyses.append(ocean.WoaTransects(config, oceanClimatologyTasks['avg'],
                                       controlConfig))

    analyses.append(ocean.GeojsonTransects(config,
                                           oceanClimatologyTasks['avg'],
                                           controlConfig))

    oceanRegionalProfiles = ocean.OceanRegionalProfiles(
        config, oceanRegionMasksTask, controlConfig)
    analyses.append(oceanRegionalProfiles)

    analyses.append(ocean.HovmollerOceanRegions(
        config, oceanRegionMasksTask, oceanRegionalProfiles, controlConfig))

    # Sea Ice Analyses
    seaIceClimatologyTask = MpasClimatologyTask(config=config,
                                                componentName='seaIce')
    seaIceTimeSeriesTask = MpasTimeSeriesTask(config=config,
                                              componentName='seaIce')

    analyses.append(seaIceClimatologyTask)
    analyses.append(sea_ice.ClimatologyMapSeaIceConc(
        config=config, mpasClimatologyTask=seaIceClimatologyTask,
        hemisphere='NH', controlConfig=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceThick(
        config=config, mpasClimatologyTask=seaIceClimatologyTask,
        hemisphere='NH', controlConfig=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceSnowDepth(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='NH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceAreaFractionRidge(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='NH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceVolumeRidge(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='NH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceAlbedo(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='NH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceProduction(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='NH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceMelting(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='NH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceAreaTendencyThermo(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='NH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceAreaTendencyTransp(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='NH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceVolumeTendencyThermo(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='NH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceVolumeTendencyTransp(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='NH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceConc(
        config=config, mpasClimatologyTask=seaIceClimatologyTask,
        hemisphere='SH', controlConfig=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceThick(
        config=config, mpasClimatologyTask=seaIceClimatologyTask,
        hemisphere='SH', controlConfig=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceSnowDepth(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='SH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceAreaFractionRidge(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='SH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceVolumeRidge(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='SH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceAlbedo(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='SH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceProduction(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='SH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapRiskIndexOutcome(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='NH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapRiskIndexOutcome(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='SH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceMelting(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='SH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceAreaTendencyThermo(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='SH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceAreaTendencyTransp(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='SH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceVolumeTendencyThermo(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='SH', control_config=controlConfig))
    analyses.append(sea_ice.ClimatologyMapSeaIceVolumeTendencyTransp(
        config=config, mpas_climatology_task=seaIceClimatologyTask,
        hemisphere='SH', control_config=controlConfig))

    analyses.append(seaIceTimeSeriesTask)
    analyses.append(sea_ice.TimeSeriesSeaIce(config, seaIceTimeSeriesTask,
                                             controlConfig))

    # Iceberg Analyses
    analyses.append(sea_ice.ClimatologyMapIcebergConc(
        config=config, mpasClimatologyTask=seaIceClimatologyTask,
        hemisphere='SH', controlConfig=controlConfig))

    # Wave Analyses
    analyses.append(ocean.ClimatologyMapWaves(
        config, oceanClimatologyTasks['avg'], oceanRegionMasksTask,
        controlConfig))

    check_for_duplicate_names(analyses)

    return analyses


def check_for_duplicate_names(analyses):
    """
    Check for duplicate taskName and subtaskName in the list of analysis tasks
    and their subtasks

    Parameters
    ----------
    analyses : list of mpas_analysis.shared.AnalysisTask
        A list of all analysis tasks
    """
    all_task_names = []
    errors = []
    for analysis in analyses:
        mainTaskName = analysis.taskName
        assert(analysis.subtaskName is None)
        fullName = (mainTaskName, None)
        if fullName in all_task_names:
            errors.append(
                f'A task named {mainTaskName} has been added more than once')
        all_task_names.append(fullName)
        for subtask in analysis.subtasks:
            taskName = subtask.taskName
            subtaskName = subtask.subtaskName
            if taskName != mainTaskName:
                errors.append(
                    f'A subtask named {taskName}: {subtaskName} has a '
                    f'different task name than its parent task: \n'
                    f'  {mainTaskName}')
                fullName = (taskName, subtaskName)
                if fullName in all_task_names:
                    errors.append(
                        f'A subtask named {taskName}: {subtaskName} has been '
                        f'added more than once')
                all_task_names.append(fullName)

    if len(errors) > 0:
        all_errors = '\n  '.join(errors)
        raise ValueError(f'Analysis tasks failed these checks:\n'
                         f'  {all_errors}')


def determine_analyses_to_generate(analyses, verbose):
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

    totalFailures = 0

    print('')

    analysesToGenerate = OrderedDict()
    # check which analysis we actually want to generate and only keep those
    for analysisTask in analyses:
        # update the dictionary with this task and perhaps its subtasks
        failureCount = add_task_and_subtasks(analysisTask, analysesToGenerate,
                                             verbose)

        totalFailures += failureCount

    if totalFailures > 0:
        print('\n{} tasks and subtasks failed during setup.'.format(
            totalFailures))
        if not verbose:
            print('To find out why these tasks are failing, use the --verbose '
                  'flag')

    print('')

    return analysesToGenerate


def add_task_and_subtasks(analysisTask, analysesToGenerate, verbose,
                          callCheckGenerate=True):

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

    totalFailures = 0

    key = (analysisTask.taskName, analysisTask.subtaskName)
    if key in analysesToGenerate.keys():
        # The task was already added
        if analysisTask._setupStatus != 'success':
            ValueError("task {} already added but this version was not set up "
                       "successfully. Typically, this indicates two tasks "
                       "with the same full name".format(
                           analysisTask.fullTaskName))
        return totalFailures

    # for each analysis task, check if we want to generate this task
    # and if the analysis task has a valid configuration
    taskTitle = analysisTask.printTaskName
    if callCheckGenerate and not analysisTask.check_generate():
        # we don't need to add this task -- it wasn't requested
        return totalFailures

    # first, we should try to add the prerequisites of this task and its
    # subtasks (if they aren't also subtasks for this task)
    prereqs = analysisTask.runAfterTasks
    for subtask in analysisTask.subtasks:
        for prereq in subtask.runAfterTasks:
            if prereq not in analysisTask.subtasks:
                prereqs.extend(subtask.runAfterTasks)

    for prereq in prereqs:
        failureCount = add_task_and_subtasks(prereq, analysesToGenerate,
                                             verbose, callCheckGenerate=False)
        totalFailures += failureCount
        if prereq._setupStatus != 'success':
            if failureCount == 0:
                 raise ValueError(f'Error: prerequisite {prereq.printTaskName} of '
                                  f'{taskTitle} did not set up successfully but also '
                                  'did not indicate a failure.  This likely indicates '
                                  'a bug like multiple tasks with the same name.')
            # a prereq failed setup_and_check
            print(f'Warning: prerequisite of {taskTitle} failed during check, '
                  'so this task will not be run')
            analysisTask._setupStatus = 'fail'
            totalFailures += 1
            return totalFailures

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
        totalFailures += 1
        return totalFailures

    # next, we should try to add the subtasks.  This is done after the current
    # analysis task has been set up in case subtasks depend on information
    # from the parent task
    for subtask in analysisTask.subtasks:
        failureCount = add_task_and_subtasks(subtask, analysesToGenerate,
                                             verbose, callCheckGenerate=False)
        totalFailures += failureCount
        if subtask._setupStatus != 'success':
            if failureCount == 0:
                 raise ValueError(f'Error: subtask {subtask.printTaskName} of '
                                  f'{taskTitle} did not set up successfully but also '
                                  'did not indicate a failure.  This likely indicates '
                                  'a bug like multiple tasks with the same name.')
            # a subtask failed setup_and_check
            print(f'Warning: subtask of {taskTitle} failed during check, '
                  'so this task will not be run')
            analysisTask._setupStatus = 'fail'
            totalFailures += 1
            return totalFailures

    analysesToGenerate[key] = analysisTask
    analysisTask._setupStatus = 'success'
    assert(totalFailures == 0)
    return totalFailures


def update_generate(config, generate):
    """
    Update the 'generate' config option using a string from the command line.

    Parameters
    ----------
    config : mpas_tools.config.MpasConfigParser
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
    config.set('output', 'generate', generateString, user=True)


def run_analysis(config, analyses):
    """
    Run all the tasks, either in serial or in parallel

    Parameters
    ----------
    config : mpas_tools.config.MpasConfigParser
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
    maxTitleLength = config.getint('plot', 'maxTitleLength')

    if len(mainRunName) > maxTitleLength:
        print('Warning: The main run name is quite long and will be '
              'truncated in some plots: \n{}\n\n'.format(mainRunName))

    configFileName = '{}/complete.{}.cfg'.format(logsDirectory, mainRunName)

    configFile = open(configFileName, 'w')
    config.write(configFile)
    configFile.close()

    parallelTaskCount = config.getint('execute', 'parallelTaskCount')

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
    logger.propagate = False

    totalTaskCount = len(analyses)
    widgets = ['Running tasks: ', progressbar.Percentage(), ' ',
               progressbar.Bar(), ' ', progressbar.ETA()]
    progress = progressbar.ProgressBar(widgets=widgets,
                                       max_value=totalTaskCount).start()

    runningProcessCount = 0

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

        progress.update(totalTaskCount - unfinishedCount)

        if unfinishedCount <= 0 and runningProcessCount == 0:
            # we're done
            break

        # launch new tasks
        runDirectly = False
        for key, analysisTask in analyses.items():
            if analysisTask._runStatus.value == AnalysisTask.READY:
                if isParallel:
                    newProcessCount = runningProcessCount + \
                        analysisTask.subprocessCount
                    if newProcessCount > parallelTaskCount and \
                            runningProcessCount > 0:
                        # this task should run next but we need to wait for
                        # more processes to finish
                        break

                    logger.info('Running {}'.format(
                        analysisTask.printTaskName))
                    if analysisTask.runDirectly:
                        analysisTask.run(writeLogFile=True)
                        runDirectly = True
                        break
                    else:
                        analysisTask._runStatus.value = AnalysisTask.RUNNING
                        analysisTask.start()
                        runningTasks[key] = analysisTask
                        runningProcessCount = newProcessCount
                        if runningProcessCount >= parallelTaskCount:
                            # don't try to run any more tasks
                            break
                else:
                    analysisTask.run(writeLogFile=False)
                    break

        if isParallel:

            if not runDirectly:
                assert(runningProcessCount > 0)
                # wait for a task to finish
                analysisTask = wait_for_task(runningTasks)
                key = (analysisTask.taskName, analysisTask.subtaskName)
                runningTasks.pop(key)
                runningProcessCount -= analysisTask.subprocessCount

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

    handler.close()
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


def wait_for_task(runningTasks, timeout=0.1):
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
                return analysisTask


def purge_output(config):
    outputDirectory = config.get('output', 'baseDirectory')
    if not os.path.exists(outputDirectory):
        print('Output directory {} does not exist.\n'
              'No purge necessary.'.format(outputDirectory))
    else:
        for subdirectory in ['plots', 'logs', 'mpasClimatology', 'mapping',
                             'timeSeries', 'html', 'mask', 'profiles',
                             'histogram']:
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


def build_config(user_config_file, shared_configs, machine_info):
    """
    Create a config parser from a user config file (either main or control)
    and a set of shared config file, also adding the username to the web_portal
    section
    """
    if not os.path.exists(user_config_file):
        raise OSError(f'A config file {user_config_file} was specified but '
                      f'the file does not exist')
    config = MpasConfigParser()
    for config_file in shared_configs:
        if config_file.endswith('.py'):
            # we'll skip config options set in python files
            continue
        config.add_from_file(config_file)
    config.add_user_config(user_config_file)

    if machine_info is not None:
        config.set('web_portal', 'username', machine_info.username)

    return config


def symlink_main_run(config, shared_configs, machine_info):
    """
    Create symlinks to the climatology and time-series directories for the
    main run that has already been computed so we don't have to recompute
    the analysis.
    """

    def link_dir(section, option):
        dest_directory = build_config_full_path(config=config,
                                                section='output',
                                                relativePathOption=option,
                                                relativePathSection=section)
        if not os.path.exists(dest_directory):

            source_directory = build_config_full_path(
                config=main_config, section='output',
                relativePathOption=option, relativePathSection=section)

            if os.path.exists(source_directory):

                dest_base = os.path.split(dest_directory)[0]

                make_directories(dest_base)

                os.symlink(source_directory, dest_directory)

    main_config_file = config.get('runs', 'mainRunConfigFile')
    main_config = build_config(main_config_file, shared_configs, machine_info)

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


def get_editable_install_dir(package_name):
    """
    Get the directory that the package is installed in if it is installed in
    editable mode, or None if it is not.

    Parameters
    ----------
    package_name : str
        The name of the package

    Returns
    -------
    install_dir : str or None
        The directory the package is installed in if in editable mode, or None
    """

    direct_url = Distribution.from_name(package_name).read_text(
        'direct_url.json')
    contents = json.loads(direct_url)
    pkg_is_editable = contents.get("dir_info", {}).get("editable", False)
    if pkg_is_editable and 'url' in contents:
        url = contents['url']
        if url.startswith('file://'):
            return url[7:]
    return None


def is_mpas_analysis_git_base():
    """
    Check if the current working directory is the base of an mpas_analysis git
    branch or a git worktree.

    Returns
    -------
    is_git_base : bool
        True if the current working directory is the base of an mpas_analysis
        git branch or a git worktree, False otherwise
    """
    mpas_analysis_dir = os.path.join(os.getcwd(), 'mpas_analysis')
    if not os.path.isdir(mpas_analysis_dir):
        # no package mpas_analysis, so can't be an mpas_analysis git base
        return False

    git_dir = os.path.join(os.getcwd(), '.git')
    if os.path.isdir(git_dir):
        # It's a git repository
        head_file = os.path.join(git_dir, 'HEAD')
    elif os.path.isfile(git_dir):
        # It's a git worktree
        with open(git_dir, 'r') as f:
            git_dir_path = f.read().strip().split(': ')[1]
        head_file = os.path.join(git_dir_path, 'HEAD')
    else:
        return False

    if not os.path.isfile(head_file):
        return False

    with open(head_file, 'r') as f:
        head_content = f.read()
        if 'ref: refs/heads/' in head_content:
            return True

    return False


def main():
    """
    Entry point for the main script ``mpas_analysis``
    """

    mpas_analysis_dir = get_editable_install_dir('mpas_analysis')
    if is_mpas_analysis_git_base() and mpas_analysis_dir is not None:
        # mpas_analysis is installed in editable mode and this is the base
        # of an mpas_analysis git branch
        if os.path.abspath(mpas_analysis_dir) != os.getcwd():
            raise OSError(
"""
The current working directory is the base of an mpas_analysis git branch,
but the package is installed in editable mode in a different directory.
Please reinstall mpas_analysis in editable mode using:
    python -m pip install --no-deps --no-build-isolation -e .
"""
            )

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-v', '--version',
                        action='version',
                        version='mpas_analysis {}'.format(
                                mpas_analysis.version.__version__),
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
    parser.add_argument("config_file", metavar="CONFIG", type=str, nargs='*',
                        help="config file")
    parser.add_argument("--plot_colormaps", dest="plot_colormaps",
                        action='store_true',
                        help="Make a plot displaying all available colormaps")
    parser.add_argument("--verbose", dest="verbose", action='store_true',
                        help="Verbose error reporting during setup-and-check "
                             "phase")
    parser.add_argument("-m", "--machine", dest="machine",
                        help="The name of the machine for loading machine-"
                             "related config options", metavar="MACH")
    parser.add_argument("--polar_regions", dest="polar_regions",
                        action='store_true',
                        help="Include config options for analysis focused on "
                             "polar regions")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    config = MpasConfigParser()

    # add default.cfg to cover default not included in the config files
    # provided on the command line
    config.add_from_package('mpas_analysis', 'default.cfg')

    # Add config options for E3SM supported machines from the mache package
    machine = args.machine

    if machine is None and 'E3SMU_MACHINE' in os.environ:
        machine = os.environ['E3SMU_MACHINE']

    if machine is None:
        machine = discover_machine()

    if machine is not None:
        print(f'Detected E3SM supported machine: {machine}')
        try:
            config.add_from_package('mache.machines', f'{machine}.cfg')
        except FileNotFoundError:

            possible_machines = []
            machine_configs = contents('mache.machines')
            for config in machine_configs:
                if config.endswith('.cfg'):
                    possible_machines.append(os.path.splitext(config)[0])

            possible_machines = '\n  '.join(sorted(possible_machines))
            raise ValueError(
                f'We could not find the machine: {machine}.\n'
                f'Possible machines are:\n  {possible_machines}')

        try:
            config.add_from_package('mpas_analysis.configuration',
                                    f'{machine}.cfg')
        except FileNotFoundError:
            # we don't have a config file for this machine, so we'll just
            # skip it.
            print(f'Warning: no MPAS-Analysis config file found for machine:'
                  f' {machine}')
            pass

    if args.polar_regions:
        config.add_from_package('mpas_analysis', 'polar_regions.cfg')

    if machine is not None:
        # set the username so we can use it in the htmlSubdirectory
        machine_info = MachineInfo(machine=machine)
        config.set('web_portal', 'username', machine_info.username)
    else:
        machine_info = None

    shared_configs = config.list_files()

    for user_config in args.config_file:
        if not os.path.exists(user_config):
            raise OSError(f'Config file {user_config} not found.')

        config.add_user_config(user_config)

    print('Using the following config files:')
    for config_file in config.list_files():
        print(f'   {config_file}')

    if args.list:
        # set this config option so we don't have issues
        config.set('diagnostics', 'baseDirectory', '')
        analyses = build_analysis_list(config, controlConfig=None)
        for analysisTask in analyses:
            print('task: {}'.format(analysisTask.taskName))
            print('    component: {}'.format(analysisTask.componentName)),
            print('    tags: {}'.format(', '.join(analysisTask.tags)))
        sys.exit(0)

    if args.plot_colormaps:
        register_custom_colormaps()
        _plot_color_gradients()
        sys.exit(0)

    if config.has_option('runs', 'controlRunConfigFile'):
        control_config_file = config.get('runs', 'controlRunConfigFile')
        control_config = build_config(control_config_file, shared_configs,
                                      machine_info)

        # replace the log directory so log files get written to this run's
        # log directory, not the control run's
        logs_dir = build_config_full_path(config, 'output', 'logsSubdirectory')

        control_config.set('output', 'logsSubdirectory', logs_dir)

        print('Comparing to control run {} rather than observations. \n'
              'Make sure that MPAS-Analysis has been run previously with the '
              'control config file.'.format(control_config.get('runs',
                                                               'mainRunName')))
    else:
        control_config = None

    if args.purge:
        purge_output(config)

    if config.has_option('runs', 'mainRunConfigFile'):
        symlink_main_run(config, shared_configs, machine_info)

    if args.generate:
        update_generate(config, args.generate)

    if control_config is not None:
        # we want to use the "generate" option from the current run, not
        # the control config file
        control_config.set('output', 'generate', config.get('output',
                                                            'generate'))

    log_dir = build_config_full_path(config, 'output', 'logsSubdirectory')
    make_directories(log_dir)

    update_time_bounds_in_config(config)

    file_cache_maxsize = config.getint('input', 'file_cache_maxsize')
    try:
        xarray.set_options(file_cache_maxsize=file_cache_maxsize)
    except ValueError:
        # xarray version doesn't support file_cache_maxsize yet...
        pass

    start_time = time.time()

    custom_config_files = list(args.config_file)
    for option in ['controlRunConfigFile', 'mainRunConfigFile']:
        if config.has_option('runs', option):
            custom_config_files.append(config.get('runs', option))

    html_base_directory = build_config_full_path(config, 'output',
                                                 'htmlSubdirectory')
    make_directories(html_base_directory)
    for config_filename in custom_config_files:
        config_filename = os.path.abspath(config_filename)
        print(f'copying {config_filename} to HTML dir.')
        basename = os.path.basename(config_filename)
        copyfile(config_filename, f'{html_base_directory}/{basename}')

    analyses = build_analysis_list(config, control_config)
    analyses = determine_analyses_to_generate(analyses, args.verbose)

    setup_duration = time.time() - start_time

    if not args.setup_only and not args.html_only:
        run_analysis(config, analyses)
        run_duration = time.time() - start_time
        m, s = divmod(setup_duration, 60)
        h, m = divmod(int(m), 60)
        print('Total setup time: {}:{:02d}:{:05.2f}'.format(h, m, s))
        m, s = divmod(run_duration, 60)
        h, m = divmod(int(m), 60)
        print('Total run time: {}:{:02d}:{:05.2f}'.format(h, m, s))

    if not args.setup_only:
        generate_html(config, analyses, control_config, custom_config_files)


if __name__ == "__main__":
    main()
