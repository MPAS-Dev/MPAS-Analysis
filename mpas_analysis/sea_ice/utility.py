'''
common utility functions for sea ice analysis tasks

Xylar Asay-Davis

Last Modified: 04/03/2017
'''

from ..shared.io import StreamsFile
from ..shared.io.utility import build_config_full_path

from ..shared.timekeeping.utility import get_simulation_start_time

from ..shared.analysis_task import setup_task


def setup_sea_ice_task(config):  # {{{
    '''
    Perform steps to set up the sea ice analysis (e.g. reading namelists and
    streams files, finding a restart file).

    Parameters
    ----------
    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    Returns
    -------
    namelist : NameList object
        for parsing namelist options

    runStreams : StreamsFile object
        for parsing the streams file related to output in the run subdirectory

    historyStreams : StreamsFile object
        for parsing the streams file related to output in the history
        subdirectory

    calendar: {'gregorian', 'gregorian_noleap'}
        The name of the calendars used in the MPAS run

    streamMap : dict
        A dictionary of MPAS stream names that map to their mpas_analysis
        counterparts.

    variableMap : dict
        A dictionary of MPAS variable names that map to their mpas_analysis
        counterparts.

    plotsDirectory : str
        the directories for writing plots

    simulationStartTime : str
        the date and time of the start of the simulation

    restartFileName : str
        the path to a restart file (which may be an MPAS-O restart file is no
        MPAS-SeaIce restart could be found)

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    04/03/2017
    '''
    # perform common setup for the task
    namelist, runStreams, historyStreams, calendar, streamMap, \
        variableMap, plotsDirectory = setup_task(config,
                                                 componentName='seaIce')

    try:
        simulationStartTime = get_simulation_start_time(runStreams)
    except IOError:
        # try the ocean stream instead
        runDirectory = build_config_full_path(config, 'input',
                                              'runSubdirectory')
        oceanStreamsFileName = build_config_full_path(
            config, 'input', 'oceanStreamsFileName')
        oceanStreams = StreamsFile(oceanStreamsFileName,
                                   streamsdir=runDirectory)
        simulationStartTime = get_simulation_start_time(oceanStreams)

    try:
        restartFileName = runStreams.readpath('restart')[0]
    except ValueError:
        # get an ocean restart file, since no sea-ice restart exists
        try:
            runDirectory = build_config_full_path(config, 'input',
                                                  'runSubdirectory')
            oceanStreamsFileName = build_config_full_path(
                config, 'input', 'oceanStreamsFileName')
            oceanStreams = StreamsFile(oceanStreamsFileName,
                                       streamsdir=runDirectory)
            restartFileName = oceanStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O or MPAS-Seaice restart file found: need '
                          'at least one restart file for seaice_timeseries '
                          'calculation')

    return namelist, runStreams, historyStreams, calendar, streamMap, \
        variableMap, plotsDirectory, simulationStartTime, restartFileName
    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python