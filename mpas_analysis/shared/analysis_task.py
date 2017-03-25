'''
common utility functions for analysis tasks

Xylar Asay-Davis

Last Modified: 03/23/2017
'''

from .io import NameList, StreamsFile
from .io.utility import build_config_full_path, make_directories

from ..ocean.variable_stream_map import oceanStreamMap, \
    oceanVariableMap

from ..sea_ice.variable_stream_map import seaIceStreamMap, \
    seaIceVariableMap


def setup_task(config, componentName):  # {{{
    '''
    Perform steps to set up the analysis (e.g. reading namelists and
    streams files).

    Parameters
    ----------
    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    componentName : {'ocean', 'seaIce', 'landIce'}
        The name of a MPAS core to be analyized

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

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    03/23/2017
    '''

    # read parameters from config file
    # the run directory contains the restart files
    runDirectory = build_config_full_path(config, 'input', 'runSubdirectory')
    # if the history directory exists, use it; if not, fall back on
    # runDirectory
    historyDirectory = build_config_full_path(
        config, 'input', '{}HistorySubdirectory'.format(componentName),
        defaultPath=runDirectory)

    namelistFileName = build_config_full_path(
        config, 'input', '{}NamelistFileName'.format(componentName))
    namelist = NameList(namelistFileName)

    streamsFileName = build_config_full_path(
        config, 'input', '{}StreamsFileName'.format(componentName))
    runStreams = StreamsFile(streamsFileName, streamsdir=runDirectory)
    historyStreams = StreamsFile(streamsFileName, streamsdir=historyDirectory)

    calendar = namelist.get('config_calendar_type')

    plotsDirectory = build_config_full_path(config, 'output',
                                            'plotsSubdirectory')
    make_directories(plotsDirectory)

    if componentName == 'ocean':
        streamMap = oceanStreamMap
        variableMap = oceanVariableMap
    elif componentName == 'seaIce':
        streamMap = seaIceStreamMap
        variableMap = seaIceVariableMap
    else:
        streamMap = None
        variableMap = None

    return namelist, runStreams, historyStreams, calendar, streamMap, \
        variableMap, plotsDirectory  # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
