'''
common utility functions for analysis tasks

Xylar Asay-Davis

Last Modified: 03/23/2017
'''

import warnings

from .io import NameList, StreamsFile
from .io.utility import build_config_full_path, make_directories

from ..ocean.variable_stream_map import oceanNamelistMap, oceanStreamMap, \
    oceanVariableMap

from ..sea_ice.variable_stream_map import seaIceNamelistMap, seaIceStreamMap, \
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

    namelistMap : dict
        A dictionary of MPAS namelist options that map to their mpas_analysis
        counterparts.

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
        namelistMap = oceanNamelistMap
        streamMap = oceanStreamMap
        variableMap = oceanVariableMap
    elif componentName == 'seaIce':
        namelistMap = seaIceNamelistMap
        streamMap = seaIceStreamMap
        variableMap = seaIceVariableMap
    else:
        namelistMap = None
        streamMap = None
        variableMap = None

    return namelist, runStreams, historyStreams, calendar, namelistMap, \
        streamMap, variableMap, plotsDirectory  # }}}


def check_analysis_enabled(namelist, analysisOptionName, namelistMap=None,
                           default=False, raiseException=True):
    '''
    Check to make sure a given analysis is turned on, issuing a warning or
    raising an exception if not.

    Parameters
    ----------
    namelist : NameList object
        for parsing namelist options

    analysisOptionName : str
        The name of a boolean namelist option indicating whether the given
        analysis member is enabled

    namelistMap : dict, optional
        A dictionary of MPAS namelist options that map to their mpas_analysis
        counterparts.

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

    Last Modified
    -------------
    04/01/2017
    '''

    try:
        if namelistMap is None:
            optionName = analysisOptionName
        else:
            optionName = namelist.find_option(namelistMap[analysisOptionName])
        enabled = namelist.getbool(optionName)
    except ValueError:
        enabled = default
        if default:
            message = 'WARNING: namelist option {} not found.\n' \
                      'This likely indicates that the simulation you are ' \
                      'analyzing was run with an\n' \
                      'older version of MPAS-O that did not support this ' \
                      'flag.  Assuming enabled.'.format(analysisOptionName)
            warnings.warn(message)

    if not enabled and raiseException:
        raise RuntimeError('*** MPAS-Analysis relies on {} = .true.\n'
                           '*** Make sure to enable this analysis '
                           'member.'.format(analysisOptionName))

    return enabled

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
