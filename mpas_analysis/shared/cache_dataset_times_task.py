from .analysis_task import AnalysisTask

from ..shared.io import StreamsFile, build_config_full_path
from ..shared.timekeeping.utility import get_simulation_start_time


class CacheDatasetTimesTask(AnalysisTask):  # {{{
    '''
    A task for caching the times in a multifile data sets of an MPAS analysis
    member for later use.  Since analysis member may be used by multiple tasks
    (indeed the ``timeSeriesStats`` is currently used by *all* tasks), it is
    important that this time information gets processed once before all other
    tasks get run in parallel.

    Authors
    -------
    Xylar Asay-Davis
    '''

    def __init__(self, config, componentName, streamName,
                 startAndEndDateSections, namelistOption=None):  # {{{
        '''
        Construct an analysis task for caching the times in the multifile
        data set in the given component and stream.  The name of the task
        includes the component and stream name.  For example, if
        ``component='ocean'`` and ``streamName='timeSeriesStats``, then the
        task name is ``cacheOceanTimeSeriesStatsTimes``.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        componentName :  {'ocean', 'seaIce'}
            The name of the component (same as the folder where the task
            resides)

        streamName : str
            The name of the stream from which the climatology data set will
            be read, used to cache times and update their bounds in the
            configuration parser.

        startAndEndDateSections : list, {'climatology', 'timeSeries', 'index'}
            The name of sections in the config file containing ``startDate``
            and ``endDate`` options as a list.

        namelistOption : str, optional
            The name of a namelist option (e.g.
            ``config_am_timeseriesstatsmonthly_enable``) that should be set to
            true of the required analysis member has been enabled.  If this
            option is ``None`` (the default), no check is performed

        Authors
        -------
        Xylar Asay-Davis
        '''

        upperComponent = componentName[0].upper() + componentName[1:]
        upperStream = streamName[0].upper() + streamName[1:]

        taskName = 'cache{}{}Times'.format(
                upperComponent, upperStream)

        # first, call the constructor from the base class (AnalysisTask).
        super(CacheDatasetTimesTask, self).__init__(
                config=config,
                taskName=taskName,
                componentName=componentName,
                tags=startAndEndDateSections)

        self.streamName = streamName
        self.namelistOption = namelistOption
        self.startAndEndDateSections = startAndEndDateSections

        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError: if startAndEndDateSections is not a list containing
            {'climatology', 'timeSeries', 'index'}

        Authors
        -------
        Xylar Asay-Davis
        '''

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(CacheDatasetTimesTask, self).setup_and_check()

        if self.namelistOption is not None:
            self.check_analysis_enabled(
                analysisOptionName=self.namelistOption,
                raiseException=True)

        for section in self.startAndEndDateSections:
            if not self.config.has_section(section):
                raise ValueError('Config file does not have a section '
                                 '{}.'.format(section))
            for option in ['startDate', 'endDate']:
                if not self.config.has_option(section, option):
                    raise ValueError('Config section {} does not have '
                                     'expected option {}.'.format(section,
                                                                  option))

        # }}}

    def run(self):  # {{{
        '''
        The main method of the task that performs the analysis task.

        Authors
        -------
        Xylar Asay-Davis
        '''
        print ""
        print "Caching times from the {} component for the {} stream " \
              "...".format(self.componentName, self.streamName)

        try:
            self.simulationStartTime = get_simulation_start_time(
                self.runStreams)
        except IOError as e:
            if self.componentName == 'ocean':
                raise e
            else:
                # try the ocean stream instead
                runDirectory = build_config_full_path(self.config, 'input',
                                                      'runSubdirectory')
                oceanStreamsFileName = build_config_full_path(
                    self.config, 'input', 'oceanStreamsFileName')
                oceanStreams = StreamsFile(oceanStreamsFileName,
                                           streamsdir=runDirectory)
                self.simulationStartTime = \
                    get_simulation_start_time(oceanStreams)

        for sectionName in self.startAndEndDateSections:
            # internally, this will also cache the times for the data set
            self.update_start_end_date(section=sectionName,
                                       streamName=self.streamName)

        # }}}

# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
