from ..shared.analysis_task import AnalysisTask

from ..shared.io import StreamsFile
from ..shared.io.utility import build_config_full_path

from ..shared.timekeeping.utility import get_simulation_start_time


class SeaIceAnalysisTask(AnalysisTask):  # {{{
    '''
    A base class for sea-ice analysis tasks that implements common setup

    Authors
    -------
    Xylar Asay-Davis
    '''

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up common to sea ice analysis, including
        finding the simulation start time and an appropriate restart file
        (from MPAS-SeaIce if possible, but from MPAS-Ocean as a fallback)

        Authors
        -------
        Xylar Asay-Davis
        '''

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(SeaIceAnalysisTask, self).setup_and_check()

        try:
            self.simulationStartTime = get_simulation_start_time(
                self.runStreams)
        except IOError:
            # try the ocean stream instead
            runDirectory = build_config_full_path(self.config, 'input',
                                                  'runSubdirectory')
            oceanStreamsFileName = build_config_full_path(
                self.config, 'input', 'oceanStreamsFileName')
            oceanStreams = StreamsFile(oceanStreamsFileName,
                                       streamsdir=runDirectory)
            self.simulationStartTime = get_simulation_start_time(oceanStreams)

        try:
            self.restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            # get an ocean restart file, since no sea-ice restart exists
            try:
                runDirectory = build_config_full_path(self.config, 'input',
                                                      'runSubdirectory')
                oceanStreamsFileName = build_config_full_path(
                    self.config, 'input', 'oceanStreamsFileName')
                oceanStreams = StreamsFile(oceanStreamsFileName,
                                           streamsdir=runDirectory)
                self.restartFileName = oceanStreams.readpath('restart')[0]
            except ValueError:
                raise IOError('No MPAS-O or MPAS-Seaice restart file found: '
                              'need at least one restart file \n'
                              'for sea ice analysis tasks')

        # }}}


# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
