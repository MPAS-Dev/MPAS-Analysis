import os
import warnings
import subprocess
from distutils.spawn import find_executable
import xarray as xr

from ..analysis_task import AnalysisTask

from ..io.utility import build_config_full_path, make_directories


class MpasTimeSeriesTask(AnalysisTask):  # {{{
    '''
    An analysis tasks for computing time series from output from the
    ``timeSeriesStatsMonthly`` analysis member.

    Attributes
    ----------

    variableList : list of str
        A list of variable names in ``timeSeriesStatsMonthly`` to be
        included in the time series

    inputFiles : list of str
        A list of input files from which to extract the time series.

    startDate, endDate : str
        The start and end dates of the time series as strings

    startYear, endYear : int
        The start and end years of the time series

    Authors
    -------
    Xylar Asay-Davis
    '''

    def __init__(self, config, componentName, taskName=None,
                 subtaskName=None):  # {{{
        '''
        Construct the analysis task for extracting time series.

        Parameters
        ----------
        config : ``MpasAnalysisConfigParser``
            Contains configuration options

        componentName : {'ocean', 'seaIce'}
            The name of the component (same as the folder where the task
            resides)

        taskName : str, optional
            The name of the task, 'mpasTimeSeriesOcean' or
            'mpasTimeSeriesSeaIce' by default (depending on ``componentName``)

        subtaskName : str, optional
            The name of the subtask (if any)


        Authors
        -------
        Xylar Asay-Davis
        '''
        self.variableList = []
        self.seasons = []

        tags = ['timeSeries']

        suffix = componentName[0].upper() + componentName[1:]
        if taskName is None:
            taskName = 'mpasTimeSeries{}'.format(suffix)

        # call the constructor from the base class (AnalysisTask)
        super(MpasTimeSeriesTask, self).__init__(
            config=config,
            taskName=taskName,
            subtaskName=subtaskName,
            componentName=componentName,
            tags=tags)

        # }}}

    def add_variables(self, variableList):  # {{{
        '''
        Add one or more variables to extract as a time series.

        Parameters
        ----------
        variableList : list of str
            A list of variable names in ``timeSeriesStatsMonthly`` to be
            included in the time series

        Authors
        -------
        Xylar Asay-Davis
        '''

        for variable in variableList:
            if variable not in self.variableList:
                self.variableList.append(variable)

        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        '''

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(MpasTimeSeriesTask, self).setup_and_check()

        config = self.config
        baseDirectory = build_config_full_path(
            config, 'output', 'timeSeriesSubdirectory')

        make_directories(baseDirectory)

        self.outputFile = '{}/{}.nc'.format(baseDirectory,
                                            self.fullTaskName)

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        # get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        startDate = config.get('timeSeries', 'startDate')
        endDate = config.get('timeSeries', 'endDate')
        streamName = 'timeSeriesStatsMonthlyOutput'
        self.inputFiles = self.historyStreams.readpath(
                streamName, startDate=startDate, endDate=endDate,
                calendar=self.calendar)

        if len(self.inputFiles) == 0:
            raise IOError('No files were found in stream {} between {} and '
                          '{}.'.format(streamName, startDate, endDate))

        self._update_time_series_bounds_from_file_names()

        # }}}

    def run_task(self):  # {{{
        '''
        Compute the requested time series

        Authors
        -------
        Xylar Asay-Davis
        '''

        if len(self.variableList) == 0:
            # nothing to do
            return

        self.logger.info('\nComputing MPAS time series from files:\n'
                         '    {} through\n    {}'.format(
                                 os.path.basename(self.inputFiles[0]),
                                 os.path.basename(self.inputFiles[-1])))

        self._compute_time_series_with_ncrcat()

        # }}}

    def _update_time_series_bounds_from_file_names(self):  # {{{
        """
        Update the start and end years and dates for time series based on the
        years actually available in the list of files.

        Authors
        -------
        Xylar Asay-Davis
        """

        config = self.config
        section = 'timeSeries'

        requestedStartYear = config.getint(section, 'startYear')
        requestedEndYear = config.getint(section, 'endYear')

        dates = sorted([fileName[-13:-6] for fileName in self.inputFiles])
        years = [int(date[0:4]) for date in dates]
        months = [int(date[5:7]) for date in dates]

        # search for the start of the first full year
        firstIndex = 0
        while(firstIndex < len(years) and months[firstIndex] != 1):
            firstIndex += 1
        startYear = years[firstIndex]

        # search for the end of the last full year
        lastIndex = len(years)-1
        while(lastIndex >= 0 and months[lastIndex] != 12):
            lastIndex -= 1
        endYear = years[lastIndex]

        if startYear != requestedStartYear or endYear != requestedEndYear:
            message = "time series start and/or end year different from " \
                      "requested\n" \
                      "requestd: {:04d}-{:04d}\n" \
                      "actual:   {:04d}-{:04d}\n".format(requestedStartYear,
                                                         requestedEndYear,
                                                         startYear,
                                                         endYear)
            warnings.warn(message)
            config.set(section, 'startYear', str(startYear))
            config.set(section, 'endYear', str(endYear))

            startDate = '{:04d}-01-01_00:00:00'.format(startYear)
            config.set(section, 'startDate', startDate)
            endDate = '{:04d}-12-31_23:59:59'.format(endYear)
            config.set(section, 'endDate', endDate)
        else:
            startDate = config.get(section, 'startDate')
            endDate = config.get(section, 'endDate')

        self.startDate = startDate
        self.endDate = endDate
        self.startYear = startYear
        self.endYear = endYear

        # }}}

    def _compute_time_series_with_ncrcat(self):
        # {{{
        '''
        Uses ncrcat to extact time series from timeSeriesMonthlyOutput files

        Raises
        ------
        OSError
            If ``ncrcat`` is not in the system path.

        Author
        ------
        Xylar Asay-Davis
        '''

        if find_executable('ncrcat') is None:
            raise OSError('ncrcat not found. Make sure the latest nco '
                          'package is installed: \n'
                          'conda install nco\n'
                          'Note: this presumes use of the conda-forge '
                          'channel.')

        if os.path.exists(self.outputFile):
            # add only input files wiht times that aren't already in the
            # output file
            dates = sorted([fileName[-13:-6] for fileName in self.inputFiles])
            inYears = [int(date[0:4]) for date in dates]
            inMonths = [int(date[5:7]) for date in dates]
            totalMonths = 12*inYears + inMonths

            with xr.open_dataset(self.outputFile) as ds:
                lastDate = str(ds.xtime_startMonthly[-1].values)

            lastYear = int(lastDate[0:4])
            lastMonth = int(lastDate[5:7])
            lastTotalMonths = 12*lastYear + lastMonth

            inputFiles = []
            for index, inputFile in enumerate(self.inputFiles):
                if totalMonths[index] > lastTotalMonths:
                    inputFiles.append(inputFile)
		
            if len(inputFiles) == 0:
                # nothing to do
                return
        else:
            inputFiles = self.inputFiles

        variableList = self.variableList + ['xtime_startMonthly',
                                            'xtime_endMonthly']

        args = ['ncrcat', '--record_append', '--no_tmp_fl',
                '-v', ','.join(variableList)]

        args.extend(inputFiles)
        args.append(self.outputFile)

        self.logger.info('running: {}'.format(' '.join(args)))
        for handler in self.logger.handlers:
            handler.flush()

        process = subprocess.Popen(args, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        if stdout:
            self.logger.info(stdout)
        if stderr:
            for line in stderr.split('\n'):
                self.logger.error(line)

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode,
                                                ' '.join(args))

        # }}}
    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
