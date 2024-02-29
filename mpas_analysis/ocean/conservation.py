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
#
# Author
# -------
# Carolyn Begeman

from distutils.spawn import find_executable
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import xarray as xr

from mpas_analysis.shared.analysis_task import AnalysisTask
from mpas_analysis.shared.constants import constants
from mpas_analysis.shared.html import write_image_xml
from mpas_analysis.shared.io import open_mpas_dataset
from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories, get_files_year_month, decode_strings
from mpas_analysis.shared.plot import timeseries_analysis_plot, savefig
from mpas_analysis.shared.timekeeping.utility import date_to_days, \
    days_to_datetime


class ConservationTask(AnalysisTask):
    """
    This task generates time series plots from output from the conservation
    analysis member. A number of different plot types are supported, as indicated
    in the `plotTypes` config option in the `conservation` section.

    Attributes
    ----------
    config :  mpas_tools.config.MpasConfigParser
        Contains configuration options

    outputFile : str
        The path to the output file produced by this analysis

    runDirectory : str
        The path to the restart files from the main simulation being analyzed

    historyDirectory : str
        The path to the history files from the main simulation being analyzed

    startYear : int
        The year to start the analysis

    endYear : int
        The year to end the analysis

    inputFiles : list of str
        The paths to all conservation AM files

    mainRunName : str
        The name of the main run from the config file

    referenceRunName : str
        The name of the reference run from the config file

    referenceInputDirectory : str
        The directory of the reference run from the config file

    plotTypes : list of str
        The plot types requested in the config file

    masterVariableList : dict of key-[list of str] pairs
        Keys are the supported plot types. Entries are lists of the variables
        that are needed to produce that plot type.

    derivedVariableList : dict of key-[list of str] pairs
        Keys are the derived variables. Entries are lists of variables in the
        AM output that are needed to derive that variable.

    xmlFileNames : list of str
        File names for xml output with full path

    filePrefixes : list of str
        File prefixes for xml files

    variableList : dict of key-[list of str] pairs
        Keys are the requested plot types. Entries are lists of the variables
        in the AM output that are needed to produce that plot type.
    """

    # Authors
    # -------
    # Carolyn Begeman

    def __init__(self, config):
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  mpas_tools.config.MpasConfigParser
            Contains configuration options
        """
        # Authors
        # -------
        # Carolyn Begeman

        super(ConservationTask, self).__init__(
            config=config,
            taskName='oceanConservation',
            componentName='ocean',
            tags=['timeSeries', 'conservation'])

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.
        Also, adds attributes to the task that will be needed later.

        Raises
        ------
        ValueError: if plot type is not supported, conservation analysis member
            is inactive or input files are missing
        """
        # Authors
        # -------
        # Carolyn Begeman

        super(ConservationTask, self).setup_and_check()

        # Check that the conservation analysis member is active
        self.check_analysis_enabled(
            analysisOptionName='config_am_conservationcheck_enable',
            raiseException=True)

        # Specify where to put analysis task output
        config = self.config
        baseDirectory = build_config_full_path(
            config, 'output', 'conservationSubdirectory')
        make_directories(baseDirectory)
        self.outputFile = f'{baseDirectory}/{self.fullTaskName}.nc'

        # get a list of conservationCheck output files from the streams file,
        # reading only those that are between the start and end dates

        # the run directory contains the restart files
        self.runDirectory = build_config_full_path(self.config, 'input',
                                                   'runSubdirectory')
        # if the history directory exists, use it; if not, fall back on
        # runDirectory
        self.historyDirectory = build_config_full_path(
            self.config, 'input',
            f'{self.componentName}HistorySubdirectory',
            defaultPath=self.runDirectory)

        self.startYear = self.config.getint('timeSeries', 'startYear')
        self.endYear = self.config.getint('timeSeries', 'endYear')
        self.inputFiles = sorted(self.historyStreams.readpath(
            'conservationCheckOutput',
            startDate=f'{self.startYear:04d}-01-01_00:00:00',
            endDate=f'{self.endYear:04d}-01-01_00:00:00',
            calendar=self.calendar))

        if len(self.inputFiles) == 0:
            raise IOError(f'No files were found matching {self.inputFiles}')

        with xr.open_dataset(self.inputFiles[0]) as ds:
            self.allVariables = list(ds.data_vars.keys())

        self.mainRunName = self.config.get('runs', 'mainRunName')
        self.referenceRunName = \
            config.get('runs', 'preprocessedReferenceRunName')
        self.referenceInputDirectory = config.get('oceanPreprocessedReference',
                                                  'baseDirectory')

        self.plotTypes = self.config.getexpression('timeSeriesConservation', 'plotTypes')

        self.masterVariableList = {'absolute_energy_error': ['absoluteEnergyError'],
                                   'total_energy_flux': ['netEnergyFlux'],
                                   'absolute_salt_error': ['absoluteSaltError'],
                                   'ice_salt_flux': ['netSaltFlux'],
                                   'total_mass_flux': ['netMassFlux'],
                                   'total_mass_change': ['netMassChange'],
                                   'land_ice_mass_change': ['landIceMassChange'],
                                   'land_ice_ssh_change': ['landIceSshChange'],
                                   'land_ice_mass_flux': ['landIceMassFlux'],
                                   'land_ice_mass_flux_components': ['accumulatedIcebergFlux',
                                                                'accumulatedLandIceFlux',
                                                                'accumulatedRemovedRiverRunoffFlux',
                                                                'accumulatedRemovedIceRunoffFlux']}

        # for each derived variable, which source variables are needed
        self.derivedVariableList = {'netMassChange': ['massChange'],
                                    'landIceMassFlux': ['accumulatedIcebergFlux',
                                                        'accumulatedLandIceFlux',
                                                        'accumulatedRemovedRiverRunoffFlux',
                                                        'accumulatedRemovedIceRunoffFlux'],
                                    'landIceSshChange': ['accumulatedIcebergFlux',
                                                         'accumulatedLandIceFlux',
                                                         'accumulatedRemovedRiverRunoffFlux',
                                                         'accumulatedRemovedIceRunoffFlux'],
                                    'landIceMassChange': ['accumulatedIcebergFlux',
                                                          'accumulatedLandIceFlux',
                                                          'accumulatedRemovedRiverRunoffFlux',
                                                          'accumulatedRemovedIceRunoffFlux']}

        # Determine the xml files for each plot and the variables each plot will use
        self.xmlFileNames = []
        self.filePrefixes = {}
        self.variableList = {}
        for plot_type in self.plotTypes:
            if plot_type not in self.masterVariableList.keys():
                raise ValueError(f'plot type {plot_type} not supported')
            filePrefix = f'conservation_{self.mainRunName}_{plot_type}_' \
                         f'years{self.startYear:04d}-{self.endYear:04d}'
            self.xmlFileNames.append(f'{self.plotsDirectory}/{filePrefix}.xml')
            self.filePrefixes[plot_type] = filePrefix
            self.variableList[plot_type] = self._add_variables(self.masterVariableList[plot_type])

    def run_task(self):
        """
        Create an output netCDF file that has all of the requested conservation AM variables
        in the requested time window. Then generate all requested conservation plot types.
        """
        # Authors
        # -------
        # Carolyn Begeman

        all_plots_variable_list = []
        for plot_type in self.plotTypes:
            for varname in self.variableList[plot_type]:
                all_plots_variable_list.append(varname)
        self._compute_time_series_with_ncrcat(all_plots_variable_list)
        for plot_type in self.plotTypes:
            self._make_plot(plot_type)

    def _add_variables(self, target_variable_list):
        """
        Add one or more variables to extract as a time series.

        Parameters
        ----------
        variableList : list of str
            A list of variable names in ``conservationCheck`` to be
            included in the time series

        Raises
        ------
        ValueError
            if this function is called before this task has been set up (so
            the list of available variables has not yet been set) or if one
            or more of the requested variables is not available in the
            ``conservationCheck`` output.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        variable_list = []
        if self.allVariables is None:
            raise ValueError('add_variables() can only be called after '
                             'setup_and_check() in ConservationTask.\n'
                             'Presumably tasks were added in the wrong order '
                             'or add_variables() is being called in the wrong '
                             'place.')

        for variable in target_variable_list:
            if variable not in self.allVariables and \
                        variable not in self.derivedVariableList.keys():
                raise ValueError(
                    f'{variable} is not available in conservationCheck'
                    'output:\n{self.allVariables}')

            if variable in self.allVariables and variable not in variable_list:
                variable_list.append(variable)
            # If it's a derived variable, add all of the variables it depends on
            if variable in self.derivedVariableList.keys() and \
                   variable not in variable_list:
                for var in self.derivedVariableList[variable]:
                    variable_list.append(var)

        return variable_list

    def _make_plot(self, plot_type):
        """
        Generate time series plots from conservation AM output.

        Parameters
        ----------
        plot_type: str
            The type of plot to generate from conservationCheck variables
        """
        config = self.config
        filePrefix = self.filePrefixes[plot_type]
        outFileName = f'{self.plotsDirectory}/{filePrefix}.png'

        titles = {}
        titles['total_energy_flux'] = 'Total energy flux'
        titles['absolute_energy_error'] = 'Energy error'
        titles['ice_salt_flux'] = 'Salt flux related to land ice and sea ice'
        titles['absolute_salt_error'] = 'Salt conservation error'
        titles['total_mass_flux'] = 'Total mass flux'
        titles['total_mass_change'] = 'Total mass anomaly'
        titles['land_ice_mass_flux'] = 'Mass flux due to land ice'
        titles['land_ice_mass_change'] = 'Mass anomaly due to land ice fluxes'
        titles['land_ice_ssh_change'] = 'SSH anomaly due to land ice fluxes'
        titles['land_ice_mass_flux_components'] = 'Mass fluxes from land ice'

        y_labels = {}
        y_labels['total_energy_flux'] = 'Energy flux (W)'
        y_labels['absolute_energy_error'] = 'Energy (J)'
        y_labels['ice_salt_flux'] = 'Salt flux (Gt/yr)'
        y_labels['absolute_salt_error'] = 'Salt (Gt)'
        y_labels['total_mass_flux'] = 'Mass flux (Gt/yr)'
        y_labels['total_mass_change'] = 'Mass (Gt)'
        y_labels['land_ice_mass_flux'] = 'Mass flux (Gt/yr)'
        y_labels['land_ice_mass_change'] = 'Mass (Gt)'
        y_labels['land_ice_ssh_change'] = 'SSH anomaly (mm)'
        y_labels['land_ice_mass_flux_components'] = 'Mass flux (Gt/yr)'

        captions = {}
        captions['total_energy_flux'] = 'Total energy flux'
        captions['absolute_energy_error'] = 'Absolute energy conservation error'
        captions['ice_salt_flux'] = 'Salt flux related to land ice and sea ice ' \
            '(sea ice salinity flux, sea ice frazil flux, and land ice frazil flux)'
        captions['absolute_salt_error'] = 'Absolute salt conservation error'
        captions['total_mass_flux'] = 'Total mass flux'
        captions['total_mass_change'] = 'Total mass anomaly'
        captions['land_ice_mass_flux'] = 'Mass flux due to land ice'
        captions['land_ice_mass_change'] = 'Mass anomaly due to land ice fluxes'
        captions['land_ice_ssh_change'] = 'SSH anomaly due to land ice fluxes. Assumes a constant ocean area.'
        captions['land_ice_mass_flux_components'] = 'Mass flux components from land ice'

        self.logger.info(f'  Open conservation file {self.outputFile}...')
        ds = open_mpas_dataset(fileName=self.outputFile,
                               calendar=self.calendar,
                               variableList=self.variableList[plot_type],
                               timeVariableNames='xtime',
                               startDate=f'{self.startYear:04d}-01-01_00:00:00',
                               endDate=f'{self.endYear:04d}-01-01_00:00:00')

        if self.referenceRunName != 'None' and self.referenceInputDirectory != 'None':
            inFilesPreprocessed = f'{self.referenceInputDirectory}/timeseries/{self.fullTaskName}.nc'
            self.logger.info('  Load in conservation for a preprocessed reference '
                             f'run {inFilesPreprocessed}...')
            ds_ref = open_mpas_dataset(fileName=inFilesPreprocessed,
                                   calendar=self.calendar,
                                   variableList=self.variableList[plot_type],
                                   timeVariableNames='xtime')
            yearEndPreprocessed = days_to_datetime(ds_ref.Time.max(),
                                                   calendar=self.calendar).year
            if self.startYear <= yearEndPreprocessed:
                timeStart = date_to_days(year=self.startYear, month=1, day=1,
                                         calendar=self.calendar)
                timeEnd = date_to_days(year=self.endYear, month=12, day=31,
                                       calendar=self.calendar)
                ds_ref_slice= \
                    ds_ref.sel(Time=slice(timeStart, timeEnd))
            else:
                self.logger.warning('Preprocessed time series ends before the '
                                    'timeSeries startYear and will not be '
                                    'plotted.')
                self.referenceRunName = 'None'

        # make the plot
        self.logger.info('  Make conservation plots...')
        xLabel = 'Time (years)'
        title = titles[plot_type]
        yLabel = y_labels[plot_type]
        lineStylesBase = ['-', '--', '-.', ':']

        # gather all the variables for this plot type
        fields = []
        legendText = []
        lineColors = []
        lineStyles = []
        for index, varname in enumerate(self.masterVariableList[plot_type]):
            variable = self._get_variable(ds, varname)
            fields.append(variable)
            legend_text = ''
            if self.referenceRunName != 'None':
                legend_text = self.mainRunName
            if len(self.masterVariableList[plot_type]) > 1:
                if len(legend_text) > 0:
                    legend_text = f'{legend_text}, '
                legend_text = f"{legend_text}{varname.replace('accumulated', '').replace('Flux', '')}"
            legendText.append(legend_text)
            lineColors.append(config.get('timeSeries', 'mainColor'))
            lineStyles.append(lineStylesBase[index])
            if self.referenceRunName != 'None':
                variable = self._get_variable(ds, varname)
                fields.append(variable)
                legend_text = self.referenceRunName
                if len(self.masterVariableList[plot_type]) > 1:
                    legend_text = f"{legend_text}, {varname.replace('accumulated', '').replace('Flux', '')}"
                legendText.append(legend_text)
                lineColors.append(config.get('timeSeries', 'controlColor'))
                lineStyles.append(lineStylesBase[index])

        lineWidths = [3 for i in fields]
        if config.has_option('timeSeries', 'movingAveragePoints'):
            movingAveragePoints = config.getint('timeSeries',
                                                'movingAveragePoints')
        else:
            movingAveragePoints = None

        if config.has_option('timeSeries', 'firstYearXTicks'):
            firstYearXTicks = config.getint('timeSeries',
                                            'firstYearXTicks')
        else:
            firstYearXTicks = None

        if config.has_option('timeSeries', 'yearStrideXTicks'):
            yearStrideXTicks = config.getint('timeSeries',
                                             'yearStrideXTicks')
        else:
            yearStrideXTicks = None

        timeseries_analysis_plot(config, fields, calendar=self.calendar,
                                 title=title, xlabel=xLabel, ylabel=yLabel,
                                 movingAveragePoints=movingAveragePoints,
                                 lineColors=lineColors,
                                 lineStyles=lineStyles[:len(fields)],
                                 lineWidths=lineWidths,
                                 legendText=legendText,
                                 firstYearXTicks=firstYearXTicks,
                                 yearStrideXTicks=yearStrideXTicks)

        # save the plot to the output file
        plt.savefig(outFileName)

        caption = captions[plot_type]
        write_image_xml(
            config=self.config,
            filePrefix=filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='Time Series',
            groupLink='timeseries',
            gallery='Conservation',
            thumbnailDescription=title,
            imageDescription=caption,
            imageCaption=caption)

    def _get_variable(self, ds, varname, mks=False):
        if varname not in self.derivedVariableList:
            variable = ds[varname]
        else:
            # Here we keep the units mks
            if varname == 'netMassChange':
                variable = self._get_variable(ds, 'massChange', mks=True)
                # mass_flux = self._get_variable(ds, 'netMassFlux')
                # # Assume that the frequency of output is monthly
                # dt = constants.sec_per_month
                # # Convert from kg/s to kg
                # derived_variable = mass_flux.cumsum(axis=0) * dt
            elif varname == 'landIceMassChange':
                land_ice_mass_flux = self._get_variable(ds, 'landIceMassFlux', mks=True)
                # Assume that the frequency of output is monthly
                dt = constants.sec_per_month
                # Convert from kg/s to kg/month
                land_ice_mass_flux = land_ice_mass_flux * dt
                # Convert from kg/month to kg
                variable = np.cumsum(land_ice_mass_flux)

            elif varname == 'landIceMassFlux':
                variable = self._get_variable(ds, 'accumulatedIcebergFlux', mks=True) + \
                           self._get_variable(ds, 'accumulatedLandIceFlux', mks=True) + \
                           self._get_variable(ds, 'accumulatedRemovedRiverRunoffFlux', mks=True) + \
                           self._get_variable(ds, 'accumulatedRemovedIceRunoffFlux', mks=True)

            elif varname == 'landIceSshChange':
                ts_files = sorted(self.historyStreams.readpath(
                    'timeSeriesStatsMonthlyOutput',
                    startDate=f'{self.startYear:04d}-01-01_00:00:00',
                    endDate=f'{self.endYear:04d}-01-01_00:00:00',
                    calendar=self.calendar))
                # Note that here we assume that the area of the ocean is constant in time
                # to save computational expense because most configurations do not allow
                # the area of the ocean to change
                ts_file = ts_files[0]
                if not os.path.exists(ts_file):
                   raise ValueError(f'Could not find timeMonthlyStats file {ts_file}')
                var = 'timeMonthly_avg_areaCellGlobal' 
                ds_ts = open_mpas_dataset(fileName=ts_file,
                                          calendar=self.calendar,
                                          variableList=[var])
                A = ds_ts[var].mean()
                land_ice_mass_change = self._get_variable(ds, 'landIceMassChange', mks=True)
                rho = self.namelist.getfloat('config_density0')
                # Convert from to kg to m
                variable = land_ice_mass_change / (rho * A)

            else:
                raise ValueError(f'Attempted to derive non-supported variable {varname}')

        removed_vars = ['accumulatedRemovedRiverRunoffFlux',
                        'accumulatedRemovedIceRunoffFlux']
        if varname in removed_vars:
            variable = -variable

        if not mks:
            # Here we do all the unit conversion from mks into whatever we want
            mass_vars = ['initialMass', 'finalMass', 'absoluteMassError',
                         'relativeMassError', 'massChange', 'landIceMassChange']
            salt_vars = ['initialSalt', 'finalSalt', 'absoluteSaltError',
                         'relativeSaltError']
            mass_flux_vars = ['netMassFlux', 'landIceMassFlux']
            salt_flux_vars = ['netSaltFlux']
            ssh_vars = ['landIceSshChange', 'sshChange']
            if (varname in mass_vars) or (varname in salt_vars):
                # Convert from kg to Gt
                variable = variable * 1e-12
            if (varname in mass_flux_vars) or (varname in salt_flux_vars):
                # Convert from kg/s to Gt/yr
                variable = variable * 1e-12 * constants.sec_per_year
            if varname in ssh_vars:
                # Convert from m to mm
                variable = variable * 1e3

        return variable

    def _compute_time_series_with_ncrcat(self, variable_list):

        """
        Uses ncrcat to extact time series from conservationCheckOutput files

        Raises
        ------
        OSError
            If ``ncrcat`` is not in the system path.
        """

        if find_executable('ncrcat') is None:
            raise OSError('ncrcat not found. Make sure the latest nco '
                          'package is installed: \n'
                          'conda install nco\n'
                          'Note: this presumes use of the conda-forge '
                          'channel.')

        inputFiles = self.inputFiles
        append = False
        if os.path.exists(self.outputFile):
            # make sure all the necessary variables are also present
            with xr.open_dataset(self.outputFile) as ds:
                if ds.sizes['Time'] == 0:
                    updateSubset = False
                else:
                    updateSubset = True
                    for variableName in variable_list:
                        if variableName not in ds.variables:
                            updateSubset = False
                            break

                if updateSubset:
                    # add only input files with times that aren't already in
                    # the output file

                    append = True

                    fileNames = sorted(self.inputFiles)
                    inYears, inMonths = get_files_year_month(
                        fileNames, self.historyStreams,
                        'conservationCheckOutput')

                    inYears = np.array(inYears)
                    inMonths = np.array(inMonths)
                    totalMonths = 12 * inYears + inMonths

                    dates = decode_strings(ds.xtime)

                    lastDate = dates[-1]

                    lastYear = int(lastDate[0:4])
                    lastMonth = int(lastDate[5:7])
                    lastTotalMonths = 12 * lastYear + lastMonth

                    inputFiles = []
                    for index, inputFile in enumerate(fileNames):
                        if totalMonths[index] > lastTotalMonths:
                            inputFiles.append(inputFile)

                    if len(inputFiles) == 0:
                        # nothing to do
                        return
                else:
                    # there is an output file but it has the wrong variables
                    # so we need ot delete it.
                    self.logger.warning('Warning: deleting file {self.outputFile}'
                                        ' because it is empty or some variables'
                                        ' were missing')
                    os.remove(self.outputFile)

        variableList = variable_list + ['xtime']

        args = ['ncrcat', '-4', '--no_tmp_fl',
                '-v', ','.join(variableList)]

        if append:
            args.append('--record_append')

        printCommand = '{} {} ... {} {}'.format(' '.join(args), inputFiles[0],
                                                inputFiles[-1],
                                                self.outputFile)
        args.extend(inputFiles)
        args.append(self.outputFile)

        self.logger.info(f'running: {printCommand}')
        for handler in self.logger.handlers:
            handler.flush()

        process = subprocess.Popen(args, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        if stdout:
            stdout = stdout.decode('utf-8')
            for line in stdout.split('\n'):
                self.logger.info(line)
        if stderr:
            stderr = stderr.decode('utf-8')
            for line in stderr.split('\n'):
                self.logger.error(line)

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode,
                                                ' '.join(args))
