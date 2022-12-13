# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2020 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2020 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2020 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr
from pyremap import LatLonGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.io.utility import build_obs_path

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask

from mpas_analysis.shared.plot import PlotClimatologyMapSubtask
from mpas_analysis.shared.interpolation.utility import add_periodic_lon

import datetime


class ClimatologyMapWaves(AnalysisTask):  # {{{
    """
    An analysis task for the computation of wave climatologies
    """
    # Authors
    # -------
    # Steven Brus
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, regionMasksTask,
                 controlConfig):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        regionMasksTask : ``ComputeRegionMasks``
            A task for computing region masks

        controlConfig :  ``MpasAnalysisConfigParser``
            Configuration options for a control run
        """
        # Authors
        # -------
        # Steven Brus
        # Xylar Asay-Davis

        fields = [{'prefix': 'significantWaveHeight',
                   'mpas': 'timeMonthly_avg_significantWaveHeight',
                   'units': r'm',
                   'titleName': 'Significant Wave Height',
                   'era5': 'swh',
                   'ss_cci': 'swh_mean'},
                  {'prefix': 'peakWavePeriod',
                   'mpas': 'timeMonthly_avg_peakWaveFrequency',
                   'units': r's',
                   'titleName': 'Peak Wave Period',
                   'era5': 'pp1d'},
                  {'prefix': 'iceFraction',
                   'mpas': 'timeMonthly_avg_iceFraction',
                   'units': r'',
                   'titleName': 'Ice Fraction'},
                  ]

        # call the constructor from the base class (AnalysisTask)
        super().__init__(
            config=config, taskName='climatologyMapWaves',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', 'waves'])

        sectionName = self.taskName

        # read in what seasons we want to plot
        seasons = config.getexpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError(f'config section {sectionName} '
                             'does not contain valid list of seasons')

        comparisonGridNames = config.getexpression(sectionName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError(f'config section {sectionName} does '
                             'not contain valid list of comparison grids')

        fieldList = config.getexpression(sectionName, 'fieldList')
        fieldsRequested = [field for field in fields
                           if field['prefix'] in fieldList]

        variableListMpas = [field['mpas'] for field in fieldsRequested]
        variableListEraObs = [field['era5'] for field in fieldsRequested]
        variableListSscciObs = [field['ss_cci'] for field in fieldsRequested
                                if 'ss_cci' in field]

        era5ObsStartYear = config.getint(sectionName, 'era5ObsStartYear')
        era5ObsEndYear = config.getint(sectionName, 'era5ObsEndYear')
        sscciObsStartYear = config.getint(sectionName, 'sscciObsStartYear')
        sscciObsEndYear = config.getint(sectionName, 'sscciObsEndYear')

        # the variableList variables  will be added to
        # mpasClimatologyTask along with the seasons.
        remapClimatologySubtask = RemapMpasWavesClimatology(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='wave',
            variableList=variableListMpas,
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=None)

        if controlConfig is None:
            observationsDirectory = build_obs_path(
                config, 'ocean', 'era5Subdirectory')

            obsFileName = f'{observationsDirectory}/ERA5_Monthly_1959-2021.nc'
            remapEraObservationsSubtask = RemapEra5ObservedWaveClimatology(
                parentTask=self, seasons=seasons, fileName=obsFileName,
                outFilePrefix='waveERA5',
                comparisonGridNames=comparisonGridNames,
                variableList=variableListEraObs,
                subtaskName='remapEraObservations')

            observationsDirectory = build_obs_path(
                config, 'ocean', 'sscciSubdirectory')

            obsFileName = f'{observationsDirectory}/SS_CCI_1991_2018.nc'
            remapSscciObservationsSubtask = RemapSscciObservedWaveClimatology(
                parentTask=self, seasons=seasons, fileName=obsFileName,
                outFilePrefix='waveSSCCI',
                comparisonGridNames=comparisonGridNames,
                variableList=variableListSscciObs,
                subtaskName='remapSscciObservations')

        else:
            remapObservationsSubtask = None

        for field in fieldsRequested:

            fieldPrefix = field['prefix']
            upperFieldPrefix = fieldPrefix[0].upper() + fieldPrefix[1:]
            configSectionName = f'{self.taskName}{upperFieldPrefix}'

            if (field['prefix'] == 'significantWaveHeight') \
                    and (controlConfig is None):
                obs_names = ['era5', 'ss_cci']
            else:
                obs_names = ['era5']

            for obs in obs_names:

                if controlConfig is None:

                    refFieldName = field[obs]
                    diffTitleLabel = None
                    if obs == 'era5':
                        obs_type = 'Reanalysis'
                        remapObservationsSubtask = \
                            remapEraObservationsSubtask
                        refTitleLabel = f'{obs_type} ({obs.upper()})' \
                                        f'{era5ObsStartYear}-{era5ObsEndYear}'
                    elif obs == 'ss_cci':
                        obs_type = 'Observations'
                        remapObservationsSubtask = \
                            remapSscciObservationsSubtask
                        refTitleLabel = f'{obs_type} ({obs.upper()})' \
                                        f'{sscciObsStartYear}-{sscciObsEndYear}'

                    galleryName = f"{field['titleName']} " \
                                  f"({obs.upper()} {obs_type})"

                else:
                    controlRunName = controlConfig.get('runs', 'mainRunName')
                    refTitleLabel = f'{field["titleName"]}' \
                                    f'(Control: {controlRunName})'
                    diffTitleLabel = 'Main - Control'
                    refFieldName = field['prefix']
                    galleryName = field['titleName']

                outFileLabel = f'{fieldPrefix}Wave{obs}'

                for comparisonGridName in comparisonGridNames:
                    for season in seasons:

                        subtaskName = f'plot{upperFieldPrefix}_{season}_' \
                                      f'{comparisonGridName}_{obs}'

                        # make a new subtask for this
                        # season and comparison grid
                        subtask = PlotClimatologyMapSubtask(
                            parentTask=self,
                            season=season,
                            comparisonGridName=comparisonGridName,
                            remapMpasClimatologySubtask=remapClimatologySubtask,
                            remapObsClimatologySubtask=remapObservationsSubtask,
                            controlConfig=controlConfig,
                            subtaskName=subtaskName)

                        subtask.set_plot_info(
                            outFileLabel=outFileLabel,
                            fieldNameInTitle=field['titleName'],
                            mpasFieldName=field['prefix'],
                            refFieldName=refFieldName,
                            refTitleLabel=refTitleLabel,
                            diffTitleLabel=diffTitleLabel,
                            unitsLabel=field['units'],
                            imageCaption=field['titleName'],
                            galleryGroup='Waves',
                            groupSubtitle=None,
                            groupLink='waves',
                            galleryName=galleryName,
                            configSectionName=configSectionName)

                        self.add_subtask(subtask)
        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.
        """
        # Authors
        # -------
        # Steven Brus
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class
        # (AnalysisTask), which will perform some common setup
        super().setup_and_check()

        if not self.namelist.getbool('config_use_active_wave'):

            raise ValueError('*** climatologyMapWaves requires '
                             'config_use_active_wave\n'
                             '    to be .true.  Otherwise, no '
                             'wave data is available \n'
                             '    for plotting.')
        # }}}
    # }}}


class RemapMpasWavesClimatology(RemapMpasClimatologySubtask):  # {{{
    """
    A subtask for reading and remapping wave climatology
    """
    # Authors
    # -------
    # Steven Brus, Xylar Asay-Davis

    def customize_masked_climatology(self, climatology, season):  # {{{
        """
        Convert peak freqency to peak period

        The descriptor for the observation grid
        Parameters
        ----------
        climatology : ``xarray.Dataset`` object
        the climatology data set

        season : str
        The name of the season to be masked

        Returns
        -------
        climatology : ``xarray.Dataset`` object
        the modified climatology data set
        """
        # Authors
        # -------
        # Steven Brus
        # Xylar Asay-Davis

        if 'timeMonthly_avg_peakWaveFrequency' in climatology:
            climatology['peakWavePeriod'] = \
                1.0/climatology['timeMonthly_avg_peakWaveFrequency']

        if 'timeMonthly_avg_significantWaveHeight' in climatology:
            climatology['significantWaveHeight'] = \
                climatology['timeMonthly_avg_significantWaveHeight']

        return climatology

        # }}}
    # }}}


class RemapEra5ObservedWaveClimatology(RemapObservedClimatologySubtask):  # {{{
    """
    A subtask for reading and remapping ERA5 wave observations
    """
    # Authors
    # -------
    # Steven Brus, Xylar Asay-Davis

    def __init__(self, parentTask, seasons, fileName, outFilePrefix,
                 variableList, subtaskName,
                 comparisonGridNames=['latlon']):
        # {{{
        '''
        Construct one analysis subtask for each plot (i.e. each season and
        comparison grid) and a subtask for computing climatologies.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent (master) task for this subtask

        seasons : list of str
           A list of seasons (keys in ``constants.monthDictionary``) over
           which the climatology should be computed.

        fileName : str
            The name of the observation file

        outFilePrefix : str
            The prefix in front of output files and mapping files, typically
            the name of the field being remapped

        variableList : list
            List of observational variables to remap

        comparisonGridNames : list of str,
            optional
            The name(s) of the comparison grid to use for remapping.

        subtaskName : str, optional
            The name of the subtask
        '''
        # Authors
        # -------
        # Steven Brus

        self.variableList = variableList

        # call the constructor from the base class (AnalysisTask)
        super().__init__(parentTask, seasons, fileName, outFilePrefix,
                         comparisonGridNames=comparisonGridNames,
                         subtaskName=subtaskName)
        # }}}

    def get_observation_descriptor(self, fileName):  # {{{
        '''
        get a MeshDescriptor for the observation grid

        Parameters
        ----------
        fileName : str
            observation file name describing the source grid

        Returns
        -------
        obsDescriptor : ``MeshDescriptor``
            The descriptor for the observation grid
        '''
        # Authors
        # -------
        # Steven Brus, Xylar Asay-Davis

        # create a descriptor of the observation grid using the lat/lon
        # coordinates
        dsObs = self.build_observational_dataset(fileName)
        obsDescriptor = LatLonGridDescriptor.read(ds=dsObs,
                                                  latVarName='latitude',
                                                  lonVarName='longitude')
        return obsDescriptor  # }}}

    def build_observational_dataset(self, fileName):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Parameters
        ----------
        fileName : str
            observation file name

        Returns
        -------
        dsObs : ``xarray.Dataset``
            The observational dataset
        '''
        # Authors
        # -------
        # Steven Brus, Xylar Asay-Davis

        sectionName = self.taskName
        climStartYear = self.config.getint(sectionName, 'era5ObsStartYear')
        climEndYear = self.config.getint(sectionName, 'era5ObsEndYear')
        timeStart = datetime.datetime(year=climStartYear, month=1, day=1)
        timeEnd = datetime.datetime(year=climEndYear, month=12, day=31)

        dsObs = xr.open_dataset(fileName)
        dsObs = dsObs.rename({'time': 'Time'})
        dsObs = dsObs.sel(Time=slice(timeStart, timeEnd))
        dsObs.coords['month'] = dsObs['Time.month']
        dsObs.coords['year'] = dsObs['Time.year']

        dsObs = dsObs[self.variableList]

        degrees = 'degree' in dsObs.longitude.units
        dsObs = add_periodic_lon(ds=dsObs, lonDim='longitude', degrees=degrees)

        return dsObs  # }}}

    # }}}


class RemapSscciObservedWaveClimatology(RemapObservedClimatologySubtask):  # {{{
    """
    A subtask for reading and remapping SS_CCI wave observations
    """
    # Authors
    # -------
    # Steven Brus, Xylar Asay-Davis

    def __init__(self, parentTask, seasons, fileName, outFilePrefix,
                 variableList, subtaskName,
                 comparisonGridNames=['latlon']):
        # {{{
        '''
        Construct one analysis subtask for each plot (i.e. each season and
        comparison grid) and a subtask for computing climatologies.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent (master) task for this subtask

        seasons : list of str
           A list of seasons (keys in ``constants.monthDictionary``) over
           which the climatology should be computed.

        fileName : str
            The name of the observation file

        outFilePrefix : str
            The prefix in front of output files and mapping files, typically
            the name of the field being remapped

        variableList : list
            List of observational variables to remap

        comparisonGridNames : list of str,
            optional
            The name(s) of the comparison grid to use for remapping.

        subtaskName : str, optional
            The name of the subtask
        '''
        # Authors
        # -------
        # Steven Brus, Xylar Asay-Davis

        self.variableList = variableList

        # call the constructor from the base class (AnalysisTask)
        super().__init__(parentTask, seasons, fileName, outFilePrefix,
                         comparisonGridNames=comparisonGridNames,
                         subtaskName=subtaskName)
        # }}}

    def get_observation_descriptor(self, fileName):  # {{{
        '''
        get a MeshDescriptor for the observation grid

        Parameters
        ----------
        fileName : str
            observation file name describing the source grid

        Returns
        -------
        obsDescriptor : ``MeshDescriptor``
            The descriptor for the observation grid
        '''
        # Authors
        # -------
        # Steven Brus, Xylar Asay-Davis

        # create a descriptor of the observation grid using the lat/lon
        # coordinates
        dsObs = self.build_observational_dataset(fileName)
        obsDescriptor = LatLonGridDescriptor.read(ds=dsObs,
                                                  latVarName='lat',
                                                  lonVarName='lon')
        return obsDescriptor  # }}}

    def build_observational_dataset(self, fileName):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Parameters
        ----------
        fileName : str
            observation file name

        Returns
        -------
        dsObs : ``xarray.Dataset``
            The observational dataset
        '''
        # Authors
        # -------
        # Steven Brus, Xylar Asay-Davis

        sectionName = self.taskName
        climStartYear = self.config.getint(sectionName, 'sscciObsStartYear')
        climEndYear = self.config.getint(sectionName, 'sscciObsEndYear')
        timeStart = datetime.datetime(year=climStartYear, month=1, day=1)
        timeEnd = datetime.datetime(year=climEndYear, month=12, day=31)

        dsObs = xr.open_dataset(fileName)
        dsObs = dsObs.rename({'time': 'Time'})
        dsObs = dsObs.sel(Time=slice(timeStart, timeEnd))
        dsObs.coords['month'] = dsObs['Time.month']
        dsObs.coords['year'] = dsObs['Time.year']

        dsObs = dsObs[self.variableList]

        degrees = 'degree' in dsObs.lon.units
        dsObs = add_periodic_lon(ds=dsObs, lonDim='lon', degrees=degrees)

        return dsObs  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
