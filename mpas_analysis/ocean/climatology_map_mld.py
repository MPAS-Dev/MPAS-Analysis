# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr
import numpy as np
from pyremap import LatLonGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.io.utility import build_obs_path

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask

from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask

from mpas_analysis.shared.mpas_xarray import mpas_xarray


class ClimatologyMapMLD(AnalysisTask):  # {{{
    """
    An analysis task for comparison of mixed layer depth (mld) against
    observations
    """
    # Authors
    # -------
    # Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

    def __init__(self, config, mpasClimatologyTask,
                 controlConfig=None):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        fieldName = 'mld'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapMLD, self).__init__(
            config=config, taskName='climatologyMapMLD',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', fieldName, 'publicObs'])

        sectionName = self.taskName

        mpasFieldName = 'timeMonthly_avg_dThreshMLD'
        iselValues = None

        # read in what seasons we want to plot
        seasons = config.getExpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of seasons'.format(sectionName))

        comparisonGridNames = config.getExpression(sectionName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of comparison grids'.format(sectionName))

        # the variable 'timeMonthly_avg_dThreshMLD' will be added to
        # mpasClimatologyTask along with the seasons.
        remapClimatologySubtask = RemapMpasClimatologySubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=fieldName,
            variableList=[mpasFieldName],
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=iselValues)

        if controlConfig is None:

            observationsDirectory = build_obs_path(
                config, 'ocean', '{}Subdirectory'.format(fieldName))

            obsFileName = "{}/holtetalley_mld_climatology_20180710.nc".format(
                observationsDirectory)

            refFieldName = 'mld'
            outFileLabel = 'mldHolteTalleyARGO'

            remapObservationsSubtask = RemapObservedMLDClimatology(
                parentTask=self, seasons=seasons, fileName=obsFileName,
                outFilePrefix=refFieldName,
                comparisonGridNames=comparisonGridNames)
            self.add_subtask(remapObservationsSubtask)
            galleryName = 'Observations: Holte-Talley ARGO'
            refTitleLabel = \
                'Observations (HolteTalley density threshold MLD)'
            diffTitleLabel = 'Model - Observations'

        else:
            remapObservationsSubtask = None
            controlRunName = controlConfig.get('runs', 'mainRunName')
            galleryName = None
            refTitleLabel = 'Control: {}'.format(controlRunName)

            refFieldName = mpasFieldName
            outFileLabel = 'mld'
            diffTitleLabel = 'Main - Control'

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(
                    self, season, comparisonGridName, remapClimatologySubtask,
                    remapObservationsSubtask, controlConfig=controlConfig)

                subtask.set_plot_info(
                    outFileLabel=outFileLabel,
                    fieldNameInTitle='MLD',
                    mpasFieldName=mpasFieldName,
                    refFieldName=refFieldName,
                    refTitleLabel=refTitleLabel,
                    diffTitleLabel=diffTitleLabel,
                    unitsLabel=r'm',
                    imageCaption='Mean Mixed-Layer Depth',
                    galleryGroup='Mixed-Layer Depth',
                    groupSubtitle=None,
                    groupLink='mld',
                    galleryName=galleryName)

                self.add_subtask(subtask)
        # }}}

    def setup_and_check(self):  # {{{
        '''
        Check if MLD capability was turned on in the run.
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ClimatologyMapMLD, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_mixedlayerdepths_enable',
            raiseException=True)

        # }}}

    # }}}


class RemapObservedMLDClimatology(RemapObservedClimatologySubtask):  # {{{
    """
    A subtask for reading and remapping MLD observations
    """
    # Authors
    # -------
    # Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

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
        # Xylar Asay-Davis

        # Load MLD observational data
        dsObs = self.build_observational_dataset(fileName)

        # create a descriptor of the observation grid using the lat/lon
        # coordinates
        obsDescriptor = LatLonGridDescriptor.read(ds=dsObs,
                                                  latVarName='lat',
                                                  lonVarName='lon')
        dsObs.close()
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
        # Xylar Asay-Davis

        # Load MLD observational data
        dsObs = xr.open_dataset(fileName)

        # Increment month value to be consistent with the model output
        dsObs.iMONTH.values += 1
        # Rename the dimensions to be consistent with other obs. data sets
        dsObs = dsObs.rename({'month': 'calmonth', 'lat': 'latCoord',
                              'lon': 'lonCoord', 'mld_dt_mean': 'mld'})
        dsObs = dsObs.rename({'iMONTH': 'Time', 'iLAT': 'lat', 'iLON': 'lon'})

        # set the coordinates now that the dimensions have the same names
        dsObs.coords['lat'] = dsObs['latCoord']
        dsObs.coords['lon'] = dsObs['lonCoord']
        dsObs.coords['Time'] = dsObs['calmonth']
        dsObs.coords['month'] = ('Time', np.array(dsObs['calmonth'], int))

        # no meaningful year since this is already a climatology
        dsObs.coords['year'] = ('Time', np.ones(dsObs.dims['Time'], int))

        dsObs = mpas_xarray.subset_variables(dsObs, ['mld', 'month'])
        return dsObs  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
