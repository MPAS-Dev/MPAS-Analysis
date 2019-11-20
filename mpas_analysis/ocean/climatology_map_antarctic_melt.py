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
from pyremap import ProjectionGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.io.utility import build_obs_path

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask, get_antarctic_stereographic_projection

from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask

from mpas_analysis.shared.constants import constants


class ClimatologyMapAntarcticMelt(AnalysisTask):  # {{{
    """
    An analysis task for comparison of Antarctic melt rates against
    observations
    """
    # Authors
    # -------
    # Xylar Asay-Davis

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

        fieldName = 'meltRate'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapAntarcticMelt, self).__init__(
            config=config, taskName='climatologyMapAntarcticMelt',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', fieldName,
                  'landIceCavities', 'antarctic'])

        sectionName = self.taskName

        mpasFieldName = 'timeMonthly_avg_landIceFreshwaterFlux'
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

        # the variable 'timeMonthly_avg_landIceFreshwaterFlux' will be added to
        # mpasClimatologyTask along with the seasons.
        remapClimatologySubtask = RemapMpasAntarcticMeltClimatology(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=fieldName,
            variableList=[mpasFieldName],
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=iselValues)

        if controlConfig is None:

            refTitleLabel = \
                'Observations (Rignot et al, 2013)'

            observationsDirectory = build_obs_path(
                config, 'ocean', 'meltSubdirectory')

            obsFileName = \
                '{}/Rignot_2013_melt_rates_6000.0x6000.0km_10.0km_' \
                'Antarctic_stereo.nc'.format(observationsDirectory)
            refFieldName = 'meltRate'
            outFileLabel = 'meltRignot'
            galleryName = 'Observations: Rignot et al. (2013)'

            remapObservationsSubtask = RemapObservedAntarcticMeltClimatology(
                parentTask=self, seasons=seasons, fileName=obsFileName,
                outFilePrefix=refFieldName,
                comparisonGridNames=comparisonGridNames)
            self.add_subtask(remapObservationsSubtask)
            diffTitleLabel = 'Model - Observations'

        else:
            remapObservationsSubtask = None
            controlRunName = controlConfig.get('runs', 'mainRunName')
            galleryName = None
            refTitleLabel = 'Control: {}'.format(controlRunName)

            refFieldName = mpasFieldName
            outFileLabel = 'melt'
            diffTitleLabel = 'Main - Control'

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(
                    self, season, comparisonGridName, remapClimatologySubtask,
                    remapObservationsSubtask, controlConfig=controlConfig)

                subtask.set_plot_info(
                    outFileLabel=outFileLabel,
                    fieldNameInTitle='Melt Rate',
                    mpasFieldName=mpasFieldName,
                    refFieldName=refFieldName,
                    refTitleLabel=refTitleLabel,
                    diffTitleLabel=diffTitleLabel,
                    unitsLabel=r'm a$^{-1}$',
                    imageCaption='Antarctic Melt Rate',
                    galleryGroup='Melt Rate',
                    groupSubtitle=None,
                    groupLink='antarctic_melt',
                    galleryName=galleryName)

                self.add_subtask(subtask)
        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class
        # (AnalysisTask), which will perform some common setup
        super(ClimatologyMapAntarcticMelt, self).setup_and_check()

        landIceFluxMode = self.namelist.get('config_land_ice_flux_mode')
        if landIceFluxMode not in ['standalone', 'coupled']:
            raise ValueError('*** climatologyMapMeltAntarctic requires '
                             'config_land_ice_flux_mode \n'
                             '    to be standalone or coupled.  Otherwise, no '
                             'melt rates are available \n'
                             '    for plotting.')
        # }}}
    # }}}


class RemapMpasAntarcticMeltClimatology(RemapMpasClimatologySubtask):  # {{{
    """
    A subtask for remapping climatologies of Antarctic melt rates and adding

    Attributes
    ----------
    landIceMask : xarray.DataArray
        A mask indicating where there is land ice on the ocean grid (thus,
        where melt rates are valid)
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def run_task(self):  # {{{
        """
        Compute climatologies of melt rates from E3SM/MPAS output

        This function has been overridden to load ``landIceMask`` from a
        restart file for later use in masking the melt rate.  It then simply
        calls the run function from
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, load the land-ice mask from the restart file
        dsLandIceMask = xr.open_dataset(self.restartFileName)
        dsLandIceMask = dsLandIceMask[['landIceMask']]
        dsLandIceMask = dsLandIceMask.isel(Time=0)
        self.landIceMask = dsLandIceMask.landIceMask > 0.

        # then, call run from the base class (RemapMpasClimatologySubtask),
        # which will perform the main function of the task
        super(RemapMpasAntarcticMeltClimatology, self).run_task()

        # }}}

    def customize_masked_climatology(self, climatology, season):  # {{{
        """
        Mask the melt rates using ``landIceMask`` and rescale it to m/yr

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
        # Xylar Asay-Davis

        fieldName = self.variableList[0]

        # scale the field to m/yr from kg/m^2/s and mask out non-land-ice areas
        climatology[fieldName] = \
            constants.sec_per_year / constants.rho_fw * \
            climatology[fieldName].where(self.landIceMask)

        return climatology  # }}}

    # }}}


class RemapObservedAntarcticMeltClimatology(RemapObservedClimatologySubtask):
    # {{{
    """
    A subtask for reading and remapping Antarctic melt-rate observations
    """
    # Authors
    # -------
    # Xylar Asay-Davis

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

        # create a descriptor of the observation grid using the x/y polar
        # stereographic coordinates
        projection = get_antarctic_stereographic_projection()
        obsDescriptor = ProjectionGridDescriptor.read(
            projection, fileName=fileName, xVarName='x', yVarName='y')
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

        return dsObs  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
