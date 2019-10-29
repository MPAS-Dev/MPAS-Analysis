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
from pyremap import LatLonGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.io.utility import build_obs_path

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask

from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask


class ClimatologyMapEKE(AnalysisTask):  # {{{
    """
    An analysis task for comparison of eddy kinetic energy (eke) against
    observations
    """
    # Authors
    # -------
    # Kevin Rosa

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
        # Xylar Asay-Davis, Kevin Rosa

        fieldName = 'eke'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapEKE, self).__init__(
            config=config, taskName='climatologyMapEKE',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', fieldName, 'publicObs'])

        mpasFieldName = 'eke'
        iselValues = {'nVertLevels': 0}

        sectionName = self.taskName

        # read in what seasons we want to plot
        seasons = config.getExpression(sectionName, 'seasons')

        # EKE observations are annual climatology so only accept annual
        # climatology **should move this to setup_and_check()
        if seasons != ['ANN']:
            raise ValueError('config section {} does not contain valid list '
                             'of seasons. For EKE, may only request annual '
                             'climatology'.format(sectionName))

        comparisonGridNames = config.getExpression(sectionName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of comparison grids'.format(sectionName))

        # the variables in variableList will be added to mpasClimatologyTask
        # along with the seasons.
        variableList = ['timeMonthly_avg_velocityZonal',
                        'timeMonthly_avg_velocityMeridional',
                        'timeMonthly_avg_velocityZonalSquared',
                        'timeMonthly_avg_velocityMeridionalSquared']
        remapClimatologySubtask = RemapMpasEKEClimatology(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=fieldName,
            variableList=variableList,
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=iselValues)

        # to compare to observations:
        if controlConfig is None:

            refTitleLabel = \
                'Observations (Surface EKE from Drifter Data)'

            observationsDirectory = build_obs_path(
                config, 'ocean', '{}Subdirectory'.format(fieldName))

            obsFileName = \
                "{}/drifter_variance_20180804.nc".format(
                    observationsDirectory)
            refFieldName = 'eke'
            outFileLabel = 'ekeDRIFTER'
            galleryName = 'Observations: EKE from Drifters'

            remapObservationsSubtask = RemapObservedEKEClimatology(
                parentTask=self, seasons=seasons, fileName=obsFileName,
                outFilePrefix=refFieldName,
                comparisonGridNames=comparisonGridNames)
            self.add_subtask(remapObservationsSubtask)
            diffTitleLabel = 'Model - Observations'

        # compare with previous run:
        else:
            remapObservationsSubtask = None
            controlRunName = controlConfig.get('runs', 'mainRunName')
            galleryName = None
            refTitleLabel = 'Control: {}'.format(controlRunName)

            refFieldName = mpasFieldName
            outFileLabel = 'eke'
            diffTitleLabel = 'Main - Control'

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(
                    self, season, comparisonGridName, remapClimatologySubtask,
                    remapObservationsSubtask, controlConfig=controlConfig)

                subtask.set_plot_info(
                    outFileLabel=outFileLabel,
                    fieldNameInTitle='EKE',
                    mpasFieldName=mpasFieldName,
                    refFieldName=refFieldName,
                    refTitleLabel=refTitleLabel,
                    diffTitleLabel=diffTitleLabel,
                    unitsLabel=r'cm$^2$/s$^2$',
                    imageCaption='Mean Surface Eddy Kinetic Energy',
                    galleryGroup='Eddy Kinetic Energy',
                    groupSubtitle=None,
                    groupLink='eke',
                    galleryName=galleryName)

                self.add_subtask(subtask)
        # }}}
    # }}}


# adds to the functionality of RemapDepthSlicesSubtask
class RemapMpasEKEClimatology(RemapMpasClimatologySubtask):  # {{{
    """
    A subtask for computing climatologies of eddy kinetic energy from means of
    velocity and velocity-squared.
    """
    # Authors
    # -------
    # Kevin Rosa

    def customize_masked_climatology(self, climatology, season):  # {{{
        """
        Construct velocity magnitude as part of the climatology

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

        # first, call the base class's version of this function so we extract
        # the desired slices.
        climatology = super(RemapMpasEKEClimatology,
                            self).customize_masked_climatology(climatology,
                                                               season)

        # calculate mpas eddy kinetic energy
        scaleFactor = 100 * 100  # m2/s2 to cm2/s2
        eke = 0.5 * scaleFactor * \
            (climatology.timeMonthly_avg_velocityZonalSquared -
             climatology.timeMonthly_avg_velocityZonal ** 2 +
             climatology.timeMonthly_avg_velocityMeridionalSquared -
             climatology.timeMonthly_avg_velocityMeridional ** 2)

        # drop unnecessary fields before re-mapping
        climatology.drop(['timeMonthly_avg_velocityZonal',
                          'timeMonthly_avg_velocityMeridional',
                          'timeMonthly_avg_velocityZonalSquared',
                          'timeMonthly_avg_velocityMeridionalSquared'])

        # this creates a new variable eke in climatology (like netcdf)
        climatology['eke'] = eke
        climatology.eke.attrs['units'] = 'cm$^[2]$ s$^{-2}$'
        climatology.eke.attrs['description'] = 'eddy kinetic energy'

        return climatology  # }}}

    # }}}


class RemapObservedEKEClimatology(RemapObservedClimatologySubtask):  # {{{
    """
    A subtask for reading and remapping EKE observations
    """
    # Authors
    # -------
    # Kevin Rosa

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
        # Kevin Rosa

        # create a descriptor of the observation grid using the lat/lon
        # coordinates
        obsDescriptor = LatLonGridDescriptor.read(fileName=fileName,
                                                  latVarName='Lat',
                                                  lonVarName='Lon')

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
        # Kevin Rosa, Xylar Asay-Davis

        dsIn = xr.open_dataset(fileName)

        scaleFactor = 100 * 100  # m2/s2 to cm2/s2
        eke = 0.5 * scaleFactor * \
            (dsIn['Up2bar'].values + dsIn['Vp2bar'].values)

        # create a new dataset for the observations. solves transpose issues.
        dsObs = xr.Dataset({'eke': (['latitude', 'longitude'], eke.T)},
                           coords={'Lat': (['latitude'], dsIn.Lat.values),
                                   'Lon': (['longitude'], dsIn.Lon.values)}
                           )
        # update attributes
        dsObs.eke.attrs['units'] = 'cm$^2$ s$^{-2}$'
        dsObs.eke.attrs['long_name'] = 'Eddy kinetic energy'

        return dsObs  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
