# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2018 Los Alamos National Security, LLC. All rights reserved.
# Copyright (c) 2018 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2018 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
'''
Analysis tasks for comparing climatology maps of Antarctic seafloor fields
against observations from Schmidtko et al. (2014, DOI: 10.1126/science.1256117)
'''
# Authors
# -------
# Xylar Asay-Davis

import xarray as xr

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.ocean.remap_depth_slices_subtask import \
    RemapDepthSlicesSubtask
from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask

from mpas_analysis.shared.io.utility import build_config_full_path

from mpas_analysis.shared.climatology import RemapObservedClimatologySubtask, \
    get_antarctic_stereographic_projection

from mpas_analysis.shared.grid import ProjectionGridDescriptor


class ClimatologyMapSchmidtko(AnalysisTask):  # {{{
    """
    An analysis task for comparison of Antarctic Bottom Water properties
    against Schmidtko et al. (2014) observations
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, controlConfig=None):  # {{{
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

        fields = \
            [{'prefix': 'temperature',
              'mpas': 'timeMonthly_avg_activeTracers_temperature',
              'units': r'$\degree$C',
              'obs': 'botTheta',
              'title': 'Potential Temperature'},
             {'prefix': 'salinity',
              'mpas': 'timeMonthly_avg_activeTracers_salinity',
              'units': r'PSU',
              'obs': 'botSalinity',
              'title': 'Salinity'},
             {'prefix': 'potentialDensity',
              'mpas': 'timeMonthly_avg_potentialDensity',
              'units': r'kg m$^{-3}$',
              'obs': 'botPotentialDensity',
              'title': 'Potential Density'}]

        tags = ['climatology', 'horizontalMap'] \
            + [field['prefix'] for field in fields]

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSchmidtko, self).__init__(
                config=config, taskName='climatologyMapSchmidtko',
                componentName='ocean', tags=tags)

        sectionName = self.taskName

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

        observationsDirectory = build_config_full_path(
            config, 'oceanObservations', 'schmidtkoSubdirectory')

        obsFileName = '{}/Schmidtko_et_al_2014_bottom_PT_S_PD_' \
                      '6000.0x6000.0km_10.0km_Antarctic_stereo.nc' \
                      ''.format(observationsDirectory)

        # the variable 'timeMonthly_avg_landIceFreshwaterFlux' will be added to
        # mpasClimatologyTask along with the seasons.
        remapClimatologySubtask = RemapDepthSlicesSubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='schmidtko',
            variableList=[field['mpas'] for field in fields],
            seasons=seasons,
            depths=['bot'],
            comparisonGridNames=comparisonGridNames)

        if controlConfig is None:
            refTitleLabel = 'Observations: Schmidtko et al. (2014)'
            diffTitleLabel = 'Model - Observations'
            groupSubtitle = refTitleLabel
        else:
            controlRunName = controlConfig.get('runs', 'mainRunName')
            refTitleLabel = 'Control: {}'.format(controlRunName)
            diffTitleLabel = 'Main - Control'
            groupSubtitle = None

        for field in fields:
            fieldPrefix = field['prefix']
            upperFieldPrefix = fieldPrefix[0].upper() + fieldPrefix[1:]
            if controlConfig is None:
                refFieldName = field['obs']
                outFileLabel = '{}Schmidtko'.format(fieldPrefix)

                remapObservationsSubtask = RemapSchmidtko(
                        parentTask=self, seasons=seasons, fileName=obsFileName,
                        outFilePrefix='{}Schmidtko'.format(fieldPrefix),
                        fieldName=refFieldName,
                        comparisonGridNames=comparisonGridNames,
                        subtaskName='remapObservations{}'.format(
                                upperFieldPrefix))

                self.add_subtask(remapObservationsSubtask)

            else:
                remapObservationsSubtask = None
                refFieldName = field['mpas']
                outFileLabel = '{}Bottom'.format(fieldPrefix)
                diffTitleLabel = 'Main - Control'

            for comparisonGridName in comparisonGridNames:
                for season in seasons:
                    subtask = PlotClimatologyMapSubtask(
                        parentTask=self,
                        season=season,
                        comparisonGridName=comparisonGridName,
                        remapMpasClimatologySubtask=remapClimatologySubtask,
                        remapObsClimatologySubtask=remapObservationsSubtask,
                        controlConfig=controlConfig,
                        depth='bot',
                        subtaskName='plot{}_{}_{}'.format(upperFieldPrefix,
                                                          season,
                                                          comparisonGridName))

                    subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle=field['title'],
                        mpasFieldName=field['mpas'],
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=field['units'],
                        imageCaption=field['title'],
                        galleryGroup='Seafloor Maps',
                        groupSubtitle=groupSubtitle,
                        groupLink='seafloor',
                        galleryName=field['title'],
                        configSectionName='climatologyMapSchmidtko{}'.format(
                            upperFieldPrefix))

                    self.add_subtask(subtask)
        # }}}

    # }}}


class RemapSchmidtko(RemapObservedClimatologySubtask):  # {{{
    """
    A subtask for reading and remapping Schmidtko et al. (2014) observations
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, seasons, fileName, outFilePrefix,
                 fieldName, comparisonGridNames=['latlon'],
                 subtaskName='remapObservations'):
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

        fieldName : str
            The name of the 3D field to remap

        comparisonGridNames : list of {'latlon', 'antarctic'}, optional
            The name(s) of the comparison grid to use for remapping.

        subtaskName : str, optional
            The name of the subtask
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.fieldName = fieldName

        # call the constructor from the base class
        # (RemapObservedClimatologySubtask)
        super(RemapSchmidtko, self).__init__(
                parentTask, seasons, fileName, outFilePrefix,
                comparisonGridNames, subtaskName)
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

        field = xr.concat([dsObs[self.fieldName]], dim='depthSlice')

        dsObs = xr.Dataset(data_vars={self.fieldName: field},
                           coords={'depthSlice': ['bot']})

        return dsObs  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
