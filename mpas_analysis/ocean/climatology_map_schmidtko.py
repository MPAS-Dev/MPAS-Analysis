# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#
'''
Analysis tasks for comparing climatology maps of Antarctic seafloor fields
against observations from Schmidtko et al. (2014, DOI: 10.1126/science.1256117)
'''
# Authors
# -------
# Xylar Asay-Davis

import xarray as xr
import numpy as np

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.ocean.remap_depth_slices_subtask import \
    RemapDepthSlicesSubtask
from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask

from mpas_analysis.shared.io.utility import build_config_full_path

from mpas_analysis.shared.climatology import RemapObservedClimatologySubtask, \
    get_antarctic_stereographic_projection

from mpas_analysis.shared.grid import ProjectionGridDescriptor

from mpas_analysis.shared.mpas_xarray import mpas_xarray


class ClimatologyMapSchmidtko(AnalysisTask):  # {{{
    """
    An analysis task for comparison of Antarctic Bottom Water properties
    against Schmidtko et al. (2014) observations
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, refConfig=None):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        fields = \
            {'temperature':
                {'mpas': 'timeMonthly_avg_activeTracers_temperature',
                 'units': r'$^\circ$C',
                 'obs': 'botTheta',
                 'title': 'Temperature'},
             'salinity':
                {'mpas': 'timeMonthly_avg_activeTracers_salinity',
                 'units': r'PSU',
                 'obs': 'botSalinity',
                 'title': 'Salinity'},
             'potentialDensity':
                {'mpas': 'timeMonthly_avg_potentialDensity',
                 'units': r'kg m$^{-3}$',
                 'obs': 'botPotentialDensity',
                 'title': 'Potential Density'}}

        tags = ['climatology', 'horizontalMap'] + list(fields.keys())
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
            variableList=[field['mpas'] for field in fields.values()],
            seasons=seasons,
            depths=['bot'],
            comparisonGridNames=comparisonGridNames)

        if refConfig is None:
            refTitleLabel = 'Observations: Schmidtko et al. (2014)'
            diffTitleLabel = 'Model - Observations'
            groupSubtitle = refTitleLabel
        else:
            refRunName = refConfig.get('runs', 'mainRunName')
            refTitleLabel = 'Ref: {}'.format(refRunName)
            diffTitleLabel = 'Main - Reference'
            groupSubtitle = None

        for fieldName in fields:
            upperFieldName = fieldName[0].upper() + fieldName[1:]
            if refConfig is None:
                refFieldName = fields[fieldName]['obs']
                outFileLabel = '{}Schmidtko'.format(fieldName)

                remapObservationsSubtask = RemapSchmidtko(
                        parentTask=self, seasons=seasons, fileName=obsFileName,
                        outFilePrefix='{}Schmidtko'.format(fieldName),
                        fieldName=refFieldName,
                        comparisonGridNames=comparisonGridNames,
                        subtaskName='remapObservations{}'.format(
                                upperFieldName))

                self.add_subtask(remapObservationsSubtask)

            else:
                remapObservationsSubtask = None
                refRunName = refConfig.get('runs', 'mainRunName')

                refFieldName = fields[fieldName]['mpas']
                outFileLabel = '{}Bottom'.format(fieldName)
                diffTitleLabel = 'Main - Reference'

            for comparisonGridName in comparisonGridNames:
                for season in seasons:
                    subtask = PlotClimatologyMapSubtask(
                        parentTask=self,
                        season=season,
                        comparisonGridName=comparisonGridName,
                        remapMpasClimatologySubtask=remapClimatologySubtask,
                        remapObsClimatologySubtask=remapObservationsSubtask,
                        refConfig=refConfig,
                        depth='bot',
                        subtaskName='plot{}_{}_{}'.format(upperFieldName,
                                                          season,
                                                          comparisonGridName))

                    subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle=fields[fieldName]['title'],
                        mpasFieldName=fields[fieldName]['mpas'],
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=fields[fieldName]['units'],
                        imageCaption=fields[fieldName]['title'],
                        galleryGroup='Seafloor Maps',
                        groupSubtitle=groupSubtitle,
                        groupLink='seafloor',
                        galleryName=fields[fieldName]['title'],
                        configSectionName='climatologyMapSchmidtko{}'.format(
                            upperFieldName))

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
