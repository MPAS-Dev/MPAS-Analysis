# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.io.utility import build_config_full_path

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask

from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask

from mpas_analysis.shared.grid import LatLonGridDescriptor


class ClimatologyMapBGC(AnalysisTask):  # {{{
    """
    An analysis task for plotting of BGC variables

    Authors
    -------
    Phillip J. Wolfram, Riley X. Brady

    """
    def __init__(self, config, mpasClimatologyTask, afieldName, ampasFieldName,
                 refConfig=None):  # {{{
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

        Authors
        -------
        Phillip J. Wolfram, Riley X. Brady
        """
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapBGC, self).__init__(
            config=config, taskName='climatologyMapBGC_' + afieldName,
            componentName='ocean',
            tags=['climatology', 'horizontalMap', 'BGC', afieldName])

        # CO2 flux has no vertical levels, throws error if you try to select
        # any. Can add any other flux-like variables to this list.
        if afieldName not in ['CO2_gas_flux']:
            iselValues = {'nVertLevels': 0}
        else:
            iselValues = None

        sectionName = 'climatologyMapBGC'

        climStartYear = config.getint('climatology',
                                      'startYear')
        climEndYear = config.getint('climatology',
                                    'endYear')
        # Read in units (since BGC units are variable)
        units = config.getExpression(sectionName + '_' + afieldName, 'units',
                                     elementType=str)

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

        # the variable mpasFieldName will be added to mpasClimatologyTask
        # along with the seasons.
        remapClimatologySubtask = RemapBGCClimatology(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=afieldName,
            variableList=[ampasFieldName],
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=iselValues)

        if refConfig is None:
            refTitleLabel = 'Observations'

            observationsDirectory = build_config_full_path(
                config, 'oceanObservations',
                '{}Subdirectory'.format(afieldName))

            obsFileName = \
                "{}/{}_1.0x1.0degree.nc" .format(
                    observationsDirectory, afieldName)

            observationsLabel = config.getExpression(sectionName + '_' +
                afieldName, 'observationsLabel', elementType=str)
            refFieldName = afieldName
            outFileLabel = afieldName + observationsLabel
            galleryName = 'Observations: ' + observationsLabel

            remapObservationsSubtask = RemapObservedBGCClimatology(
                parentTask=self, seasons=seasons, fileName=obsFileName,
                outFilePrefix=refFieldName,
                comparisonGridNames=comparisonGridNames)
            self.add_subtask(remapObservationsSubtask)
            diffTitleLabel = 'Model - Observations'
        else:
            remapObservationsSubtask = None
            refRunName = refConfig.get('runs', 'mainRunName')
            galleryName = None
            refTitleLabel = 'Ref: {}'.format(refRunName)

            refFieldName = ampasFieldName
            outFileLabel = afieldName
            diffTitleLabel = 'Main - Reference'


        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(self, season,
                                                    comparisonGridName,
                                                    remapClimatologySubtask,
                                                    remapObservationsSubtask,
                                                    refConfig)

                subtask.set_plot_info(
                    outFileLabel=afieldName,
                    fieldNameInTitle=afieldName,
                    mpasFieldName=ampasFieldName,
                    refFieldName=refFieldName,
                    refTitleLabel=refTitleLabel,
                    unitsLabel=units,
                    imageCaption='Mean ' + afieldName,
                    galleryGroup='Sea Surface Biogeochemistry',
                    groupSubtitle=None,
                    groupLink=afieldName,
                    galleryName=afieldName,
                    diffTitleLabel=diffTitleLabel)

                self.add_subtask(subtask)
        # }}}
    # }}}

class RemapBGCClimatology(RemapMpasClimatologySubtask): # {{{
    """
    Apply unit conversions to native model output to align with observations.
    """
    # Authors
    # -------
    # Riley X. Brady

    def customize_remapped_climatology(self, climatology, comparisonGridName,
                                       season): # {{{
        """
        Currently, convert gas flux from native units to mol/m2/yr

        Parameters
        ----------
        climatology : ``xarray.Dataset``
            The MPAS climatology data set that has been remapped

        comparisonGridNames : {'latlon', 'antarctic'}
            The name of the comparison grid to use for remapping

        season : str
            The name of the season to be remapped

        Returns
        -------
        climatology : ``xarray.Dataset``
            The same data set with any custom fields added or modifications
            made
        """
        # Authors
        # -------
        # Riley X. Brady

        fieldName = self.variableList[0]
        # Convert CO2 gas flux from native mmol/m2 m/s to mol/m2/yr for
        # comparison to the SOM-FFN product
        if 'CO2_gas_flux' in fieldName:
            conversion = -1 * (60 * 60 * 24 * 365.25) / 10**3
            climatology[fieldName] = conversion * climatology[fieldName]
        # Convert O2 from mmol/m3 to mL/L for comparison to WOA product
        elif 'O2' in fieldName:
            conversion = 22.391 / 10**3
            climatology[fieldName] = conversion * climatology[fieldName]
        return climatology # }}}

    # }}}

class RemapObservedBGCClimatology(RemapObservedClimatologySubtask): # {{{
    """
    A subtask for reading and remapping BGC observations
    """
    # Authors
    # -------
    # Riley X. Brady

    def get_observation_descriptor(self, fileName): # {{{
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
        # Riley X. Brady

        # create a descriptor of the observation grid using the lat/lon
        # coordinates
        obsDescriptor = LatLonGridDescriptor.read(fileName=fileName,
                                                  latVarName='lat',
                                                  lonVarName='lon')
        return obsDescriptor # }}}

    def build_observational_dataset(self, fileName): # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Parameters
        ----------
        fileName : str
            observation file name

        Returns
        ------
        dsObs : ``xarray.Dataset``
            The observational dataset
        '''
        # Authors
        # -------
        # Riley X. Brady
        # Obs are pre-processed, so nothing needed to be done here.
        dsObs = xr.open_dataset(fileName)
        return dsObs # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
