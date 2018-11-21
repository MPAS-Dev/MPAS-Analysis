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
from mpas_analysis.shared.interpolation.utility import add_periodic_lon


class ClimatologyMapBGC(AnalysisTask):  # {{{
    """
    An analysis task for plotting of BGC variables

    Authors
    -------
    Phillip J. Wolfram, Riley X. Brady, Xylar Asay-Davis

    """
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

        Authors
        -------
        Phillip J. Wolfram, Riley X. Brady, Xylar Asay-Davis
        """
        # call the constructor from the base class (AnalysisTask)

        bgcVars = config.getExpression('climatologyMapBGC', 'variables')

        super(ClimatologyMapBGC, self).__init__(
            config=config, taskName='climatologyMapBGC',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', 'BGC'] + bgcVars)

        sectionName = 'climatologyMapBGC'

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

        for fieldName in bgcVars:

            fieldSectionName = '{}_{}'.format(sectionName, fieldName)
            prefix = config.get(fieldSectionName, 'filePrefix')
            mpasFieldName = '{}{}'.format(prefix, fieldName)

            # CO2 flux and pCO2 has no vertical levels, throws error if you try
            # to select any. Can add any other flux-like variables to this
            # list.
            if fieldName not in ['CO2_gas_flux', 'pCO2surface']:
                iselValues = {'nVertLevels': 0}
            else:
                iselValues = None

            # Read in units (since BGC units are variable)
            units = config.get(fieldSectionName, 'units')

            # Pass multiple variables if working with Chlorophyll to sum
            # them to total chlorophyll
            if fieldName == 'Chl':
                prefix = 'timeMonthly_avg_ecosysTracers_'
                variableList = [prefix + 'spChl', prefix + 'diatChl',
                                prefix + 'diazChl', prefix + 'phaeoChl']
                plotField = 'Chl'
            else:
                variableList = [mpasFieldName]
                plotField = mpasFieldName

            remapClimatologySubtask = RemapBGCClimatology(
                mpasClimatologyTask=mpasClimatologyTask,
                parentTask=self,
                climatologyName=fieldName,
                variableList=variableList,
                comparisonGridNames=comparisonGridNames,
                seasons=seasons,
                iselValues=iselValues,
                subtaskName='remapMpasClimatology_{}'.format(fieldName))

            if refConfig is None:
                controlTitleLabel = 'Observations'
                preindustrial = config.getboolean(sectionName, 'preindustrial')
                if preindustrial and 'DIC' in fieldName:
                    controlTitleLabel += ' (Preindustrial)'

                observationsDirectory = build_config_full_path(
                    config, 'oceanObservations',
                    '{}Subdirectory'.format(fieldName))

                # If user wants to compare to preindustrial data, make sure
                # that we load in the right DIC field.
                if preindustrial and fieldName == 'DIC':
                    obsFileName = "{}/PI_DIC_1.0x1.0degree.nc".format(
                            observationsDirectory)
                elif fieldName == 'Chl':
                    obsFileName = \
                        "{}/Chl_SeaWIFS.nc".format(observationsDirectory)
                else:
                    obsFileName = \
                        "{}/{}_1.0x1.0degree.nc" .format(
                            observationsDirectory, fieldName)

                observationsLabel = config.get(fieldSectionName,
                                               'observationsLabel')
                controlFieldName = fieldName
                outFileLabel = fieldName + observationsLabel

                galleryLabel = config.get(fieldSectionName, 'galleryLabel')
                galleryName = '{}  (Compared to {})'.format(galleryLabel,
                                                            observationsLabel)

                remapObservationsSubtask = RemapObservedBGCClimatology(
                    parentTask=self, seasons=seasons, fileName=obsFileName,
                    outFilePrefix=controlFieldName,
                    comparisonGridNames=comparisonGridNames,
                    subtaskName='remapObservations_{}'.format(fieldName))
                self.add_subtask(remapObservationsSubtask)

                diffTitleLabel = 'Model - Observations'

                # Certain BGC observations are only available at annual
                # resolution. Need to ensure that the user is aware that their
                # seasonal or monthly climatology is being compared to ANN.
                # Currently, this is just with GLODAP.
                if observationsLabel == 'GLODAPv2':
                    diffTitleLabel += ' (Compared to ANN)'
            else:
                remapObservationsSubtask = None
                refRunName = refConfig.get('runs', 'mainRunName')
                galleryName = None
                controlTitleLabel = 'Ref: {}'.format(refRunName)

                controlFieldName = mpasFieldName
                outFileLabel = fieldName
                diffTitleLabel = 'Main - Reference'

            for comparisonGridName in comparisonGridNames:
                for season in seasons:
                    # make a new subtask for this season and comparison grid
                    subtask = PlotClimatologyMapSubtask(
                            self, season, comparisonGridName,
                            remapClimatologySubtask, remapObservationsSubtask,
                            refConfig, subtaskName = 'plot{}_{}_{}'.format(
                                fieldName, season, comparisonGridName))

                    subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle=fieldName,
                        mpasFieldName=plotField,
                        controlFieldName=controlFieldName,
                        controlTitleLabel=controlTitleLabel,
                        unitsLabel=units,
                        imageCaption='Mean ' + fieldName,
                        galleryGroup='Sea Surface Biogeochemistry',
                        groupSubtitle=None,
                        groupLink=fieldName,
                        galleryName=galleryName,
                        diffTitleLabel=diffTitleLabel,
                        configSectionName=fieldSectionName)

                    self.add_subtask(subtask)

    def setup_and_check(self):  # {{{
        '''
        Check if preindustrial flag is turned on or off.
        '''
        # Authors
        # -------
        # Riley X. Brady
        super(ClimatologyMapBGC, self).setup_and_check()

        # Clarify that the user is doing preindustrial vs. modern
        preindustrial = self.config.getboolean('climatologyMapBGC',
                                               'preindustrial')
        if preindustrial:
            print("""
                  You are comparing against all available preindustrial
                  datasets. If this is not desired, set the preindustrial
                  flag to 'False' under the ClimatologyMapBGC config section.
                  """)
        else:
            print("""
                  You are comparing against modern observations. If you desire
                  a preindustrial comparison, set the preindustrial flag to
                  'True' under the ClimatologyMapBGC config section.
                  """)

        # }}}
    # }}}


class RemapBGCClimatology(RemapMpasClimatologySubtask):  # {{{
    """
    Apply unit conversions to native model output to align with observations.
    """
    # Authors
    # -------
    # Riley X. Brady

    def customize_masked_climatology(self, climatology, season):  # {{{
        """
        Sum over all phytoplankton chlorophyll to create total chlorophyll

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
        # Riley X. Brady
        climatology = super(RemapBGCClimatology,
                            self).customize_masked_climatology(climatology,
                                                               season)

        if 'timeMonthly_avg_ecosysTracers_spChl' in climatology:
            spChl = climatology.timeMonthly_avg_ecosysTracers_spChl
            diatChl = climatology.timeMonthly_avg_ecosysTracers_diatChl
            diazChl = climatology.timeMonthly_avg_ecosysTracers_diazChl
            phaeoChl = climatology.timeMonthly_avg_ecosysTracers_phaeoChl
            climatology['Chl'] = spChl + diatChl + diazChl + phaeoChl
            climatology.Chl.attrs['units'] = 'mg m$^{-3}$'
            climatology.Chl.attrs['description'] = 'Sum of all PFT chlorophyll'
            climatology.drop(['timeMonthly_avg_ecosysTracers_spChl',
                              'timeMonthly_avg_ecosysTracers_diatChl',
                              'timeMonthly_avg_ecosysTracers_diazChl',
                              'timeMonthly_avg_ecosysTracers_phaeoChl'])

        return climatology  # }}}

    def customize_remapped_climatology(self, climatology, comparisonGridName,
                                       season):  # {{{
        """
        Convert CO2 gas flux from native units to mol/m2/yr,
        Convert dissolved O2 from native units to mL/L

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
        if fieldName == 'timeMonthly_avg_CO2_gas_flux':
            conversion = -1 * (60 * 60 * 24 * 365.25) / 10**3
            climatology[fieldName] = conversion * climatology[fieldName]
        # Convert O2 from mmol/m3 to mL/L for comparison to WOA product
        elif fieldName == 'timeMonthly_avg_ecosysTracers_O2':
            conversion = 22.391 / 10**3
            climatology[fieldName] = conversion * climatology[fieldName]
        return climatology  # }}}

    # }}}


class RemapObservedBGCClimatology(RemapObservedClimatologySubtask):  # {{{
    """
    A subtask for reading and remapping BGC observations
    """
    # Authors
    # -------
    # Riley X. Brady

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
        # Riley X. Brady, Xylar Asay-Davis

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
        ------
        dsObs : ``xarray.Dataset``
            The observational dataset
        '''
        # Authors
        # -------
        # Riley X. Brady, Xylar Asay-Davis

        dsObs = xr.open_dataset(fileName)
        degrees = 'degree' in dsObs.lon.units
        dsObs = add_periodic_lon(ds=dsObs, lonDim='lon', degrees=degrees)
        return dsObs  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
