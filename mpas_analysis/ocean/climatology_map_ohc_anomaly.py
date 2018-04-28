# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray
import numpy

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask

from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask


class ClimatologyMapOHCAnomaly(AnalysisTask):  # {{{
    """
    An analysis task for comparison of the anomaly from a reference year
    (typically the start of the simulation) of ocean heat content (OHC)
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, refYearClimatolgyTask,
                 refConfig=None):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        refYearClimatolgyTask : ``RefYearMpasClimatologyTask``
            The task that produced the climatology from the first year to be
            remapped and then subtracted from the main climatology

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        fieldName = 'deltaOHC'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapOHCAnomaly, self).__init__(
                config=config, taskName='climatologyMapOHCAnomaly',
                componentName='ocean',
                tags=['climatology', 'horizontalMap', fieldName])

        sectionName = self.taskName

        mpasFieldName = 'deltaOHC'

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

        variableList = ['timeMonthly_avg_activeTracers_temperature',
                        'timeMonthly_avg_layerThickness']

        # the variable 'timeMonthly_avg_landIceFreshwaterFlux' will be added to
        # mpasClimatologyTask along with the seasons.
        remapClimatologySubtask = RemapMpasOHCClimatology(
            mpasClimatologyTask=mpasClimatologyTask,
            refYearClimatolgyTask=refYearClimatolgyTask,
            parentTask=self,
            climatologyName=fieldName,
            variableList=variableList,
            comparisonGridNames=comparisonGridNames,
            seasons=seasons)

        self.add_subtask(remapClimatologySubtask)

        outFileLabel = 'deltaOHC'
        remapObservationsSubtask = None
        if refConfig is None:
            refTitleLabel = None
            refFieldName = None
            diffTitleLabel = 'Model - Observations'

        else:
            refRunName = refConfig.get('runs', 'mainRunName')
            refTitleLabel = 'Ref: {}'.format(refRunName)
            refFieldName = mpasFieldName
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
                        outFileLabel=outFileLabel,
                        fieldNameInTitle='OHC Anomaly',
                        mpasFieldName=mpasFieldName,
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=r'GJ m$^{-2}$',
                        imageCaption='Anomaly in Ocean Heat Content',
                        galleryGroup='OHC Anomaly',
                        groupSubtitle=None,
                        groupLink='ohc_anom',
                        galleryName=None)

                self.add_subtask(subtask)
        # }}}
    # }}}


class RemapMpasOHCClimatology(RemapMpasClimatologySubtask):  # {{{
    """
    A subtask for computing climatologies of ocean heat content from
    climatologies of temperature

    Attributes
    ----------
    landIceMask : xarray.DataArray
        A mask indicating where there is land ice on the ocean grid (thus,
        where melt rates are valid)
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, mpasClimatologyTask, refYearClimatolgyTask, parentTask,
                 climatologyName, variableList, seasons, comparisonGridNames):
        # {{{
        '''
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped

        refYearClimatolgyTask : ``RefYearMpasClimatologyTask``
            The task that produced the climatology from the first year to be
            remapped and then subtracted from the main climatology

        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        climatologyName : str
            A name that describes the climatology (e.g. a short version of
            the important field(s) in the climatology) used to name the
            subdirectories for each stage of the climatology

        variableList : list of str
            A list of variable names in ``timeSeriesStatsMonthly`` to be
            included in the climatologies

        seasons : list of str, optional
            A list of seasons (keys in ``shared.constants.monthDictionary``)
            to be computed or ['none'] (not ``None``) if only monthly
            climatologies are needed.

        comparisonGridNames : list of {'latlon', 'antarctic'}
            The name(s) of the comparison grid to use for remapping.

        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # call the constructor from the base class
        # (RemapMpasClimatologySubtask)
        super(RemapMpasOHCClimatology, self).__init__(
                 mpasClimatologyTask, parentTask, climatologyName,
                 variableList, seasons, comparisonGridNames)

        self.refYearClimatolgyTask = refYearClimatolgyTask
        self.run_after(refYearClimatolgyTask)
        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        IOError :
            If a restart file is not available from which to read mesh
            information or if no history files are available from which to
            compute the climatology in the desired time range.
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class
        # (RemapMpasClimatologySubtask), which will set up remappers and add
        # variables to mpasClimatologyTask
        super(RemapMpasOHCClimatology, self).setup_and_check()

        # don't add the variables and seasons to mpasClimatologyTask until
        # we're sure this subtask is supposed to run
        self.refYearClimatolgyTask.add_variables(self.variableList,
                                                 self.seasons)

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

        ohc = self._compute_ohc(climatology)
        refFileName = self.refYearClimatolgyTask.get_file_name(season)
        refYearClimo = xarray.open_dataset(refFileName)
        refYearClimo = refYearClimo.isel(Time=0)
        refOHC = self._compute_ohc(refYearClimo)

        climatology['deltaOHC'] = ohc - refOHC
        climatology.deltaOHC.attrs['units'] = 'GJ m$^{-2}$'
        climatology.deltaOHC.attrs['description'] = \
            'Anomaly from year {} in Ocean heat contentG in each region'.format(
                    self.refYearClimatolgyTask.startYear)

        climatology = climatology.drop(self.variableList)

        return climatology  # }}}

    def _compute_ohc(self, climatology):  # {{{
        """
        Compute the OHC from the temperature and layer thicknesses in a given
        climatology data sets.
        """
        dsRestart = xarray.open_dataset(self.restartFileName)
        dsRestart = dsRestart.isel(Time=0)

        # specific heat [J/(kg*degC)]
        cp = self.namelist.getfloat('config_specific_heat_sea_water')
        # [kg/m3]
        rho = self.namelist.getfloat('config_density0')

        unitsScalefactor = 1e-9

        nVertLevels = dsRestart.sizes['nVertLevels']

        vertIndex = xarray.DataArray.from_dict(
                {'dims': ('nVertLevels',), 'data': numpy.arange(nVertLevels)})

        mask = vertIndex < dsRestart.maxLevelCell

        temperature = \
            climatology['timeMonthly_avg_activeTracers_temperature'].where(
                    mask)
        layerThickness = \
            climatology['timeMonthly_avg_layerThickness'].where(mask)

        ohc = unitsScalefactor*rho*cp * layerThickness * temperature
        ohc = ohc.sum(dim='nVertLevels')
        return ohc  # }}}

    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
