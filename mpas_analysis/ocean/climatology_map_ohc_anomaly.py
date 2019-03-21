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

import xarray
import numpy

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask

from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask

from mpas_analysis.ocean.utility import compute_zmid


class ClimatologyMapOHCAnomaly(AnalysisTask):  # {{{
    """
    An analysis task for comparison of the anomaly from a reference year
    (typically the start of the simulation) of ocean heat content (OHC)

    Attributes
    ----------
    mpasClimatologyTask : ``MpasClimatologyTask``
        The task that produced the climatology to be remapped and plotted

    refYearClimatolgyTask : ``RefYearMpasClimatologyTask``
        The task that produced the climatology from the first year to be
        remapped and then subtracted from the main climatology
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, refYearClimatolgyTask,
                 controlConfig=None):  # {{{
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

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        fieldName = 'deltaOHC'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapOHCAnomaly, self).__init__(
            config=config, taskName='climatologyMapOHCAnomaly',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', fieldName, 'publicObs'])

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

        depthRanges = config.getExpression('climatologyMapOHCAnomaly',
                                           'depthRanges',
                                           usenumpyfunc=True)

        self.mpasClimatologyTask = mpasClimatologyTask
        self.refYearClimatolgyTask = refYearClimatolgyTask

        for minDepth, maxDepth in depthRanges:
            depthRangeString = '{:g}-{:g}m'.format(numpy.abs(minDepth),
                                                   numpy.abs(maxDepth))
            remapClimatologySubtask = RemapMpasOHCClimatology(
                mpasClimatologyTask=mpasClimatologyTask,
                refYearClimatolgyTask=refYearClimatolgyTask,
                parentTask=self,
                climatologyName='{}_{}'.format(fieldName, depthRangeString),
                variableList=variableList,
                comparisonGridNames=comparisonGridNames,
                seasons=seasons,
                minDepth=minDepth,
                maxDepth=maxDepth)

            self.add_subtask(remapClimatologySubtask)

            outFileLabel = 'deltaOHC_{}'.format(depthRangeString)
            remapObservationsSubtask = None
            if controlConfig is None:
                refTitleLabel = None
                refFieldName = None
                diffTitleLabel = 'Model - Observations'

            else:
                controlRunName = controlConfig.get('runs', 'mainRunName')
                refTitleLabel = 'Control: {}'.format(controlRunName)
                refFieldName = mpasFieldName
                diffTitleLabel = 'Main - Control'

            for comparisonGridName in comparisonGridNames:
                for season in seasons:
                    # make a new subtask for this season and comparison grid
                    subtaskName = 'plot{}_{}_{}'.format(
                        season, comparisonGridName, depthRangeString)

                    subtask = PlotClimatologyMapSubtask(
                        self, season, comparisonGridName,
                        remapClimatologySubtask, remapObservationsSubtask,
                        controlConfig, subtaskName=subtaskName)

                    subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle=r'$\Delta$OHC over {}'.format(
                            depthRangeString),
                        mpasFieldName=mpasFieldName,
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=r'GJ m$^{-2}$',
                        imageCaption='Anomaly in Ocean Heat Content over '
                        '{}'.format(depthRangeString),
                        galleryGroup='OHC Anomaly',
                        groupSubtitle=None,
                        groupLink='ohc_anom',
                        galleryName=None)

                    self.add_subtask(subtask)
        # }}}

    def setup_and_check(self):  # {{{
        '''
        Checks whether analysis is being performed only on the reference year,
        in which case the analysis will not be meaningful.

        Raises
        ------
        ValueError: if attempting to analyze only the reference year
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ClimatologyMapOHCAnomaly, self).setup_and_check()

        startYear, endYear = self.mpasClimatologyTask.get_start_and_end()
        refStartYear, refEndYear = \
            self.refYearClimatolgyTask.get_start_and_end()

        if (startYear == refStartYear) and (endYear == refEndYear):
            raise ValueError('OHC Anomaly is not meaningful and will not work '
                             'when climatology and ref year are the same.')

        # }}}
    # }}}


class RemapMpasOHCClimatology(RemapMpasClimatologySubtask):  # {{{
    """
    A subtask for computing climatologies of ocean heat content from
    climatologies of temperature

    Attributes
    ----------
    minDepth, maxDepth : float
        The minimum and maximum depths for integration
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, mpasClimatologyTask, refYearClimatolgyTask, parentTask,
                 climatologyName, variableList, seasons, comparisonGridNames,
                 minDepth, maxDepth):
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

        minDepth, maxDepth : float
            The minimum and maximum depths for integration

        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        subtaskName = 'remapMpasClimatology_{:g}-{:g}m'.format(
            numpy.abs(minDepth), numpy.abs(maxDepth))
        # call the constructor from the base class
        # (RemapMpasClimatologySubtask)
        super(RemapMpasOHCClimatology, self).__init__(
            mpasClimatologyTask, parentTask, climatologyName,
            variableList, seasons, comparisonGridNames,
            subtaskName=subtaskName)

        self.refYearClimatolgyTask = refYearClimatolgyTask
        self.run_after(refYearClimatolgyTask)
        self.minDepth = minDepth
        self.maxDepth = maxDepth
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
        if 'Time' in refYearClimo.dims:
            refYearClimo = refYearClimo.isel(Time=0)
        refOHC = self._compute_ohc(refYearClimo)

        climatology['deltaOHC'] = ohc - refOHC
        climatology.deltaOHC.attrs['units'] = 'GJ m$^{-2}$'
        climatology.deltaOHC.attrs['description'] = \
            'Anomaly from year {} in ocean heat content'.format(
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

        zMid = compute_zmid(dsRestart.bottomDepth, dsRestart.maxLevelCell,
                            dsRestart.layerThickness)

        vertIndex = xarray.DataArray.from_dict(
            {'dims': ('nVertLevels',), 'data': numpy.arange(nVertLevels)})

        temperature = climatology['timeMonthly_avg_activeTracers_temperature']
        layerThickness = climatology['timeMonthly_avg_layerThickness']

        masks = [vertIndex < dsRestart.maxLevelCell,
                 zMid <= self.minDepth,
                 zMid >= self.maxDepth]
        for mask in masks:
            temperature = temperature.where(mask)
            layerThickness = layerThickness.where(mask)

        ohc = unitsScalefactor * rho * cp * layerThickness * temperature
        ohc = ohc.sum(dim='nVertLevels')
        return ohc  # }}}

    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
