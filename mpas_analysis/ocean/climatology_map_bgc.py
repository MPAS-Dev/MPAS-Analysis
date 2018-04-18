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
import datetime

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
    Phillip J. Wolfram, Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
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
        Phillip J. Wolfram
        """
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapBGC, self).__init__(
                config=config, taskName='climatologyMapBGC_' + afieldName,
                componentName='ocean',
                tags=['climatology', 'horizontalMap', 'BGC', afieldName])

        iselValues = {'nVertLevels': 0}

        sectionName = 'climatologyMapBGC'

        climStartYear = config.getint('climatology',
                                      'startYear')
        climEndYear = config.getint('climatology',
                                    'endYear')

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
        remapClimatologySubtask = RemapMpasClimatologySubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=afieldName,
            variableList=[ampasFieldName],
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=iselValues)

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(self, season,
                                                    comparisonGridName,
                                                    remapClimatologySubtask,
                                                    None,
                                                    refConfig)

                subtask.set_plot_info(
                        outFileLabel=afieldName,
                        fieldNameInTitle=afieldName,
                        mpasFieldName=ampasFieldName,
                        refFieldName=None,
                        refTitleLabel=None,
                        unitsLabel=r'mmol m$^{-3}$',
                        imageCaption='Mean ' + afieldName,
                        galleryGroup='Sea Surface Biogeochemistry',
                        groupSubtitle=None,
                        groupLink=afieldName,
                        galleryName=afieldName,
                        diffTitleLabel=None)

                self.add_subtask(subtask)
        # }}}
    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
