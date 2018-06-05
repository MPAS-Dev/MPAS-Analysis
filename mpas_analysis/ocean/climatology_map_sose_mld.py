# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#
'''
Analysis tasks for comparing Antarctic climatology maps against observations
and reanalysis data.
'''
# Authors
# -------
# Xylar Asay-Davis

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask
from mpas_analysis.ocean.remap_sose_climatology import RemapSoseClimatology

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask

from mpas_analysis.shared.io.utility import build_config_full_path


class ClimatologyMapSoseMLD(AnalysisTask):  # {{{
    """
    An analysis task for comparison of antarctic mixed layer depth against SOSE
    fields
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask,
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
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        fieldName = 'mldSOSE'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSoseMLD, self).__init__(
                config=config, taskName='climatologyMapSoseMLD',
                componentName='ocean',
                tags=['climatology', 'horizontalMap', 'sose', 'mld',
                      'publicObs'])

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

        # the variable 'timeMonthly_avg_landIceFreshwaterFlux' will be added to
        # mpasClimatologyTask along with the seasons.
        remapClimatologySubtask = RemapMpasClimatologySubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=fieldName,
            variableList=[mpasFieldName],
            seasons=seasons,
            comparisonGridNames=comparisonGridNames,
            iselValues=iselValues)

        if refConfig is None:

            refTitleLabel = 'State Estimate (SOSE)'

            observationsDirectory = build_config_full_path(
                config, 'oceanObservations', 'soseSubdirectory')

            obsFileName = \
                '{}/SOSE_2005-2010_monthly_mld_6000.0x' \
                '6000.0km_10.0km_Antarctic_stereo.nc'.format(
                        observationsDirectory)
            refFieldName = 'mld'
            outFileLabel = 'mldSOSE'
            galleryName = 'State Estimate: SOSE'
            diffTitleLabel = 'Model - State Estimate'

            remapObservationsSubtask = RemapSoseClimatology(
                    parentTask=self, seasons=seasons, fileName=obsFileName,
                    outFilePrefix='{}SOSE'.format(refFieldName),
                    fieldName=refFieldName,
                    botFieldName=None,
                    depths=None,
                    comparisonGridNames=comparisonGridNames)

            self.add_subtask(remapObservationsSubtask)

        else:
            remapObservationsSubtask = None
            refRunName = refConfig.get('runs', 'mainRunName')
            galleryName = 'Ref: {}'.format(refRunName)
            refTitleLabel = galleryName

            refFieldName = mpasFieldName
            outFileLabel = 'mld'
            diffTitleLabel = 'Main - Reference'

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                subtask = PlotClimatologyMapSubtask(
                    parentTask=self,
                    season=season,
                    comparisonGridName=comparisonGridName,
                    remapMpasClimatologySubtask=remapClimatologySubtask,
                    remapObsClimatologySubtask=remapObservationsSubtask,
                    refConfig=refConfig)

                subtask.set_plot_info(
                    outFileLabel=outFileLabel,
                    fieldNameInTitle='Mixed Layer Depth',
                    mpasFieldName=mpasFieldName,
                    refFieldName=refFieldName,
                    refTitleLabel=refTitleLabel,
                    diffTitleLabel=diffTitleLabel,
                    unitsLabel=r'm',
                    imageCaption='Mixed Layer Depth',
                    galleryGroup='Mixed Layer Depth',
                    groupSubtitle=None,
                    groupLink='mldSose',
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
        super(ClimatologyMapSoseMLD, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_mixedlayerdepths_enable',
            raiseException=True)

        # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
