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

import xarray

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.ocean.remap_depth_slices_subtask import \
    RemapDepthSlicesSubtask
from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask
from mpas_analysis.ocean.remap_sose_climatology import RemapSoseClimatology

from mpas_analysis.shared.io.utility import build_config_full_path


class ClimatologyMapSoseVelMag(AnalysisTask):  # {{{
    """
    An analysis task for comparison of velocity magnitude around Antarctica
    against SOSE fields
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

        fieldName = 'velMagSOSE'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSoseVelMag, self).__init__(
                config=config, taskName='climatologyMapSoseVelMag',
                componentName='ocean',
                tags=['climatology', 'horizontalMap', 'sose', 'velocity',
                      'magnitude'])

        sectionName = self.taskName

        variableList = ['timeMonthly_avg_velocityMeridional',
                        'timeMonthly_avg_velocityZonal']
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

        depths = config.getExpression(sectionName, 'depths')

        if len(depths) == 0:
            raise ValueError('config section {} does not contain valid '
                             'list of depths'.format(sectionName))

        # the variable 'timeMonthly_avg_landIceFreshwaterFlux' will be added to
        # mpasClimatologyTask along with the seasons.
        remapClimatologySubtask = RemapMpasVelMagClimatology(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=fieldName,
            variableList=variableList,
            seasons=seasons,
            depths=depths,
            comparisonGridNames=comparisonGridNames,
            iselValues=iselValues)

        if refConfig is None:

            refTitleLabel = 'State Estimate (SOSE)'

            observationsDirectory = build_config_full_path(
                config, 'oceanObservations', 'soseSubdirectory')

            obsFileName = \
                '{}/SOSE_2005-2010_monthly_vel_mag_6000.0x' \
                '6000.0km_10.0km_Antarctic_stereo.nc'.format(
                        observationsDirectory)
            refFieldName = 'velMag'
            outFileLabel = 'velMagSOSE'
            galleryName = 'State Estimate: SOSE'
            diffTitleLabel = 'Model - State Estimate'

            remapObservationsSubtask = RemapSoseClimatology(
                    parentTask=self, seasons=seasons, fileName=obsFileName,
                    outFilePrefix='{}SOSE'.format(refFieldName),
                    fieldName=refFieldName,
                    botFieldName='botVelMag',
                    depths=depths,
                    comparisonGridNames=comparisonGridNames)

            self.add_subtask(remapObservationsSubtask)

        else:
            remapObservationsSubtask = None
            refRunName = refConfig.get('runs', 'mainRunName')
            galleryName = 'Ref: {}'.format(refRunName)
            refTitleLabel = galleryName

            refFieldName = 'velMag'
            outFileLabel = 'velMag'
            diffTitleLabel = 'Main - Reference'

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                for depth in depths:
                    subtask = PlotClimatologyMapSubtask(
                        parentTask=self,
                        season=season,
                        comparisonGridName=comparisonGridName,
                        remapMpasClimatologySubtask=remapClimatologySubtask,
                        remapObsClimatologySubtask=remapObservationsSubtask,
                        refConfig=refConfig,
                        depth=depth)

                    subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle='Velocity Magnitude',
                        mpasFieldName='velMag',
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=r'm s$^{-1}$',
                        imageCaption='Velocity Magnitude',
                        galleryGroup='Velocity Magnitude',
                        groupSubtitle=None,
                        groupLink='velMagSose',
                        galleryName=galleryName)

                    self.add_subtask(subtask)
        # }}}

    # }}}


class RemapMpasVelMagClimatology(RemapDepthSlicesSubtask):  # {{{
    """
    A subtask for computing climatologies of velocity magnitude from zonal
    and meridional components
    """
    # Authors
    # -------
    # Xylar Asay-Davis

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

        # first, call the base class's version of this function so we extract
        # the desired slices.
        climatology = super(RemapMpasVelMagClimatology,
                            self).customize_masked_climatology(climatology,
                                                               season)

        zonalVel = climatology.timeMonthly_avg_velocityZonal
        meridVel = climatology.timeMonthly_avg_velocityMeridional
        climatology['velMag'] = xarray.ufuncs.sqrt(zonalVel**2 + meridVel**2)
        climatology.velMag.attrs['units'] = 'm s$^{-1}$'
        climatology.velMag.attrs['description'] = 'velocity magnitude'

        climatology = climatology.drop(self.variableList)

        return climatology  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
