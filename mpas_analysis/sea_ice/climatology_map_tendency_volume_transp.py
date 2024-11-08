# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#

import xarray as xr
from pyremap import LatLon2DGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask

from mpas_analysis.shared.plot import PlotClimatologyMapSubtask

from mpas_analysis.shared.io.utility import build_obs_path


class ClimatologyMapSeaIceVolumeTendencyTransp(AnalysisTask):
    """
    An analysis task for sea ice volume tendency due to transport
    """
    # Authors
    # -------
    # Elizabeth Hunke, Xylar Asay-Davis

    def __init__(self, config, mpas_climatology_task, hemisphere,
                 control_config=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpas_climatology_task : mpas_analysis.shared.climatology.MpasClimatologyTask
            The task that produced the climatology to be remapped and plotted

        hemisphere : {'NH', 'SH'}
            The hemisphere to plot

        control_config : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Elizabeth Hunke, Xylar Asay-Davis

        task_name = f'climatologyMapSeaIceVolumeTendencyTransp{hemisphere}'

        field_name = 'VolumeTendencyTransport'
        title_long = 'Sea Ice Volume Tendency due to Transport'

        tags = ['climatology', 'horizontalMap', 'publicObs', 'tendency']
        if hemisphere == 'NH':
            tags = tags + ['arctic']
        else:
            tags = tags + ['antarctic']

        # call the constructor from the base class (AnalysisTask)
        super().__init__(config=config, taskName=task_name,
            componentName='seaIce', tags=tags)

        self.mpas_climatology_task = mpas_climatology_task

        section_name = self.taskName

        if hemisphere == 'NH':
            hemisphere_long= 'Northern'
        else:
            hemisphere_long= 'Southern'

        # read in what seasons we want to plot
        seasons = config.getexpression(section_name, 'seasons')

        if len(seasons) == 0:
            raise ValueError(f'config section {section_name} does not contain '
                             f'valid list of seasons')

        comparison_grid_names = config.getexpression(section_name,
                                                   'comparisonGrids')

        if len(comparison_grid_names) == 0:
            raise ValueError(f'config section {section_name} does not contain '
                             f'valid list of comparison grids')

        variable_list = ['timeMonthly_avg_iceVolumeTendencyTransport']

        remap_climatology_subtask = RemapMpasSeaIceTendencyClimatology(
            mpas_climatology_task=mpas_climatology_task,
            parent_task=self,
            climatology_name=f'{field_name}{hemisphere}',
            variable_list=variable_list,
            comparison_grid_names=comparison_grid_names,
            seasons=seasons)

        self.add_subtask(remap_climatology_subtask)

        for season in seasons:
            for comparison_grid_name in comparison_grid_names:

                if control_config is None:
#echmod - no observations yet
#                    if hemisphere == 'SH' and season == 'ANN':
#                        ref_title_label = 'Observations (AnIceFlux)'
#                        gallery_name = None
#                        diff_title_label = 'Model - Observations'
#                        obs_file_name = build_obs_path(
#                                config, 'seaIce',
#                                relativePathOption=f'melting{hemisphere}',
#                                relativePathSection=section_name)

#                        remap_observations_subtask = RemapAnIceFluxMeltingClimatology(
#                            parentTask=self, seasons=seasons,
#                            fileName=obs_file_name,
#                            outFilePrefix=f'{field_name}{hemisphere}',
#                            comparisonGridNames=comparison_grid_names)
#                        self.add_subtask(remap_observations_subtask)
#                    else:
                        remap_observations_subtask = None
                        gallery_name = None
                        ref_title_label = None
                        ref_field_name = None
                        diff_title_label = 'Model - Observations'

                else:
                    control_run_name = control_config.get('runs', 'mainRunName')
                    gallery_name = None
                    ref_title_label = f'Control: {control_run_name}'
                    field_name = field_name
                    diff_title_label = 'Main - Control'

                    remap_observations_subtask = None

                image_caption = f'Climatology Map of ' \
                                f'{hemisphere_long}-Hemisphere ' \
                                f'{title_long}'
                gallery_group = f'{hemisphere_long}-Hemisphere Tendencies due to Transport'
                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(
                    parentTask=self, season=season,
                    comparisonGridName=comparison_grid_name,
                    remapMpasClimatologySubtask=remap_climatology_subtask,
                    remapObsClimatologySubtask=remap_observations_subtask,
                    controlConfig=control_config)

                subtask.set_plot_info(
                    outFileLabel=f'seaice_volumetendtransp{hemisphere}',
                    fieldNameInTitle=title_long,
                    mpasFieldName=field_name,
                    refFieldName=field_name,
                    refTitleLabel=ref_title_label,
                    diffTitleLabel=diff_title_label,
                    unitsLabel=r'm yr$^{-1}$',
                    imageCaption=image_caption,
                    galleryGroup=gallery_group,
                    groupSubtitle=None,
                    groupLink=f'{hemisphere.lower()}_tendency',
                    galleryName=gallery_name,
                    extend='both',
                    prependComparisonGrid=False)

                self.add_subtask(subtask)


class RemapMpasSeaIceTendencyClimatology(RemapMpasClimatologySubtask):
    """
    A subtask for computing climatologies of sea ice tendencies
    """
    def __init__(self, mpas_climatology_task, parent_task, climatology_name,
                 variable_list, seasons, comparison_grid_names):

        """
        Construct the analysis task and add it as a subtask of the
        ``parent_task``.
        Parameters
        ----------
        mpas_climatology_task : mpas_analysis.shared.climatology.MpasClimatologyTask
            The task that produced the climatology to be remapped
        parent_task :  mpas_analysis.shared.AnalysisTask
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``
        climatology_name : str
            A name that describes the climatology (e.g. a short version of
            the important field(s) in the climatology) used to name the
            subdirectories for each stage of the climatology
        variable_list : list of str
            A list of variable names in ``timeSeriesStatsMonthly`` to be
            included in the climatologies
        seasons : list of str, optional
            A list of seasons (keys in ``shared.constants.monthDictionary``)
            to be computed or ['none'] (not ``None``) if only monthly
            climatologies are needed.
        comparison_grid_names : list of {'latlon', 'antarctic'}
            The name(s) of the comparison grid to use for remapping.
        """

        subtask_name = f'remapMpasClimatology_VolumeTendencyTransport'
        # call the constructor from the base class
        # (RemapMpasClimatologySubtask)
        super().__init__(
            mpas_climatology_task, parent_task, climatology_name,
            variable_list, seasons, comparison_grid_names)

        self.mpas_climatology_task = mpas_climatology_task
        self.variable_list = variable_list

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.
        """

        # first, call setup_and_check from the base class
        # (RemapMpasClimatologySubtask), which will set up remappers and add
        # variables to mpas_climatology_task
        super().setup_and_check()

        # don't add the variables and seasons to mpas_climatology_task until
        # we're sure this subtask is supposed to run
        self.mpas_climatology_task.add_variables(self.variable_list,
                                                     self.seasons)

    def customize_masked_climatology(self, climatology, season):
        """
        Compute the tendency
        Parameters
        ----------
        climatology : xarray.Dataset
            the climatology data set
        season : str
            The name of the season to be masked
        Returns
        -------
        climatology : xarray.Dataset
            the modified climatology data set
        """

        tendency = self._compute_tendency(climatology)

        climatology['VolumeTendencyTransport'] = tendency
        climatology.VolumeTendencyTransport.attrs['units'] = 'm/yr'
        climatology = climatology.drop_vars(self.variable_list)

        return climatology

    def _compute_tendency(self, climatology):
        """
        Compute the tendency in m/yr
        """
        ds_restart = xr.open_dataset(self.restartFileName)
        ds_restart = ds_restart.isel(Time=0)

        units_scale_factor = 60 * 60 * 24 * 365

        tendency = climatology['timeMonthly_avg_iceVolumeTendencyTransport'] * units_scale_factor
        return tendency
