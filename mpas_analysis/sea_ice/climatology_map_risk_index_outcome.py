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


class ClimatologyMapRiskIndexOutcoume(AnalysisTask):
    """
    An analysis task for evaluating the Risk Index Outcome
    for navigation in sea-ice covered water.
    (https://www.imorules.com/GUID-2C1D86CB-5D58-490F-B4D4-46C057E1D102.html)
    """
    # Authors
    # -------
    # Gennaro D'Angelo, Milena Veneziani

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
        # Gennaro D'Angelo, Milena Veneziani

        task_name = f'climatologyMapRiskIndexOutcome{hemisphere}'

        field_name = 'riskIndexOutcome'

        tags = ['climatology', 'horizontalMap', field_name]
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

        variable_list = ['timeMonthly_avg_iceAreaCell',
                         'timeMonthly_avg_iceVolumeCell']

        remap_climatology_subtask = RemapMpasRiskIndexOutcomeClimatology(
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

                image_caption = f'{season} Climatology Map of ' \
                                f'{hemisphere_long}-Hemisphere Risk Index Outcome'
                gallery_group = f'{hemisphere_long}-Hemisphere Risk Index Outcome'
                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(
                    parentTask=self, season=season,
                    comparisonGridName=comparison_grid_name,
                    remapMpasClimatologySubtask=remap_climatology_subtask,
                    remapObsClimatologySubtask=None,
                    controlConfig=control_config)

                subtask.set_plot_info(
                    outFileLabel=f'risk_index_outcome{hemisphere}',
                    fieldNameInTitle='Risk Index Outcome',
                    mpasFieldName=field_name,
                    refFieldName=field_name,
                    refTitleLabel=ref_title_label,
                    diffTitleLabel=diff_title_label,
                    unitsLabel=r'',
                    imageCaption=image_caption,
                    galleryGroup=gallery_group,
                    groupSubtitle=None,
                    groupLink=f'{hemisphere.lower()}_rio',
                    galleryName=gallery_name,
                    extend=None)

                self.add_subtask(subtask)


class RemapMpasRiskIndexOutcomeClimatology(RemapMpasClimatologySubtask):
    """
    A subtask for computing climatologies of risk index outcome from
    climatologies of sea-ice concentration and thickness.
    """
    def __init__(self, mpas_climatology_task, parent_task, climatology_name,
                 variable_list, seasons, comparison_grid_names):

        """
        Construct the analysis task and adds it as a subtask of the
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

        subtask_name = f'remapMpasClimatology_RiskIndexOutcome'
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
        Compute the risk index outcome from sea-ice concentration and thickness.

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

        rio = self._compute_risk_index_outcome(climatology)

        climatology['riskIndexOutcome'] = rio
        climatology.riskIndexOutcome.attrs['units'] = ''
        climatology = climatology.drop_vars(self.variable_list)

        return climatology

    def _compute_risk_index_outcome(self, climatology):
        """
        Compute the risk index outcome from sea-ice concentration and thickness.
        (https://www.imorules.com/GUID-2C1D86CB-5D58-490F-B4D4-46C057E1D102.html)
        """
        ds_restart = xr.open_dataset(self.restartFileName)
        ds_restart = ds_restart.isel(Time=0)

        scale_factor = 10
        pc = self.config.get(self.taskName, 'polarClass') - 1

        concentration = climatology['timeMonthly_avg_iceAreaCell']
        thickness = climatology['timeMonthly_avg_iceVolumeCell']

#       pic = ["PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "IA Super",\
#              "IA", "IB", "IC", "Not Ice Strengthened"]

        h_riv = np.array([0.5, 10, 15, 30, 50, 70, 100, 120, 170, 200, 250]) * 0.01 
        riv = np.array([[ 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1 ],\
                        [ 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 0 ],\
                        [ 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 0,-1 ],\
                        [ 3, 3, 3, 3, 2, 2, 2, 2, 1, 0,-1,-2 ],\
                        [ 3, 3, 3, 3, 2, 2, 1, 1, 0,-1,-2,-2 ],\
                        [ 3, 2, 2, 2, 2, 1, 1, 0,-1,-2,-3,-3 ],\
                        [ 3, 2, 2, 2, 1, 1, 0,-1,-2,-3,-3,-3 ],\
                        [ 3, 2, 2, 2, 2, 1, 0,-1,-2,-3,-4,-4 ],\
                        [ 3, 2, 2, 2, 1, 0,-1,-2,-3,-4,-5,-5 ],\
                        [ 3, 2, 2, 1, 0,-1,-2,-3,-4,-5,-6,-6 ],\
                        [ 3, 2, 1, 0,-1,-2,-3,-4,-5,-6,-7,-8 ],\
                        [ 3, 1, 0,-1,-2,-3,-4,-5,-6,-7,-8,-8 ]])

        riv_iceCell = np.nan*np.ones(np.shape(thickness))

        riv_mask = np.where(thickness < h_riv[0])
        riv_iceCell[riv_mask] = riv[pc, 0]

        for ind in range(len(h_riv)-1):
#           riv_mask = np.logical_and(thickness >= h_riv[ind], thickness < h_riv[ind+1])
            riv_mask = np.where(thickness >= h_riv[ind])
            riv_iceCell[riv_mask] = riv[pc, ind+1]

        riv_mask = np.where(thickness >= h_riv[-1])
        riv_iceCell[riv_mask] = riv[pc, -1]

        rio = ((1.0 - concentration) * riv[pc, 0] + concentration * riv_iceCell) * units_scale_factor

        return rio
