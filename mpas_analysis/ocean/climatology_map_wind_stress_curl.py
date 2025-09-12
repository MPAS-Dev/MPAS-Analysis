# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/main/LICENSE

import xarray as xr

from mpas_tools.ocean.streamfunction.vorticity import (
    compute_vertically_integrated_vorticity,
)

from mpas_analysis.ocean.utility import (
    vector_cell_to_edge_isotropic,
    vector_to_edge_normal,
)
from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask
from mpas_analysis.shared.plot import PlotClimatologyMapSubtask


class ClimatologyMapWindStressCurl(AnalysisTask):
    """
    An analysis task for computing and plotting maps of the wind stress curl.
    """

    def __init__(self, config, mpas_climatology_task, control_config=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config : tranche.Tranche
            Configuration options

        mpas_climatology_task : mpas_analysis.shared.climatology.MpasClimatologyTask
            The task that produced the climatology to be remapped and plotted

        control_config : tranche.Tranche, optional
            Configuration options for a control run (if any)
        """  # noqa: E501

        field_name = 'windStressCurl'
        super().__init__(
            config=config,
            taskName='climatologyMapWindStressCurl',
            componentName='ocean',
            tags=[
                'climatology', 'horizontalMap', field_name, 'publicObs'
            ],
        )

        section_name = self.taskName

        # read in what seasons we want to plot
        seasons = config.getexpression(section_name, 'seasons')
        if len(seasons) == 0:
            raise ValueError(
                f'config section {section_name} does not contain '
                f'valid list of seasons'
            )

        comparison_grid_names = config.getexpression(
            section_name, 'comparisonGrids'
        )

        if len(comparison_grid_names) == 0:
            raise ValueError(
                f'config section {section_name} does not contain '
                f'valid list of comparison grids'
            )

        variable_list = list(RemapMpasWindStressCurl.VARIABLES)
        remap_climatology_subtask = RemapMpasWindStressCurl(
            mpasClimatologyTask=mpas_climatology_task,
            parentTask=self,
            climatologyName=field_name,
            variableList=variable_list,
            seasons=seasons,
            comparisonGridNames=comparison_grid_names,
            subtaskName='remapWindStressCurl',
            vertices=True,
        )

        self.add_subtask(remap_climatology_subtask)

        out_file_label = field_name
        field_title = 'Wind Stress Curl'
        remap_observations_subtask = None

        mpas_field_name = field_name
        if control_config is None:
            ref_field_name = None
            ref_title_label = None
            diff_title_label = 'Model - Observations'

        else:
            control_run_name = control_config.get('runs', 'mainRunName')
            ref_field_name = mpas_field_name
            ref_title_label = f'Control: {control_run_name}'
            diff_title_label = 'Main - Control'

        for comparison_grid_name in comparison_grid_names:
            for season in seasons:
                # make a new subtask for this season and comparison grid
                subtask_name = f'plot{season}_{comparison_grid_name}'

                subtask = PlotClimatologyMapSubtask(
                    self, season, comparison_grid_name,
                    remap_climatology_subtask, remap_observations_subtask,
                    controlConfig=control_config, subtaskName=subtask_name)
                subtask.set_plot_info(
                    outFileLabel=out_file_label,
                    fieldNameInTitle=field_title,
                    mpasFieldName=mpas_field_name,
                    refFieldName=ref_field_name,
                    refTitleLabel=ref_title_label,
                    diffTitleLabel=diff_title_label,
                    unitsLabel=r'N m$^{-3}$',
                    imageCaption=field_title,
                    galleryGroup='Wind Stress Curl',
                    groupSubtitle=None,
                    groupLink='wsc',
                    galleryName=None,
                    configSectionName=section_name)

                self.add_subtask(subtask)


class RemapMpasWindStressCurl(RemapMpasClimatologySubtask):
    """
    A subtask for computing climatologies of the wind stress curl before
    it gets remapped to a comparison grid.
    """

    VARIABLES = (
        'timeMonthly_avg_windStressZonal',
        'timeMonthly_avg_windStressMeridional',
    )

    def setup_and_check(self):
        """
        Add the variables needed for computing wind stress curl to the
        climatology task
        """
        super().setup_and_check()

        # Add the variables and seasons, now that we have the variable list
        self.mpasClimatologyTask.add_variables(
            list(self.VARIABLES), self.seasons
        )

    def customize_masked_climatology(self, climatology, season):
        """
        Compute the wind stress curl and add it to the climatology.

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
        logger = self.logger

        ds_mesh = xr.open_dataset(self.meshFilename)
        var_list = [
            'verticesOnEdge',
            'cellsOnVertex',
            'kiteAreasOnVertex',
            'angleEdge',
            'areaTriangle',
            'dcEdge',
            'edgesOnVertex',
            'verticesOnEdge',
            'latVertex',
        ]
        ds_mesh = ds_mesh[var_list]

        ws_zonal_cell = climatology['timeMonthly_avg_windStressZonal']
        ws_meridional_cell = (
            climatology['timeMonthly_avg_windStressMeridional']
        )
        ws_zonal_edge, ws_meridional_edge = vector_cell_to_edge_isotropic(
            ds_mesh, ws_zonal_cell, ws_meridional_cell
        )
        ws_normal_edge = vector_to_edge_normal(
            ds_mesh, ws_zonal_edge, ws_meridional_edge
        )

        # despite the name, this is the curl operator
        wind_sress_curl, _ = compute_vertically_integrated_vorticity(
            ds_mesh, ws_normal_edge, logger
        )
        climatology['windStressCurl'] = wind_sress_curl
        climatology['windStressCurl'].attrs['units'] = 'N m-3'

        # drop the original wind stress variables
        climatology = climatology.drop_vars(list(self.VARIABLES))

        return climatology
