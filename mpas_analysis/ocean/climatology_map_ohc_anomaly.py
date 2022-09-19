# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
import xarray as xr
import numpy as np

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask
from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask
from mpas_analysis.ocean.utility import compute_zmid


class ClimatologyMapOHCAnomaly(AnalysisTask):
    """
    An analysis task for comparison of the anomaly from a reference year
    (typically the start of the simulation) of ocean heat content (OHC)

    Attributes
    ----------
    mpas_climatology_task : mpas_analysis.shared.climatology.MpasClimatologyTask
        The task that produced the climatology to be remapped and plotted

    ref_year_climatology_task : mpas_analysis.shared.climatology.RefYearMpasClimatologyTask
        The task that produced the climatology from the first year to be
        remapped and then subtracted from the main climatology
    """

    def __init__(self, config, mpas_climatology_task,
                 ref_year_climatology_task, control_config=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpas_climatology_task : mpas_analysis.shared.climatology.MpasClimatologyTask
            The task that produced the climatology to be remapped and plotted

        ref_year_climatology_task : mpas_analysis.shared.climatology.RefYearMpasClimatologyTask
            The task that produced the climatology from the first year to be
            remapped and then subtracted from the main climatology

        control_config : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """

        field_name = 'deltaOHC'
        # call the constructor from the base class (AnalysisTask)
        super().__init__(config=config, taskName='climatologyMapOHCAnomaly',
                         componentName='ocean',
                         tags=['climatology', 'horizontalMap', field_name,
                               'publicObs', 'anomaly'])

        self.mpas_climatology_task = mpas_climatology_task
        self.ref_year_climatology_task = ref_year_climatology_task

        section_name = self.taskName

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

        depth_ranges = config.getexpression('climatologyMapOHCAnomaly',
                                            'depthRanges',
                                            use_numpyfunc=True)

        mpas_field_name = 'deltaOHC'

        variable_list = ['timeMonthly_avg_activeTracers_temperature',
                         'timeMonthly_avg_layerThickness']

        for min_depth, max_depth in depth_ranges:
            depth_range_string = \
                f'{np.abs(min_depth):g}-{np.abs(max_depth):g}m'
            remap_climatology_subtask = RemapMpasOHCClimatology(
                mpas_climatology_task=mpas_climatology_task,
                ref_year_climatology_task=ref_year_climatology_task,
                parent_task=self,
                climatology_name=f'{field_name}_{depth_range_string}',
                variable_list=variable_list,
                comparison_grid_names=comparison_grid_names,
                seasons=seasons,
                min_depth=min_depth,
                max_depth=max_depth)

            self.add_subtask(remap_climatology_subtask)

            out_file_label = f'deltaOHC_{depth_range_string}'
            remap_observations_subtask = None
            if control_config is None:
                ref_title_label = None
                ref_field_name = None
                diff_title_label = 'Model - Observations'

            else:
                control_run_name = control_config.get('runs', 'mainRunName')
                ref_title_label = f'Control: {control_run_name}'
                ref_field_name = mpas_field_name
                diff_title_label = 'Main - Control'

            for comparison_grid_name in comparison_grid_names:
                for season in seasons:
                    # make a new subtask for this season and comparison grid
                    subtask_name = f'plot{season}_{comparison_grid_name}_{depth_range_string}'

                    subtask = PlotClimatologyMapSubtask(
                        self, season, comparison_grid_name,
                        remap_climatology_subtask, remap_observations_subtask,
                        controlConfig=control_config, subtaskName=subtask_name)

                    subtask.set_plot_info(
                        outFileLabel=out_file_label,
                        fieldNameInTitle=f'$\\Delta$OHC over {depth_range_string}',
                        mpasFieldName=mpas_field_name,
                        refFieldName=ref_field_name,
                        refTitleLabel=ref_title_label,
                        diffTitleLabel=diff_title_label,
                        unitsLabel=r'GJ m$^{-2}$',
                        imageCaption=f'Anomaly in Ocean Heat Content over {depth_range_string}',
                        galleryGroup='OHC Anomaly',
                        groupSubtitle=None,
                        groupLink='ohc_anom',
                        galleryName=None)

                    self.add_subtask(subtask)

    def setup_and_check(self):
        """
        Checks whether analysis is being performed only on the reference year,
        in which case the analysis will not be meaningful.

        Raises
        ------
        ValueError: if attempting to analyze only the reference year
        """

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super().setup_and_check()

        start_year, end_year = self.mpas_climatology_task.get_start_and_end()
        ref_start_year, ref_end_year = \
            self.ref_year_climatology_task.get_start_and_end()

        if (start_year == ref_start_year) and (end_year == ref_end_year):
            raise ValueError('OHC Anomaly is not meaningful and will not work '
                             'when climatology and ref year are the same.')


class RemapMpasOHCClimatology(RemapMpasClimatologySubtask):
    """
    A subtask for computing climatologies of ocean heat content from
    climatologies of temperature

    Attributes
    ----------
    ref_year_climatology_task : mpas_analysis.shared.climatology.RefYearMpasClimatologyTask
        The task that produced the climatology from the first year to be
        remapped and then subtracted from the main climatology

    min_depth, max_depth : float
        The minimum and maximum depths for integration
    """

    def __init__(self, mpas_climatology_task, ref_year_climatology_task,
                 parent_task, climatology_name, variable_list, seasons,
                 comparison_grid_names, min_depth, max_depth):

        """
        Construct the analysis task and adds it as a subtask of the
        ``parent_task``.

        Parameters
        ----------
        mpas_climatology_task : mpas_analysis.shared.climatology.MpasClimatologyTask
            The task that produced the climatology to be remapped

        ref_year_climatology_task : mpas_analysis.shared.climatology.RefYearMpasClimatologyTask
            The task that produced the climatology from the first year to be
            remapped and then subtracted from the main climatology

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

        min_depth, max_depth : float
            The minimum and maximum depths for integration
        """

        depth_range_string = f'{np.abs(min_depth):g}-{np.abs(max_depth):g}m'
        subtask_name = f'remapMpasClimatology_{depth_range_string}'
        # call the constructor from the base class
        # (RemapMpasClimatologySubtask)
        super().__init__(
            mpas_climatology_task, parent_task, climatology_name,
            variable_list, seasons, comparison_grid_names,
            subtaskName=subtask_name)

        self.ref_year_climatology_task = ref_year_climatology_task
        self.run_after(ref_year_climatology_task)
        self.min_depth = min_depth
        self.max_depth = max_depth

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
        self.ref_year_climatology_task.add_variables(self.variableList,
                                                     self.seasons)

    def customize_masked_climatology(self, climatology, season):
        """
        Compute the ocean heat content (OHC) anomaly from the temperature
        and layer thickness fields.

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

        ohc = self._compute_ohc(climatology)
        ref_file_name = self.ref_year_climatology_task.get_file_name(season)
        ref_year_climo = xr.open_dataset(ref_file_name)
        if 'Time' in ref_year_climo.dims:
            ref_year_climo = ref_year_climo.isel(Time=0)
        ref_ohc = self._compute_ohc(ref_year_climo)

        climatology['deltaOHC'] = ohc - ref_ohc
        climatology.deltaOHC.attrs['units'] = 'GJ m^-2'
        start_year = self.ref_year_climatology_task.startYear
        climatology.deltaOHC.attrs['description'] = \
            f'Anomaly from year {start_year} in ocean heat content'
        climatology = climatology.drop_vars(self.variableList)

        return climatology

    def _compute_ohc(self, climatology):
        """
        Compute the OHC from the temperature and layer thicknesses in a given
        climatology data sets.
        """
        ds_restart = xr.open_dataset(self.restartFileName)
        ds_restart = ds_restart.isel(Time=0)

        # specific heat [J/(kg*degC)]
        cp = self.namelist.getfloat('config_specific_heat_sea_water')
        # [kg/m3]
        rho = self.namelist.getfloat('config_density0')

        units_scale_factor = 1e-9

        n_vert_levels = ds_restart.sizes['nVertLevels']

        z_mid = compute_zmid(ds_restart.bottomDepth, ds_restart.maxLevelCell-1,
                             ds_restart.layerThickness)

        vert_index = xr.DataArray.from_dict(
            {'dims': ('nVertLevels',), 'data': np.arange(n_vert_levels)})

        temperature = climatology['timeMonthly_avg_activeTracers_temperature']
        layer_thickness = climatology['timeMonthly_avg_layerThickness']

        masks = [vert_index < ds_restart.maxLevelCell,
                 z_mid <= self.min_depth,
                 z_mid >= self.max_depth]
        for mask in masks:
            temperature = temperature.where(mask)
            layer_thickness = layer_thickness.where(mask)

        ohc = units_scale_factor * rho * cp * layer_thickness * temperature
        ohc = ohc.sum(dim='nVertLevels')
        return ohc
