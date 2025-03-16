# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#

import xarray as xr
import numpy as np
from pyremap import LatLon2DGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask

from mpas_analysis.shared.plot import PlotClimatologyMapSubtask

from mpas_analysis.shared.io.utility import build_obs_path


class ClimatologyMapRiskIndexOutcome(AnalysisTask):
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

        field_name = 'RiskIndexOutcome'

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
                             f'a valid list of seasons')

        comparison_grid_names = config.getexpression(section_name,
                                                   'comparisonGrids')

        if len(comparison_grid_names) == 0:
            raise ValueError(f'config section {section_name} does not contain '
                             f'a valid list of comparison grids')

        # read in what polar class of vessel we want to plot
        polarclass = config.getexpression(section_name, 'polarClass')

        if polarclass < 0 or polarclass > 12:
            raise ValueError(f'config section {section_name} does not contain '
                             f'a valid instance of Polar Class')

        # convert to 0-based array index
        polarclass = np.int_(polarclass) - 1

        # read in table of Risk Index Values
        #riv_csv = 'riv_MSC.1_Circ.1519_6_June_2016.csv'
        riv_csv = build_obs_path(config, 'seaIce',
                                 relativePathOption='riv{}'.format(hemisphere),
                                 relativePathSection=section_name)
        IceClassLabels = np.genfromtxt(riv_csv, delimiter=',', skip_header=1,
                                       dtype=str, usecols=(0,))
        redClassLabels = np.char.replace(IceClassLabels, ' ', '')
        redClassLabels = np.array([s[:9] for s in redClassLabels])
        riskindexvalue = np.genfromtxt(riv_csv, delimiter=',', skip_header=1)
        riskindexvalue = riskindexvalue[:, 1:]
        riskindexvalue = np.array(riskindexvalue)

        # read in reference floe thicknesses for calculation of Risk Index Values
        # Default values: [0.0, 0.5, 10, 15, 30, 50, 70, 100, 120, 170, 200, 250]/100
        # (the default values were agreed upon by Elizabeth Hunke, Andrew Roberts,
        # and Gennaro D'Angelo based on literature and IMO description)
        h_to_typeofice = config.getexpression(section_name, 'h_to_typeofice')
        h_to_typeofice = np.array(h_to_typeofice)

        # read in what type of variables we want to plot
        useIceCategories = config.getexpression(section_name, 'useIceCategories')
        if useIceCategories:
            variable_list = ['timeMonthly_avg_iceAreaCategory',
                             'timeMonthly_avg_iceVolumeCategory']
        else:
            variable_list = ['timeMonthly_avg_iceAreaCell',
                             'timeMonthly_avg_iceVolumeCell']

        remap_climatology_subtask = RemapMpasRiskIndexOutcomeClimatology(
            mpas_climatology_task=mpas_climatology_task,
            parent_task=self,
            climatology_name=f'{field_name}{hemisphere}',
            variable_list=variable_list,
            comparison_grid_names=comparison_grid_names,
            seasons=seasons,
            polarclass = polarclass,
            riskindexvalue=riskindexvalue,
            h_to_typeofice=h_to_typeofice)
        
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

                image_caption = f'Climatology Map of ' \
                                f'{hemisphere_long}-Hemisphere Risk Index Outcome, ' \
                                f'{IceClassLabels[polarclass]}'
                gallery_group = f'{hemisphere_long}-Hemisphere Risk Index Outcome, ' \
                                f'{IceClassLabels[polarclass]}'

                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(
                    parentTask=self, season=season,
                    comparisonGridName=comparison_grid_name,
                    remapMpasClimatologySubtask=remap_climatology_subtask,
                    remapObsClimatologySubtask=None,
                    controlConfig=control_config)

                subtask.set_plot_info(
                    outFileLabel=f'risk_index_outcome{hemisphere}_' \
                                 f'{redClassLabels[polarclass]}',
                    fieldNameInTitle=f'Risk Index Outcome, ' \
                                     f'{IceClassLabels[polarclass]}',
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
                    extend='min')

                self.add_subtask(subtask)


class RemapMpasRiskIndexOutcomeClimatology(RemapMpasClimatologySubtask):
    """
    A subtask for computing climatologies of Risk Index Outcome from
    climatologies of sea-ice concentration and thickness.
    """
    def __init__(self, mpas_climatology_task, parent_task, climatology_name,
                 variable_list, seasons, comparison_grid_names,
                 polarclass, riskindexvalue, h_to_typeofice):

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
        polarclass : integer
            Polar Class of vessel for which Risk Index Outcomes are computed 
        riskindexvalue : list of integers
            Risk Index Values for a vessel of given Polar Class and type of
            sea ice (values defined by the International Maritime Organization,
            IMO).
        h_to_typeofice: list of reference ice thicknesses
            Values that establish an equvalence between IMO type of ice 
            (and age of ice) and ice floe thickness
        """

        subtask_name = f'remapMpasClimatology_RiskIndexOutcome'
        # call the constructor from the base class
        # (RemapMpasClimatologySubtask)
        super().__init__(
            mpas_climatology_task, parent_task, climatology_name,
            variable_list, seasons, comparison_grid_names)

        self.mpas_climatology_task = mpas_climatology_task
        self.variable_list = variable_list
        self.polarclass = polarclass
        self.riskindexvalue = riskindexvalue
        self.h_to_typeofice = h_to_typeofice

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
        Compute the Risk Index Outcome from sea-ice concentration and (floe) thickness.

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

        climatology['RiskIndexOutcome'] = rio
        climatology.RiskIndexOutcome.attrs['units'] = ''
        climatology = climatology.drop_vars(self.variable_list)

        return climatology

    def _compute_risk_index_outcome(self, climatology):
        """
        Compute the Risk Index Outcome from sea-ice concentration and (floe) thickness,
        as outlined in the International Maritime Organization (IMO) document.
        (https://www.imorules.com/GUID-2C1D86CB-5D58-490F-B4D4-46C057E1D102.html)
        """

        # whether to use sea-ice categories for sea-ice concentration and thickness
        useIceCategories = self.config.getexpression(self.taskName,
                                                     'useIceCategories')

        ds_restart = xr.open_dataset(self.restartFileName)
        ds_restart = ds_restart.isel(Time=0)

        riv = self.riskindexvalue
        h2toi = self.h_to_typeofice
        pc = self.polarclass

        # sea-ice concentration conversion from range [0,1] to range [0,10]
        scale_factor = 10

        if useIceCategories:
            concentration = climatology['timeMonthly_avg_iceAreaCategory']
            thickness = climatology['timeMonthly_avg_iceVolumeCategory']
        else:
            concentration = climatology['timeMonthly_avg_iceAreaCell']
            thickness = climatology['timeMonthly_avg_iceVolumeCell']

        # remove 1-dimensional coordinates
        concentration = concentration.squeeze()
        thickness = thickness.squeeze()

        # correct out-of-range values
        concentration = np.clip(concentration, a_min=0, a_max=1)
        # compute sea-ice floe thickness (ice thickness averaged
        # over the fraction of the cell area covered by ice)
        with np.errstate(divide='ignore', invalid='ignore'):
            thickness = np.where(concentration > 0,
                                 np.divide(thickness, concentration), 0)

        # define RIV array and set RIV values according to sea-ice
        # floe thickness and reference values. There are as many RIV
        # values as there are cells and ice categories
        riv_iceCell = np.full_like(thickness, np.nan)
        for ind in range(len(h2toi)):
            riv_mask = np.where(thickness >= h2toi[ind])
            riv_iceCell[riv_mask] = riv[pc, ind]

        if concentration.ndim > 1:
            # Risk Index Outcome for multi-category ice. The RIO is derived
            # as a concentration-weighted sum of RIVs over each ice category
            siconc = concentration.sum(dim='nCategories')
            siconc = np.clip(siconc, a_min=0, a_max=1)
            rio_ow = (1.0 - siconc) * riv[pc, 0]
            rio_si = concentration * riv_iceCell
            rio_si = rio_si.sum(dim='nCategories')
            rio = rio_ow + rio_si
        else:
            # Risk Index Outcome for single-category ice. There are only two
            # terms per cell: one coming from the fraction of the cell covered
            # by open water and one coming from the fraction covered by sea ice
            rio = (1.0 - concentration) * riv[pc, 0] + concentration * riv_iceCell

        # out-of-range corrections
        rio = np.clip(rio, a_min=riv[pc, -1], a_max=riv[pc, 0])
        # re-scaling
        rio *= scale_factor

        return rio
