# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#

import xarray as xr
from pyremap import LatLonGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask

from mpas_analysis.sea_ice.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask

from mpas_analysis.shared.io.utility import build_obs_path


class ClimatologyMapSeaIceProduction(AnalysisTask):
    """
    An analysis task for comparison of sea ice production against
    observations
    """
    # Authors
    # -------
    # Darin Comeau, Xylar Asay-Davis

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
        # Darin Comeau, Xylar Asay-Davis

        task_name = 'climatologyMapSeaIceProduction{}'.format(hemisphere)

        field_name = 'SeaIceProduction'

        tags = ['climatology', 'horizontalMap', field_name, 'publicObs']
        if hemisphere == 'NH':
            tags = tags + ['arctic']
        else:
            tags = tags + ['antarctic']

        # call the constructor from the base class (AnalysisTask)
        super().__init__(config=config, taskName=task_name,
            componentName='seaIce',
            tags=tags)        

        # self.mpas_climatology_task = mpas_climatology_task

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

        variable_list = ['timeMonthly_avg_congelation',
                         'timeMonthly_avg_frazilFormation',
                         'timeMonthly_avg_snowiceFormation']

        remap_climatology_subtask = RemapMpasSeaIceProductionClimatology(
            mpas_climatology_task=mpas_climatology_task,
            parent_task=self,
            climatology_name='{}{}'.format(field_name, hemisphere),
            variable_list=variable_list,
            comparison_grid_names=comparison_grid_names,
            seasons=seasons)

        self.add_subtask(remap_climatology_subtask)

        if control_config is None:
            ref_title_label = 'Observations (Haumann)'
            gallery_name = 'Observations: Haumann'
            diff_title_label = 'Model - Observations'
            field_name = 'SeaIceProduction'
            obs_file_name = build_obs_path(
                    config, 'seaIce',
                    relativePathOption='production{}'.format(
                        hemisphere),
                    relativePathSection=section_name)

            remap_observations_subtask = RemapHaumannProductionClimatology(
                parentTask=self, seasons=seasons,
                fileName=obs_file_name,
                outFilePrefix='{}{}'.format(field_name,
                                            hemisphere),
                comparisonGridNames=comparison_grid_names)
            self.add_subtask(remap_observations_subtask)

        else:
            control_run_name = control_config.get('runs', 'mainRunName')
            gallery_name = None
            ref_title_label = 'Control: {}'.format(controlRunName)
            field_name = field_name
            diff_title_label = 'Main - Control'

            remap_observations_subtask = None

        for season in seasons:
            for comparison_grid_name in comparison_grid_names:
                # # make a new subtask for this season and comparison grid
                # subtask_name = f'plot{field_name}_{season}_{comparison_grid_name}'
                
                image_description = \
                    '{} Climatology Map of {}-Hemisphere Sea Ice ' \
                    'Production.'.format(season, hemisphere_long)
                image_caption = image_description
                gallery_group = \
                    '{}-Hemisphere Sea Ice Production'.format(
                        hemisphere_long)
                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(
                    self, hemisphere, season, comparison_grid_name,
                    remap_climatology_subtask, remap_observations_subtask,
                    controlConfig=control_config)

                subtask.set_plot_info(
                    outFileLabel='seaice_production{}'.format(hemisphere),
                    fieldNameInTitle='Sea Ice Production',
                    mpasFieldName=field_name,
                    refFieldName=field_name,
                    refTitleLabel=ref_title_label,
                    diffTitleLabel=diff_title_label,
                    unitsLabel=r'm yr^-1',
                    imageDescription=image_description,
                    imageCaption=image_caption,
                    galleryGroup=gallery_group,
                    groupSubtitle=None,
                    groupLink='{}_production'.format(hemisphere.lower()),
                    galleryName=gallery_name,
                    maskValue=None)

                self.add_subtask(subtask)


class RemapMpasSeaIceProductionClimatology(RemapMpasClimatologySubtask):
    """
    A subtask for computing climatologies of sea ice production from
    climatologies of individual production terms
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

        subtask_name = f'remapMpasClimatology_SeaIceProduction'
        # call the constructor from the base class
        # (RemapMpasClimatologySubtask)
        super().__init__(
            mpas_climatology_task, parent_task, climatology_name,
            variable_list, seasons, comparison_grid_names)

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
        Compute the total sea ice production from the individual production
        fields.
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

        production = self._compute_production(climatology)

        climatology['SeaIceProduction'] = production
        climatology.production.attrs['units'] = 'm y^-1'
        climatology = climatology.drop_vars(self.variable_list)

        return climatology                                                     

    def _compute_production(self, climatology):
        """
        Compute the total sea ice production in m yr^-1 from the individual
        production fields in m s^-1.
        """
        ds_restart = xr.open_dataset(self.restartFileName)
        ds_restart = ds_restart.isel(Time=0)

        units_scale_factor = 60 * 60 * 24 * 365

        congelation = climatology['timeMonthly_avg_congelation']
        frazil = climatology['timeMonthly_avg_frazilFormation']
        snowice = climatology['timeMonthly_avg_snowiceFormation']

        production = (congelation + frazil + snowice) * units_scale_factor
        return production

class RemapHaumannProductionClimatology(RemapObservedClimatologySubtask):
    """
    A subtask for reading and remapping sea ice production from Haumann
    observations
    """
    # Authors
    # -------
    # Darin Comeau, Xylar Asay-Davis

    def get_observation_descriptor(self, fileName):
        """
        get a MeshDescriptor for the observation grid

        Parameters
        ----------
        fileName : str
            observation file name describing the source grid

        Returns
        -------
        obsDescriptor : ``MeshDescriptor``
            The descriptor for the observation grid
        """
        # Authors
        # -------
        # Darin Comeau, Xylar Asay-Davis

        # create a descriptor of the observation grid using the lat/lon
        # coordinates
        obsDescriptor = LatLonGridDescriptor.read(fileName=fileName,
                                                  latVarName='latitude',
                                                  lonVarName='longitude')
        return obsDescriptor

    def build_observational_dataset(self, fileName):
        """
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Parameters
        ----------
        fileName : str
            observation file name

        Returns
        -------
        dsObs : ``xarray.Dataset``
            The observational dataset
        """
        # Authors
        # -------
        # Darin Comeau, Xylar Asay-Davis

        dsObs = xr.open_dataset(fileName)
        dsObs = dsObs.rename({'freezing': 'SeaIceProduction', 'time': 'Time', 'lat': 'latitude', 'lon': 'longitude'})
        # dsObs.coords['month'] = dsObs['Time.month']
        # dsObs.coords['year'] = dsObs['Time.year']
        dsObs = dsObs.transpose('Time', 'latitude', 'longitude')

        return dsObs
