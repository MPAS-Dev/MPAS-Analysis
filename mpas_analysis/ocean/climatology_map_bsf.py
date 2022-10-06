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
import scipy.sparse
import scipy.sparse.linalg

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask
from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask
from mpas_analysis.ocean.utility import compute_zmid


class ClimatologyMapBSF(AnalysisTask):
    """
    An analysis task for computing and plotting maps of the Barotropic /
    Subpolar Gyre streamfunction (BSF / SPGSF)

    Attributes
    ----------
    mpas_climatology_task : mpas_analysis.shared.climatology.MpasClimatologyTask
        The task that produced the climatology to be remapped and plotted
    """

    def __init__(self, config, mpas_climatology_task, control_config=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpas_climatology_task : mpas_analysis.shared.climatology.MpasClimatologyTask
            The task that produced the climatology to be remapped and plotted

        control_config : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        field_name = 'barotropicStreamfunction'
        # call the constructor from the base class (AnalysisTask)
        super().__init__(config=config, taskName='climatologyMapBSF',
                         componentName='ocean',
                         tags=['climatology', 'horizontalMap', field_name,
                               'publicObs', 'streamfunction'])

        self.mpas_climatology_task = mpas_climatology_task

        section_name = self.taskName

        # read in what seasons we want to plot
        seasons = config.getexpression(section_name, 'seasons')
        depth_ranges = config.getexpression(section_name,
                                            'depthRanges',
                                            use_numpyfunc=True)
        if len(seasons) == 0:
            raise ValueError(f'config section {section_name} does not contain '
                             f'valid list of seasons')

        comparison_grid_names = config.getexpression(section_name,
                                                     'comparisonGrids')

        if len(comparison_grid_names) == 0:
            raise ValueError(f'config section {section_name} does not contain '
                             f'valid list of comparison grids')

        mpas_field_name = field_name

        variable_list = ['timeMonthly_avg_normalVelocity',
                         'timeMonthly_avg_layerThickness']
        for min_depth, max_depth in depth_ranges:
            depth_range_string = \
                f'{np.abs(min_depth):g}-{np.abs(max_depth):g}m'
            if np.abs(max_depth) < 6000.:
                fname_title = f'Streamfunction over {depth_range_string}'
                fname_clim = f'{field_name}_{depth_range_string}'
            else:
                fname_title = f'Barotropic Streamfunction'
                fname_clim = field_name

            remap_climatology_subtask = RemapMpasBSFClimatology(
                mpas_climatology_task=mpas_climatology_task,
                parent_task=self,
                climatology_name=fname_clim,
                variable_list=variable_list,
                comparison_grid_names=comparison_grid_names,
                seasons=seasons,
                min_depth=min_depth,
                max_depth=max_depth)

            self.add_subtask(remap_climatology_subtask)

            out_file_label = fname_clim
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
                        fieldNameInTitle=fname_title,
                        mpasFieldName=mpas_field_name,
                        refFieldName=ref_field_name,
                        refTitleLabel=ref_title_label,
                        diffTitleLabel=diff_title_label,
                        unitsLabel='Sv',
                        imageCaption=fname_title,
                        galleryGroup='Horizontal Streamfunction',
                        groupSubtitle=None,
                        groupLink='bsf',
                        galleryName=None)

                    self.add_subtask(subtask)


class RemapMpasBSFClimatology(RemapMpasClimatologySubtask):
    """
    A subtask for computing climatologies of the barotropic /subpolar gyre
    streamfunction from climatologies of normal velocity and layer thickness
    """
    def __init__(self, mpas_climatology_task, parent_task,
                 climatology_name, variable_list, seasons,
                 comparison_grid_names, min_depth, max_depth):

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
            subtaskName=subtask_name, vertices=True)

        self.min_depth = min_depth
        self.max_depth = max_depth

    def customize_masked_climatology(self, climatology, season):
        """
        Compute the masked climatology from the normal velocity
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
        logger = self.logger

        ds_mesh = xr.open_dataset(self.restartFileName)
        ds_mesh = ds_mesh[['cellsOnEdge', 'cellsOnVertex', 'nEdgesOnCell',
                           'edgesOnCell', 'verticesOnCell', 'verticesOnEdge',
                           'dcEdge', 'dvEdge', 'bottomDepth', 'layerThickness',
                           'maxLevelCell']]
        ds_mesh = ds_mesh.isel(Time=0)
        ds_mesh.load()
        bsf_vertex = self._compute_barotropic_streamfunction_vertex(
            ds_mesh, climatology)
        logger.info('bsf on vertices computed.')

        climatology['barotropicStreamfunction'] = bsf_vertex
        climatology.barotropicStreamfunction.attrs['units'] = 'Sv'
        climatology.barotropicStreamfunction.attrs['description'] = \
            'barotropic streamfunction at vertices'

        climatology = climatology.drop_vars(self.variableList)

        return climatology

    def _compute_transport(self, ds_mesh, ds):

        cells_on_edge = ds_mesh.cellsOnEdge - 1
        inner_edges = np.logical_and(cells_on_edge.isel(TWO=0) >= 0,
                                     cells_on_edge.isel(TWO=1) >= 0)

        # convert from boolean mask to indices
        inner_edges = np.flatnonzero(inner_edges.values)

        cell0 = cells_on_edge.isel(nEdges=inner_edges, TWO=0)
        cell1 = cells_on_edge.isel(nEdges=inner_edges, TWO=1)
        n_vert_levels = ds.sizes['nVertLevels']

        vert_index = xr.DataArray.from_dict(
            {'dims': ('nVertLevels',), 'data': np.arange(n_vert_levels)})
        z_mid = compute_zmid(ds_mesh.bottomDepth, ds_mesh.maxLevelCell-1,
                             ds_mesh.layerThickness)
        z_mid_edge = 0.5*(z_mid.isel(nCells=cell0) +
                          z_mid.isel(nCells=cell1))
        normal_velocity = \
            ds.timeMonthly_avg_normalVelocity.isel(nEdges=inner_edges)
        layer_thickness = ds.timeMonthly_avg_layerThickness
        layer_thickness_edge = 0.5*(layer_thickness.isel(nCells=cell0) +
                                    layer_thickness.isel(nCells=cell1))
        mask_bottom = (vert_index < ds_mesh.maxLevelCell).T
        mask_bottom_edge = 0.5*(mask_bottom.isel(nCells=cell0) +
                                mask_bottom.isel(nCells=cell1))
        masks = [mask_bottom_edge,
                 z_mid_edge <= self.min_depth,
                 z_mid_edge >= self.max_depth]
        for mask in masks:
            normal_velocity = normal_velocity.where(mask)
            layer_thickness_edge = layer_thickness_edge.where(mask)
        transport = ds_mesh.dvEdge[inner_edges] * \
            (layer_thickness_edge * normal_velocity).sum(dim='nVertLevels')

        return inner_edges, transport

    def _compute_barotropic_streamfunction_vertex(self, ds_mesh, ds):
        inner_edges, transport = self._compute_transport(ds_mesh, ds)
        self.logger.info('transport computed.')

        nvertices = ds_mesh.sizes['nVertices']

        cells_on_vertex = ds_mesh.cellsOnVertex - 1
        vertices_on_edge = ds_mesh.verticesOnEdge - 1
        is_boundary_cov = cells_on_vertex == -1
        boundary_vertices = np.logical_or(is_boundary_cov.isel(vertexDegree=0),
                                          is_boundary_cov.isel(vertexDegree=1))
        boundary_vertices = np.logical_or(boundary_vertices,
                                          is_boundary_cov.isel(vertexDegree=2))

        # convert from boolean mask to indices
        boundary_vertices = np.flatnonzero(boundary_vertices.values)

        n_boundary_vertices = len(boundary_vertices)
        n_inner_edges = len(inner_edges)

        indices = np.zeros((2, 2*n_inner_edges+n_boundary_vertices), dtype=int)
        data = np.zeros(2*n_inner_edges+n_boundary_vertices, dtype=float)

        # The difference between the streamfunction at vertices on an inner
        # edge should be equal to the transport
        v0 = vertices_on_edge.isel(nEdges=inner_edges, TWO=0).values
        v1 = vertices_on_edge.isel(nEdges=inner_edges, TWO=1).values

        ind = np.arange(n_inner_edges)
        indices[0, 2*ind] = ind
        indices[1, 2*ind] = v1
        data[2*ind] = 1.

        indices[0, 2*ind+1] = ind
        indices[1, 2*ind+1] = v0
        data[2*ind+1] = -1.

        # the streamfunction should be zero at all boundary vertices
        ind = np.arange(n_boundary_vertices)
        indices[0, 2*n_inner_edges + ind] = n_inner_edges + ind
        indices[1, 2*n_inner_edges + ind] = boundary_vertices
        data[2*n_inner_edges + ind] = 1.

        rhs = np.zeros(n_inner_edges+n_boundary_vertices, dtype=float)

        # convert to Sv
        ind = np.arange(n_inner_edges)
        rhs[ind] = 1e-6*transport

        ind = np.arange(n_boundary_vertices)
        rhs[n_inner_edges + ind] = 0.

        matrix = scipy.sparse.csr_matrix(
            (data, indices),
            shape=(n_inner_edges+n_boundary_vertices, nvertices))

        solution = scipy.sparse.linalg.lsqr(matrix, rhs)
        bsf_vertex = xr.DataArray(-solution[0],
                                  dims=('nVertices',))

        return bsf_vertex
