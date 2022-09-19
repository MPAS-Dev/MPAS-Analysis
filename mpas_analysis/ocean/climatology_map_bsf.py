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


class ClimatologyMapBSF(AnalysisTask):
    """
    An analysis task for computing and plotting maps of the barotropic
    streamfunction (BSF)

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

        remap_climatology_subtask = RemapMpasBSFClimatology(
            mpasClimatologyTask=mpas_climatology_task,
            parentTask=self,
            climatologyName=field_name,
            variableList=variable_list,
            comparisonGridNames=comparison_grid_names,
            seasons=seasons)

        self.add_subtask(remap_climatology_subtask)

        out_file_label = field_name
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
                subtask_name = f'plot{season}_{comparison_grid_name}'

                subtask = PlotClimatologyMapSubtask(
                    self, season, comparison_grid_name,
                    remap_climatology_subtask, remap_observations_subtask,
                    controlConfig=control_config, subtaskName=subtask_name)

                subtask.set_plot_info(
                    outFileLabel=out_file_label,
                    fieldNameInTitle=f'Barotropic Streamfunction',
                    mpasFieldName=mpas_field_name,
                    refFieldName=ref_field_name,
                    refTitleLabel=ref_title_label,
                    diffTitleLabel=diff_title_label,
                    unitsLabel='Sv',
                    imageCaption='Barotropic Streamfunction',
                    galleryGroup='Barotropic Streamfunction',
                    groupSubtitle=None,
                    groupLink='bsf',
                    galleryName=None)

                self.add_subtask(subtask)


class RemapMpasBSFClimatology(RemapMpasClimatologySubtask):
    """
    A subtask for computing climatologies of the barotropic streamfunction
    from climatologies of normal velocity and layer thickness
    """

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
        logger = self.logger

        ds_mesh = xr.open_dataset(self.restartFileName)
        ds_mesh = ds_mesh[['cellsOnEdge', 'cellsOnVertex', 'nEdgesOnCell',
                           'edgesOnCell', 'verticesOnCell', 'verticesOnEdge',
                           'dcEdge', 'dvEdge']]
        ds_mesh.load()

        bsf_vertex = self._compute_barotropic_streamfunction_vertex(
            ds_mesh, climatology)
        logger.info('bsf on vertices computed.')
        bsf_cell = self._compute_barotropic_streamfunction_cell(
            ds_mesh, bsf_vertex)
        logger.info('bsf on cells computed.')

        climatology['barotropicStreamfunction'] = \
            bsf_cell.transpose('Time', 'nCells', 'nVertices')
        climatology.barotropicStreamfunction.attrs['units'] = 'Sv'
        climatology.barotropicStreamfunction.attrs['description'] = \
            'barotropic streamfunction at cell centers'

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

        layer_thickness = ds.timeMonthly_avg_layerThickness
        normal_velocity = \
            ds.timeMonthly_avg_normalVelocity.isel(nEdges=inner_edges)

        layer_thickness_edge = 0.5*(layer_thickness.isel(nCells=cell0) +
                                    layer_thickness.isel(nCells=cell1))
        transport = ds_mesh.dvEdge[inner_edges] * \
            (layer_thickness_edge * normal_velocity).sum(dim='nVertLevels')

        return inner_edges, transport

    def _compute_barotropic_streamfunction_vertex(self, ds_mesh, ds):
        inner_edges, transport = self._compute_transport(ds_mesh, ds)
        print('transport computed.')

        nvertices = ds_mesh.sizes['nVertices']
        ntime = ds.sizes['Time']

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

        bsf_vertex = xr.DataArray(np.zeros((ntime, nvertices)),
                                  dims=('Time', 'nVertices'))

        for tindex in range(ntime):
            rhs = np.zeros(n_inner_edges+n_boundary_vertices, dtype=float)

            # convert to Sv
            ind = np.arange(n_inner_edges)
            rhs[ind] = 1e-6*transport.isel(Time=tindex)

            ind = np.arange(n_boundary_vertices)
            rhs[n_inner_edges + ind] = 0.

            matrix = scipy.sparse.csr_matrix(
                (data, indices),
                shape=(n_inner_edges+n_boundary_vertices, nvertices))

            solution = scipy.sparse.linalg.lsqr(matrix, rhs)

            bsf_vertex[tindex, :] = -solution[0]

        return bsf_vertex

    def _compute_barotropic_streamfunction_cell(self, ds_mesh, bsf_vertex):
        """
        Interpolate the barotropic streamfunction from vertices to cells
        """
        n_edges_on_cell = ds_mesh.nEdgesOnCell
        edges_on_cell = ds_mesh.edgesOnCell - 1
        vertices_on_cell = ds_mesh.verticesOnCell - 1
        area_edge = 0.25*ds_mesh.dcEdge*ds_mesh.dvEdge

        ncells = ds_mesh.sizes['nCells']
        max_edges = ds_mesh.sizes['maxEdges']

        area_vert = xr.DataArray(np.zeros((ncells, max_edges)),
                                 dims=('nCells', 'maxEdges'))

        for ivert in range(max_edges):
            edge_indices = edges_on_cell.isel(maxEdges=ivert)
            mask = ivert < n_edges_on_cell
            area_vert[:, ivert] += 0.5*mask*area_edge.isel(nEdges=edge_indices)

        for ivert in range(max_edges-1):
            edge_indices = edges_on_cell.isel(maxEdges=ivert+1)
            mask = ivert+1 < n_edges_on_cell
            area_vert[:, ivert] += 0.5*mask*area_edge.isel(nEdges=edge_indices)

        edge_indices = edges_on_cell.isel(maxEdges=0)
        mask = n_edges_on_cell == max_edges
        area_vert[:, max_edges-1] += \
            0.5*mask*area_edge.isel(nEdges=edge_indices)

        bsf_cell = \
            ((area_vert * bsf_vertex[:, vertices_on_cell]).sum(dim='maxEdges') /
             area_vert.sum(dim='maxEdges'))

        return bsf_cell
