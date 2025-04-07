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
import os

import xarray as xr
import numpy as np
import scipy.sparse
import scipy.sparse.linalg

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask
from mpas_analysis.shared.plot import PlotClimatologyMapSubtask
from mpas_analysis.ocean.utility import compute_zmid
from mpas_analysis.shared.projection import comparison_grid_option_suffixes


class ClimatologyMapBSF(AnalysisTask):
    """
    An analysis task for computing and plotting maps of the Barotropic /
    Subpolar Gyre streamfunction (BSF / SPGSF)

    Attributes
    ----------
    mpas_climatology_task : mpas_analysis.shared.climatology.MpasClimatologyTask
        The task that produced the climatology to be remapped and plotted
    """  # noqa: E501

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
        """  # noqa: E501
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

        for min_depth, max_depth in depth_ranges:
            depth_range_string = f'{min_depth:g}_to_{max_depth:g}m'
            if max_depth >= -6000. or min_depth < 0.:
                title_depth = depth_range_string.replace('_', ' ')
                fname_title = f'Streamfunction over {title_depth}'
                fname_clim = f'{field_name}_{depth_range_string}'
            else:
                fname_title = 'Barotropic Streamfunction'
                fname_clim = field_name

            remap_climatology_subtask = RemapMpasBSFClimatology(
                mpas_climatology_task=mpas_climatology_task,
                parent_task=self,
                climatology_name=fname_clim,
                comparison_grid_names=comparison_grid_names,
                seasons=seasons,
                min_depth=min_depth,
                max_depth=max_depth)

            self.add_subtask(remap_climatology_subtask)

            out_file_label = fname_clim
            remap_observations_subtask = None
            if control_config is None:
                ref_title_label = None
                diff_title_label = 'Model - Observations'

            else:
                control_run_name = control_config.get('runs', 'mainRunName')
                ref_title_label = f'Control: {control_run_name}'
                diff_title_label = 'Main - Control'
            for comparison_grid_name in comparison_grid_names:
                grid_suffix = \
                    comparison_grid_option_suffixes[comparison_grid_name]
                config_section_name = f'{self.taskName}{grid_suffix}'
                if config.has_section(config_section_name):
                    # if this comparison grid has its own section, there is a
                    # version of the BSF that has been offset for this region
                    # and an associated colorbar/colormap
                    mpas_field_name = f'{field_name}{grid_suffix}'
                else:
                    config_section_name = self.taskName
                    mpas_field_name = field_name
                if control_config is None:
                    ref_field_name = None
                else:
                    ref_field_name = mpas_field_name

                for season in seasons:
                    # make a new subtask for this season and comparison grid
                    subtask_name = f'plot{season}_{comparison_grid_name}' \
                                   f'_{depth_range_string}'

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
                        galleryName=None,
                        configSectionName=config_section_name)

                    self.add_subtask(subtask)


class RemapMpasBSFClimatology(RemapMpasClimatologySubtask):
    """
    A subtask for computing climatologies of the barotropic /subpolar gyre
    streamfunction from climatologies of normal velocity and layer thickness

    Attributes
    ----------
    min_depth : float
        The minimum depth for integration

    max_depth : float
        The maximum depth for integration

    include_bolus : bool
        Whether to include the bolus velocity in the streamfunction

    include_submesoscale : bool
        Whether to include the submesoscale velocity in the streamfunction
    """
    def __init__(self, mpas_climatology_task, parent_task,
                 climatology_name, seasons, comparison_grid_names, min_depth,
                 max_depth):

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

        seasons : list of str, optional
            A list of seasons (keys in ``shared.constants.monthDictionary``)
            to be computed or ['none'] (not ``None``) if only monthly
            climatologies are needed.

        comparison_grid_names : list of {'latlon', 'antarctic'}
            The name(s) of the comparison grid to use for remapping.

        min_depth, max_depth : float
            The minimum and maximum depths for integration
        """  # noqa: E501

        # we'll fill this in at setup time
        variable_list = []

        depth_range_string = f'{min_depth:g}_to_{max_depth:g}m'
        subtask_name = f'remapMpasClimatology_{depth_range_string}'
        # call the constructor from the base class
        # (RemapMpasClimatologySubtask)
        super().__init__(
            mpas_climatology_task, parent_task, climatology_name,
            variable_list, seasons, comparison_grid_names,
            subtaskName=subtask_name, vertices=True)

        self.min_depth = min_depth
        self.max_depth = max_depth
        self.include_bolus = None
        self.include_submesoscale = None

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.
        """
        super().setup_and_check()

        variable_list = ['timeMonthly_avg_normalVelocity',
                         'timeMonthly_avg_layerThickness']

        # Add the bolus velocity if GM is enabled
        try:
            # the new name
            self.include_bolus = self.namelist.getbool('config_use_gm')
        except KeyError:
            # the old name
            self.include_bolus = self.namelist.getbool(
                'config_use_standardgm')
        try:
            self.include_submesoscale = \
                self.namelist.getbool('config_submesoscale_enable')
        except KeyError:
            # an old run without submesoscale
            self.include_submesoscale = False

        if self.include_bolus:
            variable_list.append('timeMonthly_avg_normalGMBolusVelocity')

        if self.include_submesoscale:
            variable_list.append('timeMonthly_avg_normalMLEvelocity')

        self.variableList = variable_list

        # Add the variables and seasons, now that we have the variable list
        self.mpasClimatologyTask.add_variables(self.variableList, self.seasons)

    def run_task(self):
        """
        Compute the requested climatologies
        """
        config = self.config
        # check if climatology exists and if all comparison grids are present
        for season in self.seasons:
            masked_climatology_filename = self.get_masked_file_name(season)
            if not os.path.exists(masked_climatology_filename):
                continue
            all_found = True
            with xr.open_dataset(masked_climatology_filename) as ds:
                for comparison_grid_name in self.comparisonDescriptors.keys():
                    grid_suffix = \
                        comparison_grid_option_suffixes[comparison_grid_name]
                    config_section_name = f'{self.taskName}{grid_suffix}'
                    if config.has_section(config_section_name):
                        mpas_field_name = \
                            f'barotropicStreamfunction{grid_suffix}'
                        if mpas_field_name not in ds:
                            all_found = False
                            break
            if not all_found:
                # if not, remove the files and recompute/remap
                os.remove(masked_climatology_filename)
                for comparison_grid_name in self.comparisonDescriptors.keys():
                    remapped_filename = self.get_remapped_file_name(
                        season, comparison_grid_name)
                    if os.path.exists(remapped_filename):
                        os.remove(remapped_filename)

        super().run_task()

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
        config = self.config

        ds_mesh = xr.open_dataset(self.restartFileName)
        ds_mesh = ds_mesh[['cellsOnEdge', 'cellsOnVertex', 'nEdgesOnCell',
                           'edgesOnCell', 'verticesOnCell', 'verticesOnEdge',
                           'edgesOnVertex', 'dcEdge', 'dvEdge', 'bottomDepth',
                           'maxLevelCell', 'latVertex', 'areaTriangle',]]
        ds_mesh.load()
        bsf_vertex = self._compute_barotropic_streamfunction_vertex(
            ds_mesh, climatology)
        logger.info('bsf on vertices computed.')

        climatology['barotropicStreamfunction'] = bsf_vertex
        climatology.barotropicStreamfunction.attrs['units'] = 'Sv'
        climatology.barotropicStreamfunction.attrs['description'] = \
            'barotropic streamfunction at vertices'

        climatology = climatology.drop_vars(self.variableList)

        # offset the BSF for specific comparison grids if defined
        for comparison_grid_name in self.comparisonDescriptors.keys():
            grid_suffix = \
                comparison_grid_option_suffixes[comparison_grid_name]
            config_section_name = f'{self.taskName}{grid_suffix}'
            if config.has_section(config_section_name):
                # if this comparison grid has its own section, there is a
                # version of the BSF that has been offset for this region
                # and an associated colorbar/colormap
                mpas_field_name = \
                    f'barotropicStreamfunction{grid_suffix}'

                lat_range = config.getexpression(
                    config_section_name, 'latitudeRangeForZeroBSF')
                climatology[mpas_field_name] = _shift_bsf(
                    bsf_vertex, lat_range, ds_mesh.cellsOnVertex - 1,
                    ds_mesh.latVertex)
                climatology[mpas_field_name].attrs['units'] = 'Sv'
                climatology[mpas_field_name].attrs['description'] = \
                    f'barotropic streamfunction at vertices, offset for ' \
                    f'{grid_suffix} plots'

        return climatology

    def _compute_vert_integ_velocity(self, ds_mesh, ds):

        cells_on_edge = ds_mesh.cellsOnEdge - 1
        inner_edges = np.logical_and(cells_on_edge.isel(TWO=0) >= 0,
                                     cells_on_edge.isel(TWO=1) >= 0)

        # convert from boolean mask to indices
        inner_edges = np.flatnonzero(inner_edges.values)

        cell0 = cells_on_edge.isel(nEdges=inner_edges, TWO=0)
        cell1 = cells_on_edge.isel(nEdges=inner_edges, TWO=1)
        n_vert_levels = ds.sizes['nVertLevels']

        layer_thickness = ds.timeMonthly_avg_layerThickness
        max_level_cell = ds_mesh.maxLevelCell - 1

        vert_index = xr.DataArray.from_dict(
            {'dims': ('nVertLevels',), 'data': np.arange(n_vert_levels)})
        z_mid = compute_zmid(ds_mesh.bottomDepth, max_level_cell,
                             layer_thickness)
        z_mid_edge = 0.5*(z_mid.isel(nCells=cell0) +
                          z_mid.isel(nCells=cell1))

        normal_velocity = ds.timeMonthly_avg_normalVelocity
        if self.include_bolus:
            normal_velocity += ds.timeMonthly_avg_normalGMBolusVelocity
        if self.include_submesoscale:
            normal_velocity += ds.timeMonthly_avg_normalMLEvelocity
        normal_velocity = normal_velocity.isel(nEdges=inner_edges)

        layer_thickness_edge = 0.5*(layer_thickness.isel(nCells=cell0) +
                                    layer_thickness.isel(nCells=cell1))
        mask_bottom = (vert_index <= max_level_cell).T
        mask_bottom_edge = np.logical_and(mask_bottom.isel(nCells=cell0),
                                          mask_bottom.isel(nCells=cell1))
        masks = [mask_bottom_edge,
                 z_mid_edge <= self.min_depth,
                 z_mid_edge >= self.max_depth]
        for mask in masks:
            normal_velocity = normal_velocity.where(mask)
            layer_thickness_edge = layer_thickness_edge.where(mask)

        vert_integ_velocity = np.zeros(ds_mesh.dims['nEdges'], dtype=float)
        inner_vert_integ_vel = (
            (layer_thickness_edge * normal_velocity).sum(dim='nVertLevels'))
        vert_integ_velocity[inner_edges] = inner_vert_integ_vel.values

        vert_integ_velocity = xr.DataArray(vert_integ_velocity,
                                           dims=('nEdges',))

        return vert_integ_velocity

    def _compute_edge_sign_on_vertex(self, ds_mesh):
        edges_on_vertex = ds_mesh.edgesOnVertex - 1
        vertices_on_edge = ds_mesh.verticesOnEdge - 1

        nvertices = ds_mesh.sizes['nVertices']
        vertex_degree = ds_mesh.sizes['vertexDegree']

        edge_sign_on_vertex = np.zeros((nvertices, vertex_degree), dtype=int)
        vertices = np.arange(nvertices)
        for iedge in range(vertex_degree):
            eov = edges_on_vertex.isel(vertexDegree=iedge)
            valid_edge = eov >= 0

            v0_on_edge = vertices_on_edge.isel(nEdges=eov, TWO=0)
            v1_on_edge = vertices_on_edge.isel(nEdges=eov, TWO=1)
            valid_edge = np.logical_and(valid_edge, v0_on_edge >= 0)
            valid_edge = np.logical_and(valid_edge, v1_on_edge >= 0)

            mask = np.logical_and(valid_edge, v0_on_edge == vertices)
            edge_sign_on_vertex[mask, iedge] = -1

            mask = np.logical_and(valid_edge, v1_on_edge == vertices)
            edge_sign_on_vertex[mask, iedge] = 1

        return edge_sign_on_vertex

    def _compute_vert_integ_vorticity(self, ds_mesh, vert_integ_velocity,
                                      edge_sign_on_vertex):

        area_vertex = ds_mesh.areaTriangle
        dc_edge = ds_mesh.dcEdge
        edges_on_vertex = ds_mesh.edgesOnVertex - 1

        vertex_degree = ds_mesh.sizes['vertexDegree']

        vert_integ_vorticity = xr.zeros_like(ds_mesh.latVertex)
        for iedge in range(vertex_degree):
            eov = edges_on_vertex.isel(vertexDegree=iedge)
            edge_sign = edge_sign_on_vertex[:, iedge]
            dc = dc_edge.isel(nEdges=eov)
            vert_integ_vel = vert_integ_velocity.isel(nEdges=eov)
            vert_integ_vorticity += (
                dc / area_vertex * edge_sign * vert_integ_vel)

        return vert_integ_vorticity

    def _compute_barotropic_streamfunction_vertex(self, ds_mesh, ds):
        edge_sign_on_vertex = self._compute_edge_sign_on_vertex(ds_mesh)
        vert_integ_velocity = self._compute_vert_integ_velocity(ds_mesh, ds)
        vert_integ_vorticity = self._compute_vert_integ_vorticity(
            ds_mesh, vert_integ_velocity, edge_sign_on_vertex)
        self.logger.info('vertically integrated vorticity computed.')

        config = self.config
        lat_range = config.getexpression(
            'climatologyMapBSF', 'latitudeRangeForZeroBSF')

        nvertices = ds_mesh.sizes['nVertices']
        vertex_degree = ds_mesh.sizes['vertexDegree']

        cells_on_vertex = ds_mesh.cellsOnVertex - 1
        edges_on_vertex = ds_mesh.edgesOnVertex - 1
        vertices_on_edge = ds_mesh.verticesOnEdge - 1
        area_vertex = ds_mesh.areaTriangle
        dc_edge = ds_mesh.dcEdge
        dv_edge = ds_mesh.dvEdge

        # one equation involving vertex degree + 1 vertices for each vertex
        # plus 2 entries for the boundary condition and Lagrange multiplier
        ndata = (vertex_degree + 1) * nvertices + 2
        indices = np.zeros((2, ndata), dtype=int)
        data = np.zeros(ndata, dtype=float)

        # the laplacian on the dual mesh of the streamfunction is the
        # vertically integrated vorticity
        vertices = np.arange(nvertices, dtype=int)
        idata = (vertex_degree + 1) * vertices + 1
        indices[0, idata] = vertices
        indices[1, idata] = vertices
        for iedge in range(vertex_degree):
            eov = edges_on_vertex.isel(vertexDegree=iedge)
            dc = dc_edge.isel(nEdges=eov)
            dv = dv_edge.isel(nEdges=eov)

            v0 = vertices_on_edge.isel(nEdges=eov, TWO=0)
            v1 = vertices_on_edge.isel(nEdges=eov, TWO=1)

            edge_sign = edge_sign_on_vertex[:, iedge]

            mask = v0 == vertices
            # the difference is v1 - v0, so we want to subtract this vertex
            # when it is v0 and add it when it is v1
            this_vert_sign = np.where(mask, -1., 1.)
            # the other vertex is obviously whichever one this is not
            other_vert_index = np.where(mask, v1, v0)
            # if there are invalid vertices, we need to make sure we don't
            # index out of bounds.  The edge_sign will mask these out
            other_vert_index = np.where(other_vert_index >= 0,
                                        other_vert_index, 0)

            idata_other = idata + iedge + 1

            indices[0, idata] = vertices
            indices[1, idata] = vertices
            indices[0, idata_other] = vertices
            indices[1, idata_other] = other_vert_index

            this_data = this_vert_sign * edge_sign * dc / (dv * area_vertex)
            data[idata] += this_data
            data[idata_other] = -this_data

        # Now, the boundary condition: To begin with, we set the BSF at the
        # frist vertext to zero
        indices[0, -2] = nvertices
        indices[1, -2] = 0
        data[-2] = 1.

        # The same in the final column
        indices[0, -1] = 0
        indices[1, -1] = nvertices
        data[-1] = 1.

        # one extra spot for the Lagrange multiplier
        rhs = np.zeros(nvertices + 1, dtype=float)

        rhs[0:-1] = vert_integ_vorticity.values

        matrix = scipy.sparse.csr_matrix(
            (data, indices),
            shape=(nvertices + 1, nvertices + 1))

        solution = scipy.sparse.linalg.spsolve(matrix, rhs)

        # drop the Lagrange multiplier and convert to Sv with the desired sign
        # convention
        bsf_vertex = xr.DataArray(-1e-6 * solution[0:-1],
                                  dims=('nVertices',))

        bsf_vertex = _shift_bsf(bsf_vertex, lat_range, cells_on_vertex,
                                ds_mesh.latVertex)

        return bsf_vertex


def _shift_bsf(bsf_vertex, lat_range, cells_on_vertex, lat_vertex):
    """
    Shift the barotropic streamfunction to be zero at the boundary over
    the given latitude range
    """
    is_boundary_cov = cells_on_vertex == -1
    boundary_vertices = is_boundary_cov.sum(dim='vertexDegree') > 0

    boundary_vertices = np.logical_and(
        boundary_vertices,
        lat_vertex >= np.deg2rad(lat_range[0])
    )
    boundary_vertices = np.logical_and(
        boundary_vertices,
        lat_vertex <= np.deg2rad(lat_range[1])
    )

    # convert from boolean mask to indices
    boundary_vertices = np.flatnonzero(boundary_vertices.values)

    mean_boundary_bsf = bsf_vertex.isel(nVertices=boundary_vertices).mean()

    bsf_shifted = bsf_vertex - mean_boundary_bsf

    return bsf_shifted
