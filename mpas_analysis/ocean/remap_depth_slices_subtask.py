# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2020 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2020 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2020 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
import xarray as xr
import numpy as np

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask

from mpas_analysis.ocean.utility import compute_zmid


class RemapDepthSlicesSubtask(RemapMpasClimatologySubtask):  # {{{
    """
    A task for creating and remapping climatologies of MPAS fields sliced
    at a given set of depths

    Attributes
    ----------
    depths : list of {None, float, 'top', 'bot'}
        A list of depths at which the climatology will be sliced in the
        vertical.

    maxLevelCell : xarray.DataArray
        The vertical index of the bottom cell in MPAS results

    verticalIndices : xarray.DataArray
        The vertical indices of slice to be plotted
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, mpasClimatologyTask, parentTask, climatologyName,
                 variableList, seasons, depths, comparisonGridNames=['latlon'],
                 iselValues=None):
        # {{{
        '''
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped

        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        climatologyName : str
            A name that describes the climatology (e.g. a short version of
            the important field(s) in the climatology) used to name the
            subdirectories for each stage of the climatology

        variableList : list of str
            A list of variable names in ``timeSeriesStatsMonthly`` to be
            included in the climatologies

        seasons : list of str, optional
            A list of seasons (keys in ``shared.constants.monthDictionary``)
            to be computed or ['none'] (not ``None``) if only monthly
            climatologies are needed.

        depths : list of {None, float, 'top', 'bot'}
            A list of depths at which the climatology will be sliced in the
            vertical.

        comparisonGridNames : list of {'latlon', 'antarctic'}, optional
            The name(s) of the comparison grid to use for remapping.

        iselValues : dict, optional
            A dictionary of dimensions and indices (or ``None``) used to
            extract a slice of the MPAS field(s).
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.depths = depths

        # call the constructor from the base class
        # (RemapMpasClimatologySubtask)
        super(RemapDepthSlicesSubtask, self).__init__(
            mpasClimatologyTask, parentTask, climatologyName, variableList,
            seasons, comparisonGridNames, iselValues)

    def run_task(self):  # {{{
        """
        Compute climatologies of T or S  from ACME/MPAS output

        This function has been overridden to load ``maxLevelCell`` from a
        restart file for later use in indexing bottom T and S.
        ``verticalIndex`` is also computed for later indexing of
        the model level. It then simply calls the run function from
        ClimatologyMapOcean.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, load the land-ice mask from the restart file
        ds = xr.open_dataset(self.restartFileName)
        ds = ds[['maxLevelCell', 'bottomDepth', 'layerThickness']]
        ds = ds.isel(Time=0)

        self.maxLevelCell = ds.maxLevelCell - 1

        depthNames = [str(depth) for depth in self.depths]

        zMid = compute_zmid(ds.bottomDepth, ds.maxLevelCell-1,
                            ds.layerThickness)

        nVertLevels = zMid.shape[1]
        zMid.coords['verticalIndex'] = \
            ('nVertLevels',
             np.arange(nVertLevels))

        zTop = zMid.isel(nVertLevels=0)
        # Each vertical layer has at most one non-NaN value so the "sum"
        # over the vertical is used to collapse the array in the vertical
        # dimension
        zBot = zMid.where(zMid.verticalIndex == self.maxLevelCell).sum(
            dim='nVertLevels')

        verticalIndices = np.zeros((len(self.depths), ds.dims['nCells']), int)

        mask = np.zeros(verticalIndices.shape, bool)

        for depthIndex, depth in enumerate(self.depths):
            depth = self.depths[depthIndex]
            if depth == 'top':
                # switch to zero-based index
                verticalIndices[depthIndex, :] = 0
                mask[depthIndex, :] = self.maxLevelCell.values >= 0
            elif depth == 'bot':
                # switch to zero-based index
                verticalIndices[depthIndex, :] = self.maxLevelCell.values
                mask[depthIndex, :] = self.maxLevelCell.values >= 0
            else:

                verticalIndex = np.abs(zMid - depth).argmin(dim='nVertLevels')

                verticalIndices[depthIndex, :] = verticalIndex.values
                mask[depthIndex, :] = np.logical_and(depth <= zTop,
                                                     depth >= zBot).values

        self.verticalIndices = \
            xr.DataArray.from_dict({'dims': ('depthSlice', 'nCells'),
                                    'coords': {'depthSlice':
                                               {'dims': ('depthSlice',),
                                                'data': depthNames}},
                                    'data': verticalIndices})
        self.verticalIndexMask = \
            xr.DataArray.from_dict({'dims': ('depthSlice', 'nCells'),
                                    'coords': {'depthSlice':
                                               {'dims': ('depthSlice',),
                                                'data': depthNames}},
                                    'data': mask})

        # then, call run from the base class (RemapMpasClimatologySubtask),
        # which will perform the main function of the task
        super(RemapDepthSlicesSubtask, self).run_task()

    def customize_masked_climatology(self, climatology, season):  # {{{
        """
        Uses ``verticalIndex`` to slice the 3D climatology field at each
        requested depth.  The resulting field has the depth appended to
        the variable name.

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

        if self.depths is None:
            return climatology

        climatology.coords['verticalIndex'] = \
            ('nVertLevels',
             np.arange(climatology.dims['nVertLevels']))

        depthNames = [str(depth) for depth in self.depths]

        climatology.coords['depthSlice'] = ('depthSlice', depthNames)

        for variableName in self.variableList:
            if 'nVertLevels' not in climatology[variableName].dims:
                continue

            # mask only the values with the right vertical index
            da = climatology[variableName].where(
                climatology.verticalIndex == self.verticalIndices)

            # Each vertical layer has at most one non-NaN value so the "sum"
            # over the vertical is used to collapse the array in the vertical
            # dimension
            climatology[variableName] = \
                da.sum(dim='nVertLevels').where(self.verticalIndexMask)

        climatology = climatology.drop_vars('verticalIndex')

        climatology = climatology.transpose('depthSlice', 'nCells')

        return climatology  # }}}

    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
