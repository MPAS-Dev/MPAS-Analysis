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
import numpy as np

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask

from mpas_analysis.ocean.utility import (
    compute_zinterface,
    compute_zmid
)


class RemapDepthSlicesSubtask(RemapMpasClimatologySubtask):
    """
    A task for creating and remapping climatologies of MPAS fields sliced
    at a given set of depths

    Attributes
    ----------
    depths : list of {None, float, 'top', 'bot'}
        A list of depths at which the climatology will be sliced in the
        vertical.

    dsSlice : xarray.Dataset
        A dataset containing information needed to index variables at the
        designated depths
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, mpasClimatologyTask, parentTask, climatologyName,
                 variableList, seasons, depths, comparisonGridNames=['latlon'],
                 iselValues=None, subtaskName='remapDepthSlices'):

        """
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

        subtaskName : str, optional
            The name of the subtask
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.depths = depths
        self.dsSlice = xr.Dataset()

        # call the constructor from the base class
        # (RemapMpasClimatologySubtask)
        super().__init__(
            mpasClimatologyTask, parentTask, climatologyName, variableList,
            seasons, comparisonGridNames, iselValues,
            subtaskName=subtaskName)

    def run_task(self):
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

        # first, load the land-ice mask from the mesh file
        ds = xr.open_dataset(self.meshFilename)
        ds = ds[['maxLevelCell', 'bottomDepth', 'layerThickness']]
        ds = ds.isel(Time=0)

        depthNames = [str(depth) for depth in self.depths]

        bottomDepth = ds.bottomDepth
        layerThickness = ds.layerThickness
        maxLevelCell = ds.maxLevelCell - 1
        self.dsSlice['maxLevelCell'] = maxLevelCell

        zMid = compute_zmid(bottomDepth, maxLevelCell, layerThickness)

        zInterface = compute_zinterface(
            bottomDepth, maxLevelCell, layerThickness)

        horizontalMask = maxLevelCell >= 0

        nVertLevels = ds.sizes['nVertLevels']
        zMid.coords['verticalIndex'] = \
            ('nVertLevels',
             np.arange(nVertLevels))

        nVertLevelsP1 = zInterface.sizes['nVertLevelsP1']
        zInterface.coords['verticalIndex'] = \
            ('nVertLevelsP1',
             np.arange(nVertLevelsP1))

        zLevelTop = zMid.isel(nVertLevels=0)
        # Each vertical layer has at most one non-NaN value so the "sum"
        # over the vertical is used to collapse the array in the vertical
        # dimension
        zLevelBot = zMid.where(zMid.verticalIndex == maxLevelCell).sum(
            dim='nVertLevels')

        zInterfaceTop = zInterface.isel(nVertLevelsP1=0)
        zInterfaceBot = zInterface.where(
            zInterface.verticalIndex == maxLevelCell + 1).sum(
                dim='nVertLevelsP1')

        levelIndices = np.zeros((len(self.depths), ds.sizes['nCells']), int)
        levelMask = np.zeros(levelIndices.shape, bool)
        interfaceIndices = np.zeros(levelIndices.shape, int)
        interfaceMask = np.zeros(levelIndices.shape, bool)

        for depthIndex, depth in enumerate(self.depths):
            depth = self.depths[depthIndex]
            if depth == 'top':
                levelIndices[depthIndex, :] = 0
                levelMask[depthIndex, :] = horizontalMask.values
                interfaceIndices[depthIndex, :] = 0
                interfaceMask[depthIndex, :] = horizontalMask.values
            elif depth == 'bot':
                # switch to zero-based index
                levelIndices[depthIndex, :] = maxLevelCell.values
                levelMask[depthIndex, :] = horizontalMask.values
                interfaceIndices[depthIndex, :] = maxLevelCell.values + 1
                interfaceMask[depthIndex, :] = horizontalMask.values
            else:
                levelDiff = np.abs(zMid - depth).where(horizontalMask,
                                                       drop=True)
                levelIndex = levelDiff.argmin(dim='nVertLevels')

                levelIndices[depthIndex, horizontalMask.values] = \
                    levelIndex.values
                levelMask[depthIndex, :] = np.logical_and(
                    depth <= zLevelTop, depth >= zLevelBot).values

                interfaceDiff = np.abs(zInterface - depth).where(
                    horizontalMask, drop=True)
                interfaceIndex = interfaceDiff.argmin(dim='nVertLevelsP1')

                interfaceIndices[depthIndex, horizontalMask.values] = \
                    interfaceIndex.values
                interfaceMask[depthIndex, :] = np.logical_and(
                    depth <= zInterfaceTop, depth >= zInterfaceBot).values

        self.dsSlice.coords['depthSlice'] = ('depthSlice', depthNames)

        self.dsSlice['levelIndices'] = (('depthSlice', 'nCells'),
                                        levelIndices)
        self.dsSlice['levelIndexMask'] = (('depthSlice', 'nCells'),
                                          levelMask)
        self.dsSlice['interfaceIndices'] = (('depthSlice', 'nCells'),
                                            interfaceIndices)
        self.dsSlice['interfaceIndexMask'] = (('depthSlice', 'nCells'),
                                              interfaceMask)

        # then, call run from the base class (RemapMpasClimatologySubtask),
        # which will perform the main function of the task
        super().run_task()

    def customize_masked_climatology(self, climatology, season):
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

        if 'nVertLevels' in climatology.dims:
            climatology.coords['levelIndex'] = \
                ('nVertLevels',
                 np.arange(climatology.sizes['nVertLevels']))
        if 'nVertLevelsP1' in climatology.dims:
            climatology.coords['interfaceIndex'] = \
                ('nVertLevelsP1',
                 np.arange(climatology.sizes['nVertLevelsP1']))

        depthNames = [str(depth) for depth in self.depths]

        climatology.coords['depthSlice'] = ('depthSlice', depthNames)

        levelIndices = self.dsSlice.levelIndices
        levelIndexMask = self.dsSlice.levelIndexMask
        interfaceIndices = self.dsSlice.interfaceIndices
        interfaceIndexMask = self.dsSlice.interfaceIndexMask

        # iterate over all variables since some new ones may have been
        # added by a subclass
        for variableName in climatology.data_vars:
            if 'nVertLevels' in climatology[variableName].dims:
                # mask only the values with the right vertical index
                da = climatology[variableName].where(
                    climatology.levelIndex == levelIndices)

                # Each vertical layer has at most one non-NaN value so the
                # "sum" over the vertical is used to collapse the array in the
                # vertical dimension
                climatology[variableName] = \
                    da.sum(dim='nVertLevels').where(levelIndexMask)
            elif 'nVertLevelsP1' in climatology[variableName].dims:
                da = climatology[variableName].where(
                    climatology.interfaceIndex == interfaceIndices)

                climatology[variableName] = \
                    da.sum(dim='nVertLevelsP1').where(interfaceIndexMask)

        if 'levelIndex' in climatology.coords:
            climatology = climatology.drop_vars('levelIndex')
        if 'interfaceIndex' in climatology.coords:
            climatology = climatology.drop_vars('interfaceIndex')

        climatology = climatology.transpose('depthSlice', 'nCells')

        return climatology
