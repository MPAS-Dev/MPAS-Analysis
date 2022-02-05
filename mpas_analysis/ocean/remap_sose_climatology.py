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
'''
Analysis tasks for comparing Antarctic climatology maps against observations
and reanalysis data.
'''
# Authors
# -------
# Xylar Asay-Davis

import xarray as xr

from pyremap import ProjectionGridDescriptor

from mpas_analysis.shared.climatology import RemapObservedClimatologySubtask
from mpas_analysis.shared.projection import get_pyproj_projection


class RemapSoseClimatology(RemapObservedClimatologySubtask):
    # {{{
    """
    A subtask for reading and remapping SOSE fields to the comparison grid
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, seasons, fileName, outFilePrefix,
                 fieldName, botFieldName=None, depths=None,
                 comparisonGridNames=['latlon'],
                 subtaskName='remapObservations'):
        # {{{
        '''
        Construct one analysis subtask for each plot (i.e. each season and
        comparison grid) and a subtask for computing climatologies.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent (master) task for this subtask

        seasons : list of str
           A list of seasons (keys in ``constants.monthDictionary``) over
           which the climatology should be computed.

        fileName : str
            The name of the observation file

        outFilePrefix : str
            The prefix in front of output files and mapping files, typically
            the name of the field being remapped

        fieldName : str
            The name of the 3D field to remap

        botFieldName : str, optional
            The name of the same field as ``fieldName`` but sampled at the
            sea floor

        depths : list of {None, float, 'top', 'bot'}, optional
            A list of depths at which the climatology will be sliced in the
            vertical.

        comparisonGridNames : list of {'latlon', 'antarctic'}, optional
            The name(s) of the comparison grid to use for remapping.

        subtaskName : str, optional
            The name of the subtask
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.fieldName = fieldName
        self.botFieldName = botFieldName
        self.depths = depths

        if depths is not None and 'bot' in depths:
            assert(botFieldName is not None)

        # call the constructor from the base class
        # (RemapObservedClimatologySubtask)
        super(RemapSoseClimatology, self).__init__(
            parentTask, seasons, fileName, outFilePrefix,
            comparisonGridNames, subtaskName)
        # }}}

    def get_observation_descriptor(self, fileName):  # {{{
        '''
        get a MeshDescriptor for the observation grid

        Parameters
        ----------
        fileName : str
            observation file name describing the source grid

        Returns
        -------
        obsDescriptor : ``MeshDescriptor``
            The descriptor for the observation grid
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # create a descriptor of the observation grid using the x/y polar
        # stereographic coordinates
        projection = get_pyproj_projection(comparison_grid_name='antarctic')
        obsDescriptor = ProjectionGridDescriptor.read(
            projection, fileName=fileName, xVarName='x', yVarName='y')
        return obsDescriptor  # }}}

    def build_observational_dataset(self, fileName):  # {{{
        '''
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
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # Load MLD observational data
        dsObs = xr.open_dataset(fileName)

        varList = [self.fieldName, 'month', 'year']

        if self.botFieldName is not None:
            varList.append(self.botFieldName)
        dsObs = dsObs[varList]

        if self.depths is not None:
            field = dsObs[self.fieldName]
            slices = []
            for depth in self.depths:
                if depth == 'top':
                    slices.append(field.sel(method='nearest', z=0.).drop_vars(
                        'z'))
                elif depth == 'bot':
                    slices.append(dsObs[self.botFieldName])
                else:
                    level = field.sel(method='nearest', z=depth).drop_vars(
                        'z')
                    slices.append(level)

            depthNames = [str(depth) for depth in self.depths]
            field = xr.concat(slices, dim='depthSlice')

            dsObs = xr.Dataset(data_vars={self.fieldName: field},
                               coords={'depthSlice': depthNames})

        return dsObs  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
