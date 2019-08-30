# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
'''
Analysis tasks for comparing Polar (and global) climatology maps against
WOA18 climatological data.
'''
# Authors
# -------
# Milena Veneziani

import xarray as xr
import numpy as np

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.ocean.remap_depth_slices_subtask import \
    RemapDepthSlicesSubtask
from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask

from mpas_analysis.shared.io.utility import build_obs_path

from mpas_analysis.shared.climatology import RemapObservedClimatologySubtask

from mpas_analysis.shared.grid import LatLonGridDescriptor

from mpas_analysis.shared.mpas_xarray import mpas_xarray


class ClimatologyMapWoa(AnalysisTask):  # {{{
    """
    An analysis task for comparison of polar and global temperature and
    salinity against WOA18 climatology fields
    """
    # Authors
    # -------
    # Milena Veneziani

    def __init__(self, config, mpasClimatologyTask,
                 controlConfig=None):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Milena Veneziani

        fields = \
            [{'prefix': 'temperature',
              'mpas': 'timeMonthly_avg_activeTracers_temperature',
              'units': r'$\degree$C',
              'titleName': 'Potential Temperature',
              'obsFieldName': 't_an'},
             {'prefix': 'salinity',
              'mpas': 'timeMonthly_avg_activeTracers_salinity',
              'units': r'PSU',
              'titleName': 'Salinity',
              'obsFieldName': 's_an'}]

        tags = ['climatology', 'horizontalMap', 'woa', 'publicObs'] + \
            [field['prefix'] for field in fields]

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapWoa, self).__init__(
            config=config, taskName='climatologyMapWoa',
            componentName='ocean',
            tags=tags)

        sectionName = self.taskName

        fieldList = config.getExpression(sectionName, 'fieldList')
        fields = [field for field in fields if field['prefix'] in fieldList]

        # read in what seasons we want to plot
        seasons = config.getExpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of seasons'.format(sectionName))

        comparisonGridNames = config.getExpression(sectionName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of comparison grids'.format(sectionName))

        depths = config.getExpression(sectionName, 'depths')

        if len(depths) == 0:
            raise ValueError('config section {} does not contain valid '
                             'list of depths'.format(sectionName))

        variableList = [field['mpas'] for field in fields]

        remapMpasSubtask = RemapDepthSlicesSubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='WOA',
            variableList=variableList,
            seasons=seasons,
            depths=depths,
            comparisonGridNames=comparisonGridNames,
            iselValues=None)

        for field in fields:
            fieldPrefix = field['prefix']
            upperFieldPrefix = fieldPrefix[0].upper() + fieldPrefix[1:]
            sectionName = '{}{}'.format(self.taskName, upperFieldPrefix)

            if controlConfig is None:

                refTitleLabel = 'WOA18 Climatology'

                observationsDirectory = build_obs_path(
                        config, 'ocean', 'woa18Subdirectory')

                refFieldName = field['obsFieldName']
                outFileLabel = '{}WOA18'.format(fieldPrefix)
                galleryName = 'WOA18 Climatology'
                diffTitleLabel = 'Model - Climatology'
                # Since we have a WOA18 annual climatology file and
                # another file containing the 12 WOA18 monthly climatologies,
                # do the remapping for each season separately
                for season in seasons:
                    if season == 'ANN':
                        obsFileName = \
                            '{}/woa18_decav_04_TS_ann_20190829.nc'.format(
                                observationsDirectory)
                    else:
                        obsFileName = \
                            '{}/woa18_decav_04_TS_mon_20190829.nc'.format(
                                observationsDirectory)

                    subtaskName = 'remapObservations{}_{}'.format(
                            upperFieldPrefix, season)

                    remapObsSubtask = RemapWoaClimatology(
                            parentTask=self, seasons=[season],
                            fileName=obsFileName,
                            outFilePrefix='{}WOA18'.format(refFieldName),
                            fieldName=refFieldName,
                            depths=depths,
                            comparisonGridNames=comparisonGridNames,
                            subtaskName=subtaskName)

                    self.add_subtask(remapObsSubtask)

            else:
                remapObsSubtask = None
                controlRunName = controlConfig.get('runs', 'mainRunName')
                galleryName = 'Control: {}'.format(controlRunName)
                refTitleLabel = galleryName

                refFieldName = field['mpas']
                outFileLabel = fieldPrefix
                diffTitleLabel = 'Main - Control'

            for comparisonGridName in comparisonGridNames:
                for season in seasons:
                    for depth in depths:

                        subtaskName = 'plot{}_{}_{}_depth_{}'.format(
                                                            upperFieldPrefix,
                                                            season,
                                                            comparisonGridName,
                                                            depth)

                        subtask = PlotClimatologyMapSubtask(
                            parentTask=self,
                            season=season,
                            comparisonGridName=comparisonGridName,
                            remapMpasClimatologySubtask=remapMpasSubtask,
                            remapObsClimatologySubtask=remapObsSubtask,
                            controlConfig=controlConfig,
                            depth=depth,
                            subtaskName=subtaskName)

                        subtask.set_plot_info(
                            outFileLabel=outFileLabel,
                            fieldNameInTitle=field['titleName'],
                            mpasFieldName=field['mpas'],
                            refFieldName=refFieldName,
                            refTitleLabel=refTitleLabel,
                            diffTitleLabel=diffTitleLabel,
                            unitsLabel=field['units'],
                            imageCaption=field['titleName'],
                            galleryGroup=field['titleName'],
                            groupSubtitle=None,
                            groupLink='{}Woa'.format(fieldPrefix),
                            galleryName=galleryName,
                            configSectionName='climatologyMapWoa{}'.format(
                                upperFieldPrefix))

                        self.add_subtask(subtask)
        # }}}

    # }}}


class RemapWoaClimatology(RemapObservedClimatologySubtask):
    # {{{
    """
    A subtask for reading and remapping WOA fields to the comparison grid
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani

    def __init__(self, parentTask, seasons, fileName, outFilePrefix,
                 fieldName, depths,
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

        depths : list of {None, float, 'top', 'bot'}
            A list of depths at which the climatology will be sliced in the
            vertical.

        comparisonGridNames : list of {'latlon', 'arctic', 'antarctic'},
            optional
            The name(s) of the comparison grid to use for remapping.

        subtaskName : str, optional
            The name of the subtask
        '''
        # Authors
        # -------
        # Xylar Asay-Davis, Milena Veneziani

        self.fieldName = fieldName
        self.depths = depths

        # call the constructor from the base class
        # (RemapObservedClimatologySubtask)
        super(RemapWoaClimatology, self).__init__(
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
        # Xylar Asay-Davis, Milena Veneziani

        # Load WOA18 climatological data
        dsObs = self.build_observational_dataset(fileName)

        # create a descriptor of the observation grid using Lat/Lon
        # coordinates
        obsDescriptor = LatLonGridDescriptor.read(ds=dsObs,
                                                  latVarName='lat',
                                                  lonVarName='lon')
        dsObs.close()
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
        # Xylar Asay-Davis, Milena Veneziani

        # Load WOA18 climatological data
        dsObs = xr.open_dataset(fileName)

        # Rename coordinates to be consistent with other datasets
        if 'month' in dsObs.coords:
            # need to add a dummy year coordinate
            dsObs = dsObs.rename({'month': 'Time'})
            dsObs.coords['month'] = dsObs['Time']
            dsObs.coords['year'] = ('Time', np.ones(dsObs.dims['Time'], int))

        slices = []
        field = dsObs[self.fieldName]
        for depth in self.depths:
            if depth == 'top':
                slices.append(field.sel(method='nearest', depth=0.).drop(
                    'depth'))
            else:
                slices.append(field.sel(method='nearest', depth=-depth).drop(
                    'depth'))

        depthNames = [str(depth) for depth in self.depths]
        field = xr.concat(slices, dim='depthSlice')

        dsObs = xr.Dataset(data_vars={self.fieldName: field},
                           coords={'depthSlice': depthNames})

        return dsObs  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
