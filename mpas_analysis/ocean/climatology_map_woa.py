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
"""
Analysis tasks for comparing Polar (and global) climatology maps against
WOA23 climatological data.
"""
# Authors
# -------
# Milena Veneziani

import xarray as xr
import numpy as np
from pyremap import LatLonGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.ocean.remap_depth_slices_subtask import \
    RemapDepthSlicesSubtask
from mpas_analysis.shared.plot import PlotClimatologyMapSubtask

from mpas_analysis.shared.io.utility import build_obs_path

from mpas_analysis.shared.climatology import RemapObservedClimatologySubtask


class ClimatologyMapWoa(AnalysisTask):
    """
    An analysis task for comparison of polar and global temperature and
    salinity against WOA23 climatology fields
    """
    # Authors
    # -------
    # Milena Veneziani

    def __init__(self, config, mpasClimatologyTask,
                 controlConfig=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        controlconfig : mpas_tools.config.MpasConfigParser, optional
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
              'obsFieldName': 'pt_an'},
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

        fieldList = config.getexpression(sectionName, 'fieldList')
        fields = [field for field in fields if field['prefix'] in fieldList]

        # read in what seasons we want to plot
        seasons = config.getexpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of seasons'.format(sectionName))

        comparisonGridNames = config.getexpression(sectionName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of comparison grids'.format(sectionName))

        depths = config.getexpression(sectionName, 'depths')

        if len(depths) == 0:
            raise ValueError('config section {} does not contain valid '
                             'list of depths'.format(sectionName))

        variableList = [field['mpas'] for field in fields]

        shallowVsDeepColormapDepth = config.getfloat(
            sectionName, 'shallowVsDeepColormapDepth')

        shallow = []
        for depth in depths:
            if depth == 'top':
                shallow.append(True)
            elif depth == 'bot':
                shallow.append(False)
            else:
                shallow.append(depth >= shallowVsDeepColormapDepth)

        remapMpasSubtask = RemapDepthSlicesSubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='WOA',
            variableList=variableList,
            seasons=seasons,
            depths=depths,
            comparisonGridNames=comparisonGridNames,
            iselValues=None)

        remapObsSubtask = None

        if controlConfig is None:
            observationsDirectory = build_obs_path(
                    config, 'ocean', 'woa23Subdirectory')

            refFieldNames = [field['obsFieldName'] for field in fields]

            obsFileName = f'{observationsDirectory}/' \
                          f'woa23_decav_04_pt_s_mon_ann.20241101.nc'

            remapObsSubtask = RemapWoaClimatology(
                parentTask=self, seasons=seasons,
                fileName=obsFileName,
                outFilePrefix='woa23',
                fieldNames=refFieldNames,
                depths=depths,
                comparisonGridNames=comparisonGridNames,
                subtaskName='remapObservationsMon')

            self.add_subtask(remapObsSubtask)

        for field in fields:
            fieldPrefix = field['prefix']
            upperFieldPrefix = fieldPrefix[0].upper() + fieldPrefix[1:]
            sectionName = f'{self.taskName}{upperFieldPrefix}'

            if controlConfig is None:

                refTitleLabel = 'WOA23 Climatology'

                refFieldName = field['obsFieldName']
                outFileLabel = f'{fieldPrefix}WOA23'
                galleryName = 'WOA23 Climatology'
                diffTitleLabel = 'Model - Climatology'

            else:
                controlRunName = controlConfig.get('runs', 'mainRunName')
                galleryName = f'Control: {controlRunName}'
                refTitleLabel = galleryName

                refFieldName = field['mpas']
                outFileLabel = f'{fieldPrefix}WOA23'
                diffTitleLabel = 'Main - Control'

            for comparisonGridName in comparisonGridNames:
                for depthIndex, depth in enumerate(depths):
                    for season in seasons:

                        subtaskName = f'plot{upperFieldPrefix}_{season}_' \
                                      f'{comparisonGridName}_depth_{depth}'

                        subtask = PlotClimatologyMapSubtask(
                            parentTask=self,
                            season=season,
                            comparisonGridName=comparisonGridName,
                            remapMpasClimatologySubtask=remapMpasSubtask,
                            remapObsClimatologySubtask=remapObsSubtask,
                            controlConfig=controlConfig,
                            depth=depth,
                            subtaskName=subtaskName)

                        configSectionName = \
                            f'climatologyMapWoa{upperFieldPrefix}'

                        # if available, use a separate color map for shallow
                        # and deep
                        if depth is not None:
                            if shallow[depthIndex]:
                                suffix = 'Shallow'
                            else:
                                suffix = 'Deep'
                            testSectionName = '{}{}'.format(configSectionName,
                                                            suffix)
                            if config.has_section(testSectionName):
                                configSectionName = testSectionName

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
                            groupLink=f'{fieldPrefix}_woa',
                            galleryName=galleryName,
                            configSectionName=configSectionName)

                        self.add_subtask(subtask)


class RemapWoaClimatology(RemapObservedClimatologySubtask):

    """
    A subtask for reading and remapping WOA fields to the comparison grid
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani

    def __init__(self, parentTask, seasons, fileName, outFilePrefix,
                 fieldNames, depths,
                 comparisonGridNames=['latlon'],
                 subtaskName='remapObservations'):

        """
        An analysis task for remapping WOA fields (either annual or monthly
        mean) to the comparison grid(s), depths and seasons provided

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent (main) task for this subtask

        seasons : list of str
           A list of seasons (keys in ``constants.monthDictionary``) over
           which the climatology should be computed.

        fileName : str
            The name of the observation file

        outFilePrefix : str
            The prefix in front of output files and mapping files, typically
            the name of the field being remapped

        fieldNames : list of str
            The names of the 3D fields to remap

        depths : list of {None, float, 'top'}
            A list of depths at which the climatology will be sliced in the
            vertical.

        comparisonGridNames : list of {'latlon', 'arctic', 'antarctic'},
            optional
            The name(s) of the comparison grid to use for remapping.

        subtaskName : str, optional
            The name of the subtask
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Milena Veneziani

        self.fieldNames = fieldNames
        self.depths = depths

        # call the constructor from the base class
        # (RemapObservedClimatologySubtask)
        super(RemapWoaClimatology, self).__init__(
            parentTask, seasons, fileName, outFilePrefix,
            comparisonGridNames, subtaskName)

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
        # Xylar Asay-Davis, Milena Veneziani

        # Load WOA23 climatological data
        dsObs = self.build_observational_dataset(fileName)

        # create a descriptor of the observation grid using Lat/Lon
        # coordinates
        obsDescriptor = LatLonGridDescriptor.read(ds=dsObs,
                                                  latVarName='lat',
                                                  lonVarName='lon')
        dsObs.close()
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
        # Xylar Asay-Davis, Milena Veneziani

        # Load WOA23 climatological data
        dsObs = xr.open_dataset(fileName)

        # must have capital Time
        dsObs = dsObs.rename({'time': 'Time'})
        # make sure month is a coord
        dsObs = dsObs.set_coords('month')
        # add a dummy year to the dataset
        dsObs.coords['year'] = ('Time', np.ones(dsObs.sizes['Time'], int))

        data_vars = {}
        for fieldName in self.fieldNames:
            slices = []
            field = dsObs[fieldName]
            for depth in self.depths:
                if depth == 'top':
                    slices.append(field.sel(method='nearest',
                                            depth=0.).drop_vars('depth'))
                else:
                    slices.append(field.sel(method='nearest',
                                            depth=-depth).drop_vars('depth'))

            depthNames = [str(depth) for depth in self.depths]
            field = xr.concat(slices, dim='depthSlice')
            data_vars[fieldName] = field

        coords = {'depthSlice': depthNames}

        dsObs = xr.Dataset(data_vars=data_vars, coords=coords)

        return dsObs
