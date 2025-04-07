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
Analysis tasks for comparing Global climatology maps against Argo data.
"""
# Authors
# -------
# Luke Van Roekel

import xarray as xr
import numpy as np

from pyremap import LatLonGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.ocean.remap_depth_slices_subtask import \
    RemapDepthSlicesSubtask

from mpas_analysis.shared.plot import PlotClimatologyMapSubtask

from mpas_analysis.shared.io.utility import build_obs_path

from mpas_analysis.shared.climatology import RemapObservedClimatologySubtask


class ClimatologyMapArgoTemperature(AnalysisTask):
    """
    An analysis task for comparison of potential temperature against Argo
    observations
    """

    # Authors
    # -------
    # Luke Van Roekel, Xylar Asay-Davis

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
        # Luke Van Roekel, Xylar Asay-Davis

        fieldName = 'temperatureArgo'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapArgoTemperature, self).__init__(
            config=config, taskName='climatologyMapArgoTemperature',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', 'argo', 'temperature',
                  'publicObs'])

        sectionName = self.taskName

        mpasFieldName = 'timeMonthly_avg_activeTracers_temperature'
        iselValues = None

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

        # the variable 'timeMonthly_avg_landIceFreshwaterFlux' will be added to
        # mpasClimatologyTask along with the seasons.
        remapClimatologySubtask = RemapDepthSlicesSubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=fieldName,
            variableList=[mpasFieldName],
            seasons=seasons,
            depths=depths,
            comparisonGridNames=comparisonGridNames,
            iselValues=iselValues)

        if controlConfig is None:

            refTitleLabel = 'Roemmich-Gilson Argo Climatology: Potential ' \
                            'Temperature'

            observationsDirectory = build_obs_path(
                config, 'ocean', 'argoSubdirectory')

            obsFileName = \
                '{}/ArgoClimatology_TS_20180710.nc'.format(
                    observationsDirectory)
            refFieldName = 'theta'
            outFileLabel = 'tempArgo'
            galleryName = 'Roemmich-Gilson Climatology: Argo'
            diffTitleLabel = 'Model - Argo'

            remapObservationsSubtask = RemapArgoClimatology(
                parentTask=self, seasons=seasons, fileName=obsFileName,
                outFilePrefix='{}Argo'.format(refFieldName),
                fieldName=refFieldName,
                depths=depths,
                comparisonGridNames=comparisonGridNames)

            self.add_subtask(remapObservationsSubtask)

        else:
            remapObservationsSubtask = None
            controlRunName = controlConfig.get('runs', 'mainRunName')
            galleryName = 'Control: {}'.format(controlRunName)
            refTitleLabel = galleryName

            refFieldName = mpasFieldName
            outFileLabel = 'temp'
            diffTitleLabel = 'Main - Control'

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                for depth in depths:
                    subtask = PlotClimatologyMapSubtask(
                        parentTask=self,
                        season=season,
                        comparisonGridName=comparisonGridName,
                        remapMpasClimatologySubtask=remapClimatologySubtask,
                        remapObsClimatologySubtask=remapObservationsSubtask,
                        controlConfig=controlConfig,
                        depth=depth)

                    subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle='Potential Temperature',
                        mpasFieldName=mpasFieldName,
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=r'$\degree$C',
                        imageCaption='Model potential temperature compared '
                                     'with Argo observations',
                        galleryGroup='Argo Potential Temperature',
                        groupSubtitle=None,
                        groupLink='tempArgo',
                        galleryName=galleryName)

                    self.add_subtask(subtask)


class ClimatologyMapArgoSalinity(AnalysisTask):
    """
    An analysis task for comparison of global salinity against Argo
    observations
    """

    # Authors
    # -------
    # Xylar Asay-Davis, Luke Van Roekel

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
        # Xylar Asay-Davis, Luke Van Roekel

        fieldName = 'salinityArgo'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapArgoSalinity, self).__init__(
            config=config, taskName='climatologyMapArgoSalinity',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', 'argo', 'salinity'])

        sectionName = self.taskName

        mpasFieldName = 'timeMonthly_avg_activeTracers_salinity'
        iselValues = None

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

        # the variable 'timeMonthly_avg_landIceFreshwaterFlux' will be added to
        # mpasClimatologyTask along with the seasons.
        remapClimatologySubtask = RemapDepthSlicesSubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=fieldName,
            variableList=[mpasFieldName],
            seasons=seasons,
            depths=depths,
            comparisonGridNames=comparisonGridNames,
            iselValues=iselValues)

        if controlConfig is None:

            refTitleLabel = 'Roemmich-Gilson Argo Climatology: Salinity'

            observationsDirectory = build_obs_path(
                config, 'ocean', 'argoSubdirectory')

            obsFileName = \
                '{}/ArgoClimatology_TS_20180710.nc'.format(
                    observationsDirectory)
            refFieldName = 'salinity'
            outFileLabel = 'salinArgo'
            galleryName = 'Roemmich-Gilson Climatology: Argo'
            diffTitleLabel = 'Model - Argo'

            remapObservationsSubtask = RemapArgoClimatology(
                parentTask=self, seasons=seasons, fileName=obsFileName,
                outFilePrefix='{}Argo'.format(refFieldName),
                fieldName=refFieldName,
                depths=depths,
                comparisonGridNames=comparisonGridNames)

            self.add_subtask(remapObservationsSubtask)

        else:
            remapObservationsSubtask = None
            controlRunName = controlConfig.get('runs', 'mainRunName')
            galleryName = None
            refTitleLabel = 'Control: {}'.format(controlRunName)

            refFieldName = mpasFieldName
            outFileLabel = 'salin'
            diffTitleLabel = 'Main - Control'

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                for depth in depths:
                    subtask = PlotClimatologyMapSubtask(
                        parentTask=self,
                        season=season,
                        comparisonGridName=comparisonGridName,
                        remapMpasClimatologySubtask=remapClimatologySubtask,
                        remapObsClimatologySubtask=remapObservationsSubtask,
                        controlConfig=controlConfig,
                        depth=depth)

                    subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle='Salinity',
                        mpasFieldName=mpasFieldName,
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=r'PSU',
                        imageCaption='Model Salinity compared with Argo '
                                     'observations',
                        galleryGroup='Argo Salinity',
                        groupSubtitle=None,
                        groupLink='salinArgo',
                        galleryName=galleryName)

                    self.add_subtask(subtask)


class RemapArgoClimatology(RemapObservedClimatologySubtask):
    """
    A subtask for reading and remapping SOSE fields to the comparison grid
    """

    # Authors
    # -------
    # Xylar Asay-Davis, Luke Van Roekel

    def __init__(self, parentTask, seasons, fileName, outFilePrefix,
                 fieldName, depths,
                 comparisonGridNames=['latlon'],
                 subtaskName='remapObservations'):

        """
        Construct one analysis subtask for each plot (i.e. each season and
        comparison grid) and a subtask for computing climatologies.

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

        fieldName : str
            The name of the 3D field to remap

        depths : list of {None, float, 'top'}
            A list of depths at which the climatology will be sliced in the
            vertical.

        comparisonGridNames : list of {'latlon', 'antarctic'}, optional
            The name(s) of the comparison grid to use for remapping.

        subtaskName : str, optional
            The name of the subtask
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Luke Van Roekel

        self.fieldName = fieldName
        self.depths = depths

        # call the constructor from the base class
        # (RemapObservedClimatologySubtask)
        super(RemapArgoClimatology, self).__init__(
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
        # Xylar Asay-Davis, Luke Van Roekel

        # Load Argo observational Data
        dsObs = self.build_observational_dataset(fileName)

        # create a descriptor of the observation grid using Lat/Lon
        # coordinates
        obsDescriptor = LatLonGridDescriptor.read(ds=dsObs,
                                                  lat_var_name='latCoord',
                                                  lon_var_name='lonCoord')
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
        # Xylar Asay-Davis, Luke Van Roekel

        # Load Argo observational data
        dsObs = xr.open_dataset(fileName)

        # Rename coordinates to be consistent with other datasets
        dsObs = dsObs.rename({'month': 'calmonth', 'LATITUDE': 'latCoord',
                              'LONGITUDE': 'lonCoord', 'DEPTH': 'depth'})
        dsObs.coords['LATITUDE'] = dsObs['latCoord']
        dsObs.coords['LONGITUDE'] = dsObs['lonCoord']
        dsObs.coords['DEPTH'] = dsObs['depth']
        dsObs.coords['month'] = ('Time', np.array(dsObs['calmonth'], int))

        # no meaningful year since this is already a climatology
        dsObs.coords['year'] = ('Time', np.ones(dsObs.sizes['Time'], int))
        dsObs = dsObs[[self.fieldName, 'month']]

        slices = []
        field = dsObs[self.fieldName]
        for depth in self.depths:
            if depth == 'top':
                slices.append(field.sel(method='nearest', depth=0.).drop_vars(
                    'depth'))
            else:
                slices.append(
                    field.sel(method='nearest', depth=depth).drop_vars('depth'))

        depthNames = [str(depth) for depth in self.depths]
        field = xr.concat(slices, dim='depthSlice')

        dsObs = xr.Dataset(data_vars={self.fieldName: field},
                           coords={'depthSlice': depthNames})

        return dsObs
