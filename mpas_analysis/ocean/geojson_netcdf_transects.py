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
import json
from pathlib import Path

import numpy as np
import xarray as xr

from mpas_analysis.ocean.compute_transects_subtask import (
    ComputeTransectsSubtask,
    TransectsObservations
)
from mpas_analysis.ocean.plot_transect_subtask import PlotTransectSubtask
from mpas_analysis.shared import AnalysisTask


class GeojsonNetcdfTransects(AnalysisTask):
    """
    Plot model output at transects defined by lat/lon points in a geojson or
    NetCDF file
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, controlConfig=None):
        """
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        config : tranche.Tranche
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted
            as a transect

        controlconfig : tranche.Tranche, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        tags = ['climatology', 'transect', 'geojson', 'netcdf']

        # call the constructor from the base class (AnalysisTask)
        super().__init__(
            config=config, taskName='geojsonNetcdfTransects',
            componentName='ocean',
            tags=tags)

        sectionName = self.taskName

        geojsonOrNetcdfFiles = config.getexpression(sectionName,
                                                    'geojsonOrNetcdfFiles')
        if len(geojsonOrNetcdfFiles) == 0:
            return

        seasons = config.getexpression(sectionName, 'seasons')

        horizontalResolution = config.get(sectionName, 'horizontalResolution')

        verticalComparisonGridName = config.get(sectionName,
                                                'verticalComparisonGridName')

        if verticalComparisonGridName in ['mpas', 'obs']:
            verticalComparisonGrid = None
        else:
            verticalComparisonGrid = config.getexpression(
                sectionName, 'verticalComparisonGrid', allow_numpy=True)

        availableVariables = config.getexpression(
            sectionName, 'availableVariables')

        prefixes = config.getexpression(sectionName, 'variables')
        fields = [field for field in availableVariables
                  if field['prefix'] in prefixes]

        geojsonFileNames = {}
        netcdfFileNames = {}
        transectNames = []
        for fileName in geojsonOrNetcdfFiles:
            ext = Path(fileName).suffix
            if ext == '.nc':
                transectName = Path(fileName).stem
                netcdfFileNames[transectName] = fileName
            elif ext == '.geojson':
                with open(fileName) as filePointer:
                    jsonFile = json.load(filePointer)

                    for feature in jsonFile['features']:
                        if feature['geometry']['type'] != 'LineString':
                            continue
                        transectName = feature['properties']['name']

                        geojsonFileNames[transectName] = fileName
            else:
                raise ValueError(f'Unexptect file extension: {ext}')

            if transectName in transectNames:
                raise ValueError(f'Transect name {transectName} is repeated.')

            transectNames.append(transectName)

        variableList = [field['mpas'] for field in fields]

        computeGeojsonTransectsSubtask = None
        computeNetcdfTransectsSubtask = None

        if geojsonFileNames:
            transectCollectionName = 'geojson_transects'
            if horizontalResolution != 'obs' and \
                    horizontalResolution != 'mpas':
                transectCollectionName = \
                    f'{transectCollectionName}_{horizontalResolution}km'

            geojsonObservations = GeojsonTransectsObservations(
                config, geojsonFileNames, horizontalResolution,
                transectCollectionName)

            computeGeojsonTransectsSubtask = ComputeTransectsSubtask(
                mpasClimatologyTask=mpasClimatologyTask,
                parentTask=self,
                climatologyName='geojson',
                transectCollectionName=transectCollectionName,
                variableList=variableList,
                seasons=seasons,
                obsDatasets=geojsonObservations,
                verticalComparisonGridName=verticalComparisonGridName,
                verticalComparisonGrid=verticalComparisonGrid,
                subtaskName='remapGeojson')

        if netcdfFileNames:
            transectCollectionName = 'netcdf_transects'
            if horizontalResolution != 'obs' and \
                    horizontalResolution != 'mpas':
                transectCollectionName = \
                    f'{transectCollectionName}_{horizontalResolution}km'

            netcdfObservations = NetcdfTransectsObservations(
                config, netcdfFileNames, horizontalResolution,
                transectCollectionName)

            computeNetcdfTransectsSubtask = ComputeTransectsSubtask(
                mpasClimatologyTask=mpasClimatologyTask,
                parentTask=self,
                climatologyName='netcdf',
                transectCollectionName=transectCollectionName,
                variableList=variableList,
                seasons=seasons,
                obsDatasets=netcdfObservations,
                verticalComparisonGridName=verticalComparisonGridName,
                verticalComparisonGrid=verticalComparisonGrid,
                subtaskName='remapNetcdf')

        for field in fields:
            for transectName in geojsonFileNames:
                for season in seasons:
                    self._add_plot_subtasks(
                        field=field, season=season,
                        transectName=transectName,
                        computeSubtask=computeGeojsonTransectsSubtask,
                        galleryGroup='Geojson Transects',
                        groupLink='geojson',
                        controlConfig=controlConfig)

            for transectName in netcdfFileNames:
                for season in seasons:
                    self._add_plot_subtasks(
                        field=field, season=season,
                        transectName=transectName,
                        computeSubtask=computeNetcdfTransectsSubtask,
                        galleryGroup='NetCDF Transects',
                        groupLink='nctransects',
                        controlConfig=controlConfig)

    def _add_plot_subtasks(self, field, season, transectName, computeSubtask,
                           galleryGroup, groupLink, controlConfig):
        """
        Add a sbutask for plotting the given transect, field and season
        """
        config = self.config
        sectionName = self.taskName
        verticalBounds = config.getexpression(sectionName, 'verticalBounds')

        fieldPrefix = field['prefix']

        outFileLabel = fieldPrefix

        if controlConfig is None:
            refFieldName = None
            refTitleLabel = None
            diffTitleLabel = None
        else:
            refFieldName = field['mpas']
            controlRunName = controlConfig.get('runs', 'mainRunName')
            refTitleLabel = f'Control: {controlRunName}'
            diffTitleLabel = 'Main - Control'

        fieldPrefixUpper = fieldPrefix[0].upper() + fieldPrefix[1:]
        fieldNameInTitle = field['titleName']
        transectNameInTitle = transectName.replace('_', ' ')
        fieldNameInTitle = \
            f'{fieldNameInTitle} from {transectNameInTitle}'

        configSectionName = f'geojsonNetcdf{fieldPrefixUpper}Transects'

        # make a new subtask for this season and comparison grid
        subtask = PlotTransectSubtask(
            parentTask=self,
            season=season,
            transectName=transectName,
            fieldName=fieldPrefix,
            computeTransectsSubtask=computeSubtask,
            plotObs=False,
            controlConfig=controlConfig)

        subtask.set_plot_info(
            outFileLabel=outFileLabel,
            fieldNameInTitle=fieldNameInTitle,
            mpasFieldName=field['mpas'],
            refFieldName=refFieldName,
            refTitleLabel=refTitleLabel,
            diffTitleLabel=diffTitleLabel,
            unitsLabel=field['units'],
            imageCaption=fieldNameInTitle,
            galleryGroup=galleryGroup,
            groupSubtitle=None,
            groupLink=groupLink,
            galleryName=field['titleName'],
            configSectionName=configSectionName,
            verticalBounds=verticalBounds)

        self.add_subtask(subtask)


class GeojsonTransectsObservations(TransectsObservations):
    """
    A class for loading and manipulating geojson transects
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def build_observational_dataset(self, fileName, transectName):
        """
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Parameters
        ----------
        fileName : str
            observation file name

        transectName : str
            transect name

        Returns
        -------
        dsObs : ``xarray.Dataset``
            The observational dataset
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        with open(fileName) as filePointer:
            jsonFile = json.load(filePointer)

            for feature in jsonFile['features']:
                if feature['properties']['name'] != transectName:
                    continue
                assert feature['geometry']['type'] == 'LineString'

                coordinates = feature['geometry']['coordinates']
                lon, lat = zip(*coordinates)
                break

        dsObs = xr.Dataset()
        dsObs['lon'] = (('nPoints',), np.array(lon))
        dsObs.lon.attrs['units'] = 'degrees'
        dsObs['lat'] = (('nPoints',), np.array(lat))
        dsObs.lat.attrs['units'] = 'degrees'

        return dsObs


class NetcdfTransectsObservations(TransectsObservations):
    """
    A class for loading and manipulating netcdf transects
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def build_observational_dataset(self, fileName, transectName):
        """
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Parameters
        ----------
        fileName : str
            observation file name

        transectName : str
            transect name

        Returns
        -------
        dsObs : ``xarray.Dataset``
            The observational dataset
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        dsObs = xr.open_dataset(fileName)
        # drop all variables besides lon and lat
        dsObs = dsObs[['lon', 'lat']]

        return dsObs
