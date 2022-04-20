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
from collections import OrderedDict
import json
import xarray
import numpy

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.ocean.compute_transects_subtask import \
    ComputeTransectsSubtask, TransectsObservations

from mpas_analysis.ocean.plot_transect_subtask import PlotTransectSubtask


class GeojsonTransects(AnalysisTask):
    """
    Plot model output at transects defined by lat/lon points in a geojson file
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
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted
            as a transect

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        tags = ['climatology', 'transect', 'geojson']

        # call the constructor from the base class (AnalysisTask)
        super(GeojsonTransects, self).__init__(
            config=config, taskName='geojsonTransects',
            componentName='ocean',
            tags=tags)

        sectionName = self.taskName

        geojsonFiles = config.getexpression(sectionName, 'geojsonFiles')
        if len(geojsonFiles) == 0:
            return

        seasons = config.getexpression(sectionName, 'seasons')

        horizontalResolution = config.get(sectionName, 'horizontalResolution')

        verticalComparisonGridName = config.get(sectionName,
                                                'verticalComparisonGridName')

        if verticalComparisonGridName in ['mpas', 'obs']:
            verticalComparisonGrid = None
        else:
            verticalComparisonGrid = config.getexpression(
                sectionName, 'verticalComparisonGrid', use_numpyfunc=True)

        verticalBounds = config.getexpression(sectionName, 'verticalBounds')

        fields = config.getexpression(sectionName, 'fields')

        obsFileNames = OrderedDict()
        for fileName in geojsonFiles:
            with open(fileName) as filePointer:
                jsonFile = json.load(filePointer)

                for feature in jsonFile['features']:
                    if feature['geometry']['type'] != 'LineString':
                        continue
                    transectName = feature['properties']['name']

                    obsFileNames[transectName] = fileName

        transectCollectionName = 'geojson_transects'
        if horizontalResolution != 'obs':
            transectCollectionName = '{}_{}km'.format(transectCollectionName,
                                                      horizontalResolution)

        transectsObservations = GeojsonTransectsObservations(
            config, obsFileNames, horizontalResolution,
            transectCollectionName)

        computeTransectsSubtask = ComputeTransectsSubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='geojson',
            transectCollectionName=transectCollectionName,
            variableList=[field['mpas'] for field in fields],
            seasons=seasons,
            obsDatasets=transectsObservations,
            verticalComparisonGridName=verticalComparisonGridName,
            verticalComparisonGrid=verticalComparisonGrid)

        plotObs = False
        if controlConfig is None:

            refTitleLabel = None

            diffTitleLabel = None

        else:
            controlRunName = controlConfig.get('runs', 'mainRunName')
            refTitleLabel = 'Control: {}'.format(controlRunName)

            diffTitleLabel = 'Main - Control'

        for field in fields:
            fieldPrefix = field['prefix']
            for transectName in obsFileNames:
                for season in seasons:
                    outFileLabel = fieldPrefix
                    if controlConfig is None:
                        refFieldName = None
                    else:
                        refFieldName = field['mpas']

                    fieldPrefixUpper = fieldPrefix[0].upper() + fieldPrefix[1:]
                    fieldNameInTytle = '{} from {}'.format(
                        field['titleName'],
                        transectName.replace('_', ' '))

                    # make a new subtask for this season and comparison grid
                    subtask = PlotTransectSubtask(self, season, transectName,
                                                  fieldPrefix,
                                                  computeTransectsSubtask,
                                                  plotObs, controlConfig)

                    subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle=fieldNameInTytle,
                        mpasFieldName=field['mpas'],
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=field['units'],
                        imageCaption=fieldNameInTytle,
                        galleryGroup='Geojson Transects',
                        groupSubtitle=None,
                        groupLink='geojson',
                        galleryName=field['titleName'],
                        configSectionName='geojson{}Transects'.format(
                            fieldPrefixUpper),
                        verticalBounds=verticalBounds)

                    self.add_subtask(subtask)


class GeojsonTransectsObservations(TransectsObservations):
    """
    A class for loading and manipulating geojson transects

    Attributes
    ----------

    obsDatasets : OrderedDict
        A dictionary of observational datasets
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
                assert(feature['geometry']['type'] == 'LineString')

                coordinates = feature['geometry']['coordinates']
                lon, lat = zip(*coordinates)
                break

        dsObs = xarray.Dataset()
        dsObs['lon'] = (('nPoints',), numpy.array(lon))
        dsObs.lon.attrs['units'] = 'degrees'
        dsObs['lat'] = (('nPoints',), numpy.array(lat))
        dsObs.lat.attrs['units'] = 'degrees'

        return dsObs
