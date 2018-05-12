from __future__ import absolute_import, division, print_function, \
    unicode_literals

from collections import OrderedDict
import json
import xarray
import numpy

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.ocean.compute_transects_subtask import \
    ComputeTransectsSubtask, TransectsObservations

from mpas_analysis.ocean.plot_transect_subtask import PlotTransectSubtask


class GeojsonTransects(AnalysisTask):  # {{{
    """
    Plot model output at transects defined by lat/lon points in a geojson file

    Attributes
    ----------

    Authors
    -------
    Xylar Asay-Davis
    """

    def __init__(self, config, mpasClimatologyTask, refConfig=None):

        # {{{
        '''
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted
            as a transect

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any)

        Authors
        -------
        Xylar Asay-Davis
        '''

        tags = ['climatology', 'transect', 'geojson']

        # call the constructor from the base class (AnalysisTask)
        super(GeojsonTransects, self).__init__(
                config=config, taskName='geojsonTransects',
                componentName='ocean',
                tags=tags)

        sectionName = self.taskName

        seasons = config.getExpression(sectionName, 'seasons')

        horizontalResolution = config.get(sectionName, 'horizontalResolution')

        verticalComparisonGridName = config.get(sectionName,
                                                'verticalComparisonGridName')

        if verticalComparisonGridName in ['mpas', 'obs']:
            verticalComparisonGrid = None
        else:
            verticalComparisonGrid = config.getExpression(
                    sectionName, 'verticalComparisonGrid', usenumpyfunc=True)

        fields = config.getExpression(sectionName, 'fields')

        geojsonFiles = config.getExpression(sectionName, 'geojsonFiles')
        if len(geojsonFiles):
            print('Waring: No geojson files to plot in {}'.format(
                    self.taskName))

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
                config, obsFileNames,  horizontalResolution,
                transectCollectionName)

        computeTransectsSubtask = ComputeTransectsSubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='geojson',
            transectCollectionName=transectCollectionName,
            variableList=[field['mpas'] for field in fields.values()],
            seasons=seasons,
            obsDatasets=transectsObservations,
            verticalComparisonGridName=verticalComparisonGridName,
            verticalComparisonGrid=verticalComparisonGrid)

        plotObs = False
        if refConfig is None:

            refTitleLabel = None

            diffTitleLabel = None

        else:
            refRunName = refConfig.get('runs', 'mainRunName')
            refTitleLabel = 'Ref: {}'.format(refRunName)

            diffTitleLabel = 'Main - Reference'

        for fieldName in fields:
            for transectName in obsFileNames:
                for season in seasons:
                    outFileLabel = fieldName
                    if refConfig is None:
                        refFieldName = None
                    else:
                        refFieldName = fields[fieldName]['mpas']

                    fieldNameUpper = fieldName[0].upper() + fieldName[1:]
                    fieldNameInTytle = '{} from {}'.format(
                            fieldNameUpper,
                            transectName.replace('_', ' '))

                    # make a new subtask for this season and comparison grid
                    subtask = PlotTransectSubtask(self, season, transectName,
                                                  fieldName,
                                                  computeTransectsSubtask,
                                                  plotObs, refConfig)

                    subtask.set_plot_info(
                            outFileLabel=outFileLabel,
                            fieldNameInTitle=fieldNameInTytle,
                            mpasFieldName=fields[fieldName]['mpas'],
                            refFieldName=refFieldName,
                            refTitleLabel=refTitleLabel,
                            diffTitleLabel=diffTitleLabel,
                            unitsLabel=fields[fieldName]['units'],
                            imageCaption='{} {}'.format(fieldNameInTytle,
                                                        season),
                            galleryGroup='Geojson Transects',
                            groupSubtitle=None,
                            groupLink='geojson',
                            galleryName=fieldNameUpper,
                            configSectionName='geojson{}Transects'.format(
                                    fieldNameUpper))

                    self.add_subtask(subtask)
        # }}}

    # }}}


class GeojsonTransectsObservations(TransectsObservations):  # {{{
    """
    A class for loading and manipulating geojson transects

    Attributes
    ----------

    obsDatasets : OrderedDict
        A dictionary of observational datasets

    Authors
    -------
    Xylar Asay-Davis
    """

    def build_observational_dataset(self, fileName, transectName):  # {{{
        '''
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

        Authors
        -------
        Xylar Asay-Davis
        '''

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

        return dsObs  # }}}
    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
