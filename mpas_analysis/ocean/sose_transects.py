from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr
import numpy

from collections import OrderedDict

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.ocean.compute_transects_subtask import \
    ComputeTransectsSubtask, TransectsObservations

from mpas_analysis.ocean.plot_transect_subtask import PlotTransectSubtask

from mpas_analysis.shared.io.utility import build_config_full_path


class SoseTransects(AnalysisTask):  # {{{
    """
    Plot model output at WOCE transects and compare it against WOCE
    observations

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

        tags = ['climatology', 'transect', 'sose']

        # call the constructor from the base class (AnalysisTask)
        super(SoseTransects, self).__init__(
                config=config, taskName='soseTransects',
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

        observationsDirectory = build_config_full_path(
            config, 'oceanObservations', 'soseSubdirectory')

        longitudes = config.getExpression(sectionName, 'longitudes',
                                          usenumpyfunc=True)

        fields = \
            {'temperature':
                {'mpas': 'timeMonthly_avg_activeTracers_temperature',
                 'obs': 'theta',
                 'obsFile': 'SOSE_2005-2010_monthly_pot_temp_SouthernOcean'
                            '_0.167x0.167degree.nc',
                 'units': r'$\degree$C'},
             'salinity':
                {'mpas': 'timeMonthly_avg_activeTracers_salinity',
                 'obs': 'salinity',
                 'obsFile': 'SOSE_2005-2010_monthly_salinity_SouthernOcean'
                            '_0.167x0.167degree.nc',
                 'units': r'PSU'}}

        transectCollectionName = 'SOSE_transects'
        if horizontalResolution != 'obs':
            transectCollectionName = '{}_{}km'.format(transectCollectionName,
                                                      horizontalResolution)

        fieldNameDict = {'temperature':'temperatureTransect',
                         'salinity':'salinityTransect'}

        obsFileNames = OrderedDict()
        for fieldName in fields:
            for lon in longitudes:
                transectName = '{}_SOSE_lon_{}'.format(fieldNameDict[fieldName],
                                                       lon)

                fileName = '{}/{}'.format(observationsDirectory,
                                          fields[fieldName]['obsFile'])
                obsFileNames[transectName] = fileName

        transectsObservations = SoseTransectsObservations(
                config, obsFileNames, horizontalResolution,
                transectCollectionName)

        computeTransectsSubtask = ComputeTransectsSubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='SOSE',
            transectCollectionName=transectCollectionName,
            variableList=[field['mpas'] for field in fields.values()],
            seasons=seasons,
            obsDatasets=transectsObservations,
            verticalComparisonGridName=verticalComparisonGridName,
            verticalComparisonGrid=verticalComparisonGrid)

        plotObs = refConfig is None
        if plotObs:

            refTitleLabel = 'State Estimate (SOSE)'

            diffTitleLabel = 'Model - State Estimate'

        else:
            refRunName = refConfig.get('runs', 'mainRunName')
            refTitleLabel = 'Ref: {}'.format(refRunName)

            diffTitleLabel = 'Main - Reference'

        for fieldName in fields:
            for lon in longitudes:
                transectName = '{}_SOSE_lon_{}'.format(fieldNameDict[fieldName],
                                                       lon)
                for season in seasons:
                    outFileLabel = ''
                    if plotObs:
                        refFieldName = fields[fieldName]['obs']
                    else:
                        refFieldName = fields[fieldName]['mpas']

                    fieldNameUpper = fieldName[0].upper() + fieldName[1:]
                    fieldNameInTytle = r'{} at {}$\degree$ Longitude'.format(
                            fieldNameUpper, lon)

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
                            galleryGroup='SOSE Transects',
                            groupSubtitle=None,
                            groupLink='sose_transects',
                            galleryName=fieldNameUpper,
                            configSectionName='sose{}Transects'.format(
                                    fieldNameUpper))

                    self.add_subtask(subtask)
        # }}}

    # }}}


class SoseTransectsObservations(TransectsObservations):  # {{{
    """
    A class for loading and manipulating SOSE transect observations

    Attributes
    ----------

    obsDatasets : OrderedDict
        A dictionary of observational datasets

    Authors
    -------
    Xylar Asay-Davis
    """

    def __init__(self, config, obsFileName, horizontalResolution,
                 transectCollectionName):  # {{{
        '''
        Construct the object, setting the observations dictionary to None.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        obsFileNames : OrderedDict
            The names of transects and the file names of the corresponding
            observations for a transect

        horizontalResolution : str
            'obs' for the obs as they are or a size in km if subdivision is
            desired.

        transectCollectionName : str
            A name that describes the collection of transects (e.g. the name
            of the collection of observations) used to name the
            destination "mesh" for regridding

        Authors
        -------
        Xylar Asay-Davis
        '''
        # call the constructor from the base class (TransectsObservations)
        super(SoseTransectsObservations, self).__init__(
            config, obsFileName, horizontalResolution,
            transectCollectionName)

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

        config = self.config
        minLat = config.getfloat('soseTransects', 'minLat')
        maxLat = config.getfloat('soseTransects', 'maxLat')

        dsObs = xr.open_dataset(fileName)

        lon = float(transectName.rsplit('_', 1)[-1])

        dsObs = dsObs.sel(method=str('nearest'), lon=lon)
        lon = dsObs.lon.values

        lat = dsObs.lat.values
        mask = numpy.logical_and(lat >= minLat, lat <= maxLat)
        indices = numpy.argwhere(mask)
        dsObs = dsObs.isel(lat=slice(indices[0][0], indices[-1][0]))

        # do some dropping and renaming so we end up wiht the right coordinates
        # and dimensions
        dsObs = dsObs.rename({'lat': 'nPoints', 'z': 'nz'})
        dsObs['lat'] = dsObs.nPoints
        dsObs['z'] = dsObs.nz
        dsObs['lon'] = ('nPoints', lon*numpy.ones(dsObs.sizes['nPoints']))
        dsObs = dsObs.drop(['nPoints', 'nz'])

        return dsObs  # }}}
    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
