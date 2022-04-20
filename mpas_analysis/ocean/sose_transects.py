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
import xarray as xr
import numpy
import os

from collections import OrderedDict

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.ocean.compute_transects_subtask import \
    TransectsObservations

from mpas_analysis.ocean.compute_transects_with_vel_mag import \
    ComputeTransectsWithVelMag

from mpas_analysis.ocean.plot_transect_subtask import PlotTransectSubtask

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories, build_obs_path
from mpas_analysis.shared.io import write_netcdf


class SoseTransects(AnalysisTask):
    """
    Plot model output at WOCE transects and compare it against WOCE
    observations
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

        tags = ['climatology', 'transect', 'sose', 'publicObs', 'antarctic']

        # call the constructor from the base class (AnalysisTask)
        super(SoseTransects, self).__init__(
            config=config, taskName='soseTransects',
            componentName='ocean',
            tags=tags)

        sectionName = self.taskName

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

        longitudes = sorted(config.getexpression(sectionName, 'longitudes',
                                                 use_numpyfunc=True))

        fields = \
            [{'prefix': 'temperature',
              'mpas': 'timeMonthly_avg_activeTracers_temperature',
              'units': r'$\degree$C',
              'titleName': 'Potential Temperature',
              'obsFilePrefix': 'pot_temp',
              'obsFieldName': 'theta'},
             {'prefix': 'salinity',
              'mpas': 'timeMonthly_avg_activeTracers_salinity',
              'units': r'PSU',
              'titleName': 'Salinity',
              'obsFilePrefix': 'salinity',
              'obsFieldName': 'salinity'},
             {'prefix': 'potentialDensity',
              'mpas': 'timeMonthly_avg_potentialDensity',
              'units': r'kg m$^{-3}$',
              'titleName': 'Potential Density',
              'obsFilePrefix': 'pot_den',
              'obsFieldName': 'potentialDensity'},
             {'prefix': 'zonalVelocity',
              'mpas': 'timeMonthly_avg_velocityZonal',
              'units': r'm s$^{-1}$',
              'titleName': 'Zonal Velocity',
              'obsFilePrefix': 'zonal_vel',
              'obsFieldName': 'zonalVel'},
             {'prefix': 'meridionalVelocity',
              'mpas': 'timeMonthly_avg_velocityMeridional',
              'units': r'm s$^{-1}$',
              'titleName': 'Meridional Velocity',
              'obsFilePrefix': 'merid_vel',
              'obsFieldName': 'meridVel'},
             {'prefix': 'velocityMagnitude',
              'mpas': 'velMag',
              'units': r'm s$^{-1}$',
              'titleName': 'Velocity Magnitude',
              'obsFilePrefix': None,
              'obsFieldName': 'velMag'}]

        fieldList = config.getexpression(sectionName, 'fieldList')
        fields = [field for field in fields if field['prefix'] in fieldList]

        variableList = [field['mpas'] for field in fields
                        if field['mpas'] != 'velMag']

        transectCollectionName = 'SOSE_transects'
        if horizontalResolution not in ['obs', 'mpas']:
            transectCollectionName = '{}_{}km'.format(transectCollectionName,
                                                      horizontalResolution)

        transectsObservations = SoseTransectsObservations(
            config, horizontalResolution,
            transectCollectionName, fields)

        computeTransectsSubtask = ComputeTransectsWithVelMag(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='SOSE_transects',
            transectCollectionName=transectCollectionName,
            variableList=variableList,
            seasons=seasons,
            obsDatasets=transectsObservations,
            verticalComparisonGridName=verticalComparisonGridName,
            verticalComparisonGrid=verticalComparisonGrid)

        plotObs = controlConfig is None
        if plotObs:

            refTitleLabel = 'State Estimate (SOSE)'

            diffTitleLabel = 'Model - State Estimate'

        else:
            controlRunName = controlConfig.get('runs', 'mainRunName')
            refTitleLabel = 'Control: {}'.format(controlRunName)

            diffTitleLabel = 'Main - Control'

        for field in fields:
            fieldPrefix = field['prefix']
            fieldPrefixUpper = fieldPrefix[0].upper() + fieldPrefix[1:]
            for lon in longitudes:
                transectName = 'lon_{}'.format(lon)

                for season in seasons:
                    outFileLabel = 'SOSE_{}_'.format(fieldPrefix)
                    if plotObs:
                        refFieldName = field['obsFieldName']
                    else:
                        refFieldName = field['mpas']

                    fieldNameInTytle = r'{} at {}$\degree$ Lon.'.format(
                        field['titleName'], lon)

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
                        imageCaption='{} {}'.format(fieldNameInTytle,
                                                    season),
                        galleryGroup='SOSE Transects',
                        groupSubtitle=None,
                        groupLink='sose_transects',
                        galleryName=field['titleName'],
                        configSectionName='sose{}Transects'.format(
                            fieldPrefixUpper),
                        verticalBounds=verticalBounds)

                    self.add_subtask(subtask)


class SoseTransectsObservations(TransectsObservations):
    """
    A class for loading and manipulating SOSE transect data

    Attributes
    ----------

    fields : list of dict
        dictionaries defining the fields with associated SOSE data
    """

    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, horizontalResolution,
                 transectCollectionName, fields):

        """
        Construct the object, setting the observations dictionary to None.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        horizontalResolution : str
            'obs' for the obs as they are or a size in km if subdivision is
            desired.

        transectCollectionName : str
            A name that describes the collection of transects (e.g. the name
            of the collection of observations) used to name the
            destination "mesh" for regridding

        fields : list of dict
            dictionaries defining the fields with associated SOSE data
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # this will be constructed later
        obsFileNames = None

        # call the constructor from the base class (TransectsObservations)
        super(SoseTransectsObservations, self).__init__(
            config, obsFileNames, horizontalResolution,
            transectCollectionName)

        self.fields = fields

    def get_observations(self):
        """
        Read in and set up the observations.

        Returns
        -------
        obsDatasets : OrderedDict
            The observational dataset
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, combine the SOSE observations into a single data set
        if self.obsFileNames is None:
            self.combine_observations()

        # then, call the method from the parent class
        return super(SoseTransectsObservations, self).get_observations()

    def combine_observations(self):
        """
        Combine SOSE oservations into a single file
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config

        longitudes = sorted(config.getexpression('soseTransects', 'longitudes',
                                                 use_numpyfunc=True))

        observationsDirectory = build_obs_path(
            config, 'ocean', 'soseSubdirectory')

        outObsDirectory = build_config_full_path(
            config=config, section='output',
            relativePathOption='climatologySubdirectory',
            relativePathSection='oceanObservations')

        make_directories(outObsDirectory)

        combinedFileName = '{}/{}.nc'.format(outObsDirectory,
                                             self.transectCollectionName)
        obsFileNames = OrderedDict()
        for lon in longitudes:
            transectName = 'lon_{}'.format(lon)
            obsFileNames[transectName] = combinedFileName

        self.obsFileNames = obsFileNames

        if os.path.exists(combinedFileName):
            return

        print('Preprocessing SOSE transect data...')

        minLat = config.getfloat('soseTransects', 'minLat')
        maxLat = config.getfloat('soseTransects', 'maxLat')

        dsObs = None
        for field in self.fields:
            prefix = field['obsFilePrefix']
            fieldName = field['obsFieldName']
            if prefix is None:
                continue
            print('  {}'.format(field['prefix']))

            fileName = '{}/SOSE_2005-2010_monthly_{}_SouthernOcean' \
                       '_0.167x0.167degree_20180710.nc'.format(
                observationsDirectory, prefix)

            dsLocal = xr.open_dataset(fileName)

            lat = dsLocal.lat.values
            mask = numpy.logical_and(lat >= minLat, lat <= maxLat)
            indices = numpy.argwhere(mask)
            dsLocal = dsLocal.isel(lat=slice(indices[0][0],
                                             indices[-1][0]))
            dsLocal.load()

            if fieldName == 'zonalVel':
                # need to average in longitude
                nLon = dsLocal.sizes['lon']
                lonIndicesP1 = numpy.mod(numpy.arange(nLon) + 1, nLon)
                dsLocal = 0.5 * (dsLocal + dsLocal.isel(lon=lonIndicesP1))

            if fieldName == 'meridVel':
                # need to average in latitude
                nLat = dsLocal.sizes['lat']
                latIndicesP1 = numpy.mod(numpy.arange(nLat) + 1, nLat)
                dsLocal = 0.5 * (dsLocal + dsLocal.isel(lat=latIndicesP1))

            dsLocal = dsLocal.sel(lon=longitudes, method=str('nearest'))

            if dsObs is None:
                dsObs = dsLocal
            else:
                dsLocal['lon'] = dsObs.lon
                dsLocal['lat'] = dsObs.lat
                dsObs[fieldName] = dsLocal[fieldName]
                dsLocal.close()

        if 'zonalVel' in dsObs and 'meridVel' in dsObs:
            # compute the velocity magnitude
            print('  velMag')
            description = 'Monthly velocity magnitude climatologies ' \
                          'from 2005-2010 average of the Southern Ocean ' \
                          'State Estimate (SOSE)'
            dsObs['velMag'] = numpy.sqrt(
                dsObs.zonalVel ** 2 + dsObs.meridVel ** 2)
            dsObs.velMag.attrs['units'] = 'm s$^{-1}$'
            dsObs.velMag.attrs['description'] = description

        # make a copy of the top set of data at z=0
        dsObs = xr.concat((dsObs.isel(z=0), dsObs), dim='z')
        z = dsObs.z.values
        z[0] = 0.
        dsObs['z'] = ('z', z)
        write_netcdf(dsObs, combinedFileName)

        print('  Done.')

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

        lon = float(transectName.rsplit('_', 1)[-1])

        dsObs = dsObs.sel(method=str('nearest'), lon=lon)
        lon = dsObs.lon.values

        # do some dropping and renaming so we end up with the right coordinates
        # and dimensions
        dsObs = dsObs.rename({'lat': 'nPoints', 'z': 'nz'})
        dsObs['lat'] = dsObs.nPoints
        dsObs['z'] = dsObs.nz
        dsObs['lon'] = ('nPoints', lon * numpy.ones(dsObs.sizes['nPoints']))
        dsObs = dsObs.drop_vars(['nPoints', 'nz'])

        return dsObs
