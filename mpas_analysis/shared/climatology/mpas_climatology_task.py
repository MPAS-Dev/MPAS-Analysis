import xarray as xr
import os
import warnings
import sys
import subprocess
from distutils.spawn import find_executable

from ..analysis_task import AnalysisTask

from ..constants import constants

from ..io.utility import build_config_full_path, make_directories
from ..io import write_netcdf

from .climatology import get_remapper
from .comparison_descriptors import get_comparison_descriptor

from ..grid import MpasMeshDescriptor

from ..mpas_xarray import mpas_xarray


class MpasClimatology(AnalysisTask):  # {{{
    '''
    An analysis tasks for computing climatologies from output from the
    ``timeSeriesStatsMonthly`` analysis member.

    Attributes
    ----------

    taskSuffix :  str
        The suffix to append to the task name, typically a short name for
        the field(s) being analyzed.  For clarity, the taskSuffix should
        start with a capital letter.

    variableList : list of str
        A list of variable names in ``timeSeriesStatsMonthly`` to be
        included in the climatologies

    iselValues : dict
        A dictionary of dimensions and indices (or ``None``) used to extract
        a slice of the MPAS field.

    seasons : list of str
        A list of seasons (keys in ``shared.constants.monthDictionary``)
        over which the climatology should be computed or ['none'] if only
        monthly climatologies are needed.

    inputFiles : list of str
        A list of input files used to compute the climatologies.

    comparisonGridNames : list of {``None``, 'latlon', 'antarctic'}
        The name(s) of the comparison grid to use for remapping.

    restartFileName : str
        If ``comparisonGridName`` is not ``None``, the name of a restart
        file from which the MPAS mesh can be read.

    ncclimoModel : {'mpaso', 'mpascice'}
        The name of the component expected by ``ncclimo``

    startDate, endDate : str
        The start and end dates of the climatology as strings

    startYear, endYear : int
        The start and end years of the climatology

    fillValue : float
        The fill value used in MPAS output (but currently not written to the
        ``_FillValue`` attribute)

    Authors
    -------
    Xylar Asay-Davis
    '''

    def __init__(self, config, variableList, taskSuffix,
                 componentName, comparisonGridNames=None,
                 seasons=['none'], tags=None, iselValues=None):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        variableList : list of str
            A list of variable names in ``timeSeriesStatsMonthly`` to be
            included in the climatologies

        taskSuffix :  str
            The suffix to append to the task name, typically a short name for
            the field(s) being analyzed.  For clarity, the taskSuffix should
            start with a capital letter.

        componentName :  {'ocean', 'seaIce'}
            The name of the component (same as the folder where the task
            resides)

        comparisonGridNames : list of {'latlon', 'antarctic'}, optional
            The name(s) of the comparison grid to use for remapping.

        seasons : list of str, optional
            A list of seasons (keys in ``shared.constants.monthDictionary``)
            to be computed or ['none'] (not ``None``) if only monthly
            climatologies are needed.

        tags :  list of str, optional
            Tags used to describe the task (e.g. 'timeSeries', 'climatology',
            horizontalMap', 'index', 'transect').  These are used to determine
            which tasks are generated (e.g. 'all_transect' or 'no_climatology'
            in the 'generate' flags)

        iselValues : dict, optional
            A dictionary of dimensions and indices (or ``None``) used to
            extract a slice of the MPAS field(s).

        Authors
        -------
        Xylar Asay-Davis
        '''
        self.variableList = variableList
        self.seasons = seasons
        self.taskSuffix = taskSuffix
        self.comparisonGridNames = comparisonGridNames
        self.iselValues = iselValues

        # this is a stopgap until MPAS implements the _FillValue attribute
        # correctly
        self.fillValue = -9.99999979021476795361e+33

        if 'climatology' not in tags:
            tags.append('climatology')

        if componentName == 'ocean':
            self.ncclimoModel = 'mpaso'
        elif componentName == 'seaIce':
            self.ncclimoModel = 'mpascice'
        else:
            raise ValueError('component {} is not supported by ncclimo.\n'
                             'Check with Charlie Zender and Xylar Asay-Davis\n'
                             'about getting it added'.format(componentName))

        # call the constructor from the base class (AnalysisTask)
        super(MpasClimatology, self).__init__(
            config=config,
            taskName='mpasClimatology{}'.format(taskSuffix),
            componentName=componentName,
            tags=tags)

        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        '''
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(MpasClimatology, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        self.restartFileName = None
        if self.comparisonGridNames is not None:
            try:
                self.restartFileName = self.runStreams.readpath('restart')[0]
            except ValueError:
                raise IOError('No MPAS restart file found: need at least one '
                              'restart file to perform remapping of '
                              'climatologies.')

        # get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        startDate = self.config.get('climatology', 'startDate')
        endDate = self.config.get('climatology', 'endDate')
        streamName = \
            self.historyStreams.find_stream(self.streamMap['timeSeriesStats'])
        self.inputFiles = self.historyStreams.readpath(
                streamName, startDate=startDate, endDate=endDate,
                calendar=self.calendar)

        if len(self.inputFiles) == 0:
            raise IOError('No files were found in stream {} between {} and '
                          '{}.'.format(streamName, startDate, endDate))

        self._update_climatology_bounds_from_file_names()

        # }}}

    def run(self):  # {{{
        '''
        Compute the requested climatologies

        Authors
        -------
        Xylar Asay-Davis
        '''

        print '\nComputing {} climatologies from files:\n' \
              '    {} through\n    {}'.format(
                  self.taskSuffix,
                  os.path.basename(self.inputFiles[0]),
                  os.path.basename(self.inputFiles[-1]))

        config = self.config

        mpasMeshName = config.get('input', 'mpasMeshName')

        if self.seasons[0] is 'none':
            seasonsToCheck = ['{:02d}'.format(month) for month in range(1, 13)]
        else:
            seasonsToCheck = self.seasons

        allExist = True
        for season in seasonsToCheck:

            climatologyFileName, climatologyDirectory = \
                self.get_ncclimo_file_name(season, 'unmasked',
                                           returnDir=True)

            if not os.path.exists(climatologyFileName):
                allExist = False
                break

        if not allExist:
            self._compute_climatologies_with_ncclimo(
                    inDirectory=self.historyDirectory,
                    outDirectory=climatologyDirectory)

        if self.comparisonGridNames is not None:

            parallel = self.config.getint('execute', 'parallelTaskCount') > 1
            if parallel:
                # avoid writing the same mapping file from multiple processes
                mappingFilePrefix = 'map_{}'.format(self.taskName)
            else:
                mappingFilePrefix = 'map'

            mpasDescriptor = MpasMeshDescriptor(
                    self.restartFileName,
                    meshName=mpasMeshName)

            dsMask = xr.open_dataset(self.inputFiles[0])
            dsMask = mpas_xarray.subset_variables(dsMask, self.variableList)
            iselValues = {'Time': 0}
            if self.iselValues is not None:
                for dim in self.iselValues:
                    # we've already hyperslabbed this dimension in ncclimo
                    iselValues[dim] = 0
            # select only Time=0 and possibly only the desired vertical
            # slice
            dsMask = dsMask.isel(**iselValues)

            firstGrid = True
            for comparisonGridName in self.comparisonGridNames:
                comparisonDescriptor = \
                    get_comparison_descriptor(config, comparisonGridName)

                mpasRemapper = get_remapper(
                    config=config, sourceDescriptor=mpasDescriptor,
                    comparisonDescriptor=comparisonDescriptor,
                    mappingFilePrefix=mappingFilePrefix,
                    method=config.get('climatology',
                                      'mpasInterpolationMethod'))

                for season in self.seasons:
                    if firstGrid:
                        self._mask_climatologies(season, dsMask,
                                                 comparisonDescriptor)

                    maskedClimatologyFileName = self.get_ncclimo_file_name(
                            season, 'masked', comparisonDescriptor)

                    remappedFileName = self.get_ncclimo_file_name(
                            season, 'remapped', comparisonDescriptor)

                    if not os.path.exists(remappedFileName):
                        self._remap(inFileName=maskedClimatologyFileName,
                                    outFileName=remappedFileName,
                                    remapper=mpasRemapper,
                                    comparisonGridName=comparisonGridName)

                firstGrid = False
        # }}}

    def get_ncclimo_file_name(self, season, stage, comparisonDescriptor=None,
                              returnDir=False):  # {{{
        """
        Given config options, the name of a field and a string identifying the
        months in a seasonal climatology, returns the full path for MPAS
        climatology files before and after remapping.

        Parameters
        ----------
        season : str
            One of the seasons in ``constants.monthDictionary``

        mpasMeshName : str
            The name of the MPAS mesh

        stage : {'unmasked', 'masked', 'remapped'}
            The stage of the masking and remapping process

        comparisonDescriptor : MeshDescriptor, optional
            The comparison mesh descriptor, used to get the mesh name

        returnDir : bool, optional
            Return the directory as well

        Returns
        -------
        fileName : str
            The path to the climatology file for the specified season.

        Authors
        -------
        Xylar Asay-Davis
        """

        climatologyBaseDirectory = build_config_full_path(
            self.config, 'output', 'mpasClimatologySubdirectory')

        mpasMeshName = self.config.get('input', 'mpasMeshName')

        climatologyBaseDirectory = '{}/{}'.format(climatologyBaseDirectory,
                                                  stage)

        if stage in ['unmasked', 'masked']:
            directory = '{}/{}_{}'.format(
                    climatologyBaseDirectory, self.taskSuffix, mpasMeshName)
        elif stage == 'remapped':
            directory = '{}/{}_{}_to_{}'.format(
                    climatologyBaseDirectory, self.taskSuffix, mpasMeshName,
                    comparisonDescriptor.meshName)
        else:
            raise ValueError('Unsupported stage {}'.format(stage))

        make_directories(directory)

        monthValues = sorted(constants.monthDictionary[season])
        startMonth = monthValues[0]
        endMonth = monthValues[-1]

        suffix = '{:04d}{:02d}_{:04d}{:02d}_climo'.format(
                self.startYear, startMonth, self.endYear, endMonth)

        if season in constants.abrevMonthNames:
            season = '{:02d}'.format(monthValues[0])
        fileName = '{}/{}_{}_{}.nc'.format(directory, self.ncclimoModel,
                                           season, suffix)

        if returnDir:
            return fileName, directory
        else:
            return fileName  # }}}

    def _update_climatology_bounds_from_file_names(self):  # {{{
        """
        Update the start and end years and dates for climatologies based on the
        years actually available in the list of files.

        Authors
        -------
        Xylar Asay-Davis
        """

        config = self.config

        requestedStartYear = config.getint('climatology', 'startYear')
        requestedEndYear = config.getint('climatology', 'endYear')

        dates = sorted([fileName[-13:-6] for fileName in self.inputFiles])
        years = [int(date[0:4]) for date in dates]
        months = [int(date[5:7]) for date in dates]

        # search for the start of the first full year
        firstIndex = 0
        while(firstIndex < len(years) and months[firstIndex] != 1):
            firstIndex += 1
        startYear = years[firstIndex]

        # search for the end of the last full year
        lastIndex = len(years)-1
        while(lastIndex >= 0 and months[lastIndex] != 12):
            lastIndex -= 1
        endYear = years[lastIndex]

        if startYear != requestedStartYear or endYear != requestedEndYear:
            message = "climatology start and/or end year different from " \
                      "requested\n" \
                      "requestd: {:04d}-{:04d}\n" \
                      "actual:   {:04d}-{:04d}\n".format(requestedStartYear,
                                                         requestedEndYear,
                                                         startYear,
                                                         endYear)
            warnings.warn(message)
            config.set('climatology', 'startYear', str(startYear))
            config.set('climatology', 'endYear', str(endYear))

            startDate = '{:04d}-01-01_00:00:00'.format(startYear)
            config.set('climatology', 'startDate', startDate)
            endDate = '{:04d}-12-31_23:59:59'.format(endYear)
            config.set('climatology', 'endDate', endDate)

        else:
            startDate = config.get('climatology', 'startDate')
            endDate = config.get('climatology', 'endDate')

        self.startDate = startDate
        self.endDate = endDate
        self.startYear = startYear
        self.endYear = endYear

        # }}}

    def _compute_climatologies_with_ncclimo(self, inDirectory, outDirectory,
                                            remapper=None,
                                            remappedDirectory=None):  # {{{
        '''
        Uses ncclimo to compute monthly, seasonal and/or annual climatologies.

        Parameters
        ----------
        inDirectory : str
            The run directory containing timeSeriesStatsMonthly output

        outDirectory : str
            The output directory where climatologies will be written

        remapper : ``shared.intrpolation.Remapper`` object, optional
            If present, a remapper that defines the source and desitnation
            grids for remapping the climatologies.

        remappedDirectory : str, optional
            If present, the path where remapped climatologies should be
            written. By default, remapped files are stored in the same
            directory as the climatologies on the source grid.  Has no effect
            if ``remapper`` is ``None``.

        Raises
        ------
        OSError
            If ``ncclimo`` is not in the system path.

        Author
        ------
        Xylar Asay-Davis
        '''

        if find_executable('ncclimo') is None:
            raise OSError('ncclimo not found. Make sure the latest nco '
                          'package is installed: \n'
                          'conda install nco\n'
                          'Note: this presumes use of the conda-forge '
                          'channel.')

        parallelMode = self.config.get('execute', 'ncclimoParallelMode')

        args = ['ncclimo',
                '--clm_md=mth',
                '-a', 'sdd',
                '-m', self.ncclimoModel,
                '-p', parallelMode,
                '-v', ','.join(self.variableList),
                '--seasons={}'.format(','.join(self.seasons)),
                '-s', '{:04d}'.format(self.startYear),
                '-e', '{:04d}'.format(self.endYear),
                '-i', inDirectory,
                '-o', outDirectory]

        if remapper is not None:
            args.extend(['-r', remapper.mappingFileName])
            if remappedDirectory is not None:
                args.extend(['-O', remappedDirectory])

        if self.iselValues is not None:
            ncksOptions = ['-O', '--no_tmp_fl']

            for dim in self.iselValues:
                val = self.iselValues[dim]
                ncksOptions.extend(['-d', '{},{},{}'.format(dim, val, val)])

            args.extend(['-n', ' '.join(ncksOptions)])

        # make sure any output is flushed before we add output from the
        # subprocess
        sys.stdout.flush()
        sys.stderr.flush()

        subprocess.check_call(args)  # }}}

    def _mask_climatologies(self, season, dsMask, comparisonDescriptor):  # {{{
        '''
        For each season, creates a masked version of the climatology

        Parameters
        ----------
        season : str
            The name of the season to be masked

        dsMask : ``xarray.Dataset`` object
            A data set (from the first input file) that can be used to
            determine the mask in MPAS output files.

        comparisonDescriptor : MeshDescriptor, optional
            The comparison mesh descriptor, used to get the mesh name

        Author
        ------
        Xylar Asay-Davis
        '''

        climatologyFileName = self.get_ncclimo_file_name(
                season, 'unmasked', comparisonDescriptor)

        maskedClimatologyFileName = self.get_ncclimo_file_name(
                season, 'masked', comparisonDescriptor)

        if not os.path.exists(maskedClimatologyFileName):
            # slice and mask the data set
            climatology = xr.open_dataset(climatologyFileName)
            iselValues = {'Time': 0}
            if self.iselValues is not None:
                iselValues.update(self.iselValues)
            # select only Time=0 and possibly only the desired vertical
            # slice
            climatology = climatology.isel(**iselValues)

            # mask the data set
            for variableName in self.variableList:
                climatology[variableName] = \
                    climatology[variableName].where(
                        dsMask[variableName] != self.fillValue)

            write_netcdf(climatology, maskedClimatologyFileName)
        # }}}

    def _remap(self, inFileName, outFileName, remapper, comparisonGridName):
        # {{{
        """
        Performs remapping either using ``ncremap`` or the native python code,
        depending on the requested setting and the comparison grid

        Parameters
        ----------
        inFileName : str
            The name of the input file to be remapped.

        outFileName : str
            The name of the output file to which the remapped data set should
            be written.

        remapper : ``Remapper`` object
            A remapper that can be used to remap files or data sets to a
            comparison grid.

        comparisonGridNames : {'latlon', 'antarctic'}
            The name of the comparison grid to use for remapping.

        Authors
        -------
        Xylar Asay-Davis
        """
        if remapper.mappingFileName is None:
            # no remapping is needed
            return

        useNcremap = self.config.getboolean('climatology', 'useNcremap')

        if comparisonGridName == 'antarctic':
            # ncremap doesn't support polar stereographic grids
            useNcremap = False

        renormalizationThreshold = self.config.getfloat(
            'climatology', 'renormalizationThreshold')

        if useNcremap:
            remapper.remap_file(inFileName=inFileName,
                                outFileName=outFileName,
                                overwrite=True,
                                renormalize=renormalizationThreshold)
        else:

            climatologyDataSet = xr.open_dataset(inFileName)

            remappedClimatology = remapper.remap(climatologyDataSet,
                                                 renormalizationThreshold)
            write_netcdf(remappedClimatology, outFileName)
        # }}}

    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
