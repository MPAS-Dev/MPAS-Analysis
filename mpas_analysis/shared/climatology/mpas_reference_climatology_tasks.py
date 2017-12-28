
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray
import os

from .mpas_climatology_task import MpasClimatologyTask

from .remap_mpas_climatology_subtask import RemapMpasClimatologySubtask


class MpasReferenceClimatologyTask(MpasClimatologyTask):  # {{{
    '''
    An analysis tasks for making sure climatologies from a reference run have
    been computed.

    Authors
    -------
    Xylar Asay-Davis
    '''

    def __init__(self, config, componentName):  # {{{
        '''
        Construct the analysis task

        Parameters
        ----------
        config : ``MpasAnalysisConfigParser``
            Contains configuration options

        componentName : {'ocean', 'seaIce'}
            The name of the component (same as the folder where the task
            resides)

        Authors
        -------
        Xylar Asay-Davis
        '''
        suffix = componentName[0].upper() + componentName[1:]
        taskName = 'refMpasClimatology{}'.format(suffix)

        # call the constructor from the base class (MpasClimatologyTask)
        super(MpasReferenceClimatologyTask, self).__init__(
            config=config,
            componentName=componentName,
            taskName=taskName)

        # }}}

    def setup_and_check(self):  # {{{
        '''
        Check if the necessary files exist.

        Raises
        ------
        IOError :
            If a climatology file is missing

        Authors
        -------
        Xylar Asay-Davis
        '''

        # first, call setup_and_check from the parent class
        # (MpasClimatologyTask)
        super(MpasReferenceClimatologyTask, self).setup_and_check()

        if len(self.seasons) == 0:
            return

        if self.seasons[0] is 'none':
            seasonsToCheck = ['{:02d}'.format(month) for month in range(1, 13)]
        else:
            seasonsToCheck = self.seasons

        for season in seasonsToCheck:

            climatologyFileName, climatologyDirectory = \
                self.get_file_name(season, returnDir=True)

            if not os.path.exists(climatologyFileName):
                raise IOError('Missing reference climatology file {}'.format(
                        climatologyFileName))

        # make sure all the necessary variables are also present
        climatologyFileName = self.get_file_name(seasonsToCheck[0],
                                                 returnDir=False)
        ds = xarray.open_dataset(climatologyFileName)

        for variableName in self.variableList:
            if variableName not in ds.variables:
                raise IOError('Missing variable {} in reference climatology '
                              'file {}'.format(variableName,
                                               climatologyFileName))
        # }}}

    def run_task(self):  # {{{
        '''
        Do nothing, since the climatologies already exist

        Authors
        -------
        Xylar Asay-Davis
        '''
        pass  # }}}

    # }}}


class RemapMpasReferenceClimatologySubtask(RemapMpasClimatologySubtask):  # {{{
    '''
    An analysis tasks for making sure climatologies from a reference run have
    been masked and remapped.

    Attributes
    ----------

    Authors
    -------
    Xylar Asay-Davis
    '''

    def __init__(self, mpasClimatologyTask, parentTask, climatologyName,
                 variableList, seasons, comparisonGridNames=['latlon'],
                 iselValues=None):
        # {{{
        '''
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped

        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        climatologyName : str
            A name that describes the climatology (e.g. a short version of
            the important field(s) in the climatology) used to name the
            subdirectories for each stage of the climatology

        variableList : list of str
            A list of variable names in ``timeSeriesStatsMonthly`` to be
            included in the climatologies

        seasons : list of str, optional
            A list of seasons (keys in ``shared.constants.monthDictionary``)
            to be computed or ['none'] (not ``None``) if only monthly
            climatologies are needed.

        comparisonGridNames : list of {'latlon', 'antarctic'}, optional
            The name(s) of the comparison grid to use for remapping.

        iselValues : dict, optional
            A dictionary of dimensions and indices (or ``None``) used to
            extract a slice of the MPAS field(s).

        Authors
        -------
        Xylar Asay-Davis
        '''

        # call the constructor from the base class
        # (RemapMpasClimatologySubtask)
        super(RemapMpasReferenceClimatologySubtask, self).__init__(
            mpasClimatologyTask, parentTask, climatologyName,
            variableList, seasons, comparisonGridNames, iselValues,
            subtaskName='remapMpasRefClimatology')

        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        IOError :
            If a remapped climatology file is missing

        Authors
        -------
        Xylar Asay-Davis
        '''

        # first, call setup_and_check from the parent class (AnalysisTask),
        # skipping setup_and_check from RemapMpasClimatologySubtask because
        # we don't want to create remappers.
        super(RemapMpasClimatologySubtask, self).setup_and_check()

        self._setup_file_names()

        for comparisonGridName in self.comparisonGridNames:

            for season in self.seasons:

                remappedFileName = self.get_file_name(
                        season, 'remapped', comparisonGridName)

                if not os.path.exists(remappedFileName):
                    raise OSError('Remapped climatology file {} from the '
                                  'reference run does not exist'.format(
                                          remappedFileName))

        # }}}

    def run_task(self):  # {{{
        '''
        Do nothing because the files already exist

        Authors
        -------
        Xylar Asay-Davis
        '''

        pass  # }}}
    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
