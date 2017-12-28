
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os
import xarray

from .mpas_time_series_task import MpasTimeSeriesTask

from ..io.utility import build_config_full_path


class MpasReferenceTimeSeriesTask(MpasTimeSeriesTask):  # {{{
    '''
    Ensure that the time series file for a reference task exists

    Authors
    -------
    Xylar Asay-Davis
    '''

    def __init__(self, config, componentName, section='timeSeries'):  # {{{
        '''
        Construct the analysis task for extracting time series.

        Parameters
        ----------
        config : ``MpasAnalysisConfigParser``
            Contains configuration options

        componentName : {'ocean', 'seaIce'}
            The name of the component (same as the folder where the task
            resides)

        section : str, optional
            The section of the config file from which to read the start and
            end times for the time series, also added as a tag

        Authors
        -------
        Xylar Asay-Davis
        '''

        suffix = section[0].upper() + section[1:] + \
            componentName[0].upper() + componentName[1:]

        taskName = 'mpasRef{}'.format(suffix)
        self.outputPrefix = 'mpas{}'.format(suffix)

        # call the constructor from the base class (MpasTimeSeriesTask)
        super(MpasReferenceTimeSeriesTask, self).__init__(
            config=config,
            componentName=componentName,
            taskName=taskName,
            section=section)
        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        '''

        # first, call setup_and_check from the base class (AnalysisTask),
        # skipping MpasTimeSeriesTask because we don't to check for the
        # input files, only the output file
        super(MpasTimeSeriesTask, self).setup_and_check()

        config = self.config
        baseDirectory = build_config_full_path(
            config, 'output', 'timeSeriesSubdirectory')

        self.outputFile = '{}/{}.nc'.format(baseDirectory,
                                            self.outputPrefix)

        if not os.path.exists(self.outputFile):
            raise OSError('Missing reference time series file {}'.format(
                    self.outputFile))

        ds = xarray.open_dataset(self.outputFile)

        for variableName in self.variableList:
            if variableName not in ds.variables:
                raise IOError('Missing variable {} in reference climatology '
                              'file {}'.format(variableName,
                                               self.outputFile))

        # }}}

    def run_task(self):  # {{{
        '''
        Nothing to do, since the output file already exists

        Authors
        -------
        Xylar Asay-Davis
        '''
        pass  # }}}

    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
