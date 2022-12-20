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

from io import StringIO

from mpas_tools.config import MpasConfigParser

from mpas_analysis.shared.climatology import MpasClimatologyTask
from mpas_analysis.shared.timekeeping.utility import get_simulation_start_time


class RefYearMpasClimatologyTask(MpasClimatologyTask):
    """
    An analysis tasks for computing a reference-year climatology for computing
    climatology anomalies from the``timeSeriesStatsMonthly`` analysis member.

    Attributes
    ----------

    anomalyRefYear : int
        The reference year over which to compute the climatology
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, componentName, taskName=None):
        """
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Contains configuration options

        componentName : {'ocean', 'seaIce'}
            The name of the component (same as the folder where the task
            resides)

        taskName : str, optional
            the name of the task, defaults to mpasClimatology<ComponentName>
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if taskName is None:
            suffix = componentName[0].upper() + componentName[1:]
            taskName = 'refYearMpasClimatology{}'.format(suffix)

        # make a deep copy of the config so we can change the start and end
        # years and dates without causing trouble for other tasks

        new_config = config.copy()

        # call the constructor from the base class (AnalysisTask)
        super(RefYearMpasClimatologyTask, self).__init__(
            config=new_config,
            componentName=componentName,
            taskName=taskName)

        self.tags.append('anomaly')

    def get_start_and_end(self):
        """
        Get the start and end years and dates for the climatology.  This
        function is provided to allow a custom method for setting the start
        and end years of the climatology.

        Returns
        -------
        startYear, endYear : int
           The start and end years of the climatology

        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config
        if config.has_option('climatology', 'anomalyRefYear'):
            anomalyRefYear = config.getint('climatology', 'anomalyRefYear')
        else:
            # we want the start year of the simulation, which we can't get yet
            anomalyRefDate = get_simulation_start_time(self.runStreams)
            anomalyRefYear = int(anomalyRefDate[0:4])

        startYear = anomalyRefYear
        endYear = anomalyRefYear

        startDate = '{:04d}-01-01_00:00:00'.format(startYear)
        endDate = '{:04d}-12-31_23:59:59'.format(endYear)

        # using "user=True" so these changes take priority over user config
        # options
        config.set('climatology', 'startYear', str(startYear), user=True)
        config.set('climatology', 'startDate', startDate, user=True)
        config.set('climatology', 'endYear', str(endYear), user=True)
        config.set('climatology', 'endDate', endDate, user=True)

        return startYear, endYear


