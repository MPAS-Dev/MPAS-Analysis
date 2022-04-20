# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2020 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2020 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2020 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
# Authors
# -------
# Xylar Asay-Davis, Milena Veneziani, Luke Van Roekel, Greg Streletz

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, FixedLocator
import numpy as np
from functools import partial

from mpas_analysis.shared.timekeeping.utility import days_to_datetime, \
    date_to_days


def plot_xtick_format(calendar, minDays, maxDays, maxXTicks, yearStride=None):
    '''
    Formats tick labels and positions along the x-axis for time series
    / index plots

    Parameters
    ----------
    calendar : str
        the calendar to use for formatting the time axis

    minDays : float
        start time for labels

    maxDays : float
        end time for labels

    maxXTicks : int
        the maximum number of tick marks to display, used to sub-sample ticks
        if there are too many

    yearStride : int, optional
        the number of years to skip over between ticks
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    ax = plt.gca()

    start = days_to_datetime(np.amin(minDays), calendar=calendar)
    end = days_to_datetime(np.amax(maxDays), calendar=calendar)

    if yearStride is not None or end.year - start.year > maxXTicks / 2:
        if yearStride is None:
            yearStride = 1
        else:
            maxXTicks = None
        major = [date_to_days(year=year, calendar=calendar)
                 for year in np.arange(start.year, end.year + 1, yearStride)]
        formatterFun = partial(_date_tick, calendar=calendar,
                               includeMonth=False)
    else:
        # add ticks for months
        major = []
        for year in range(start.year, end.year + 1):
            for month in range(1, 13):
                major.append(date_to_days(year=year, month=month,
                                          calendar=calendar))
        formatterFun = partial(_date_tick, calendar=calendar,
                               includeMonth=True)

    ax.xaxis.set_major_locator(FixedLocator(major, maxXTicks))
    ax.xaxis.set_major_formatter(FuncFormatter(formatterFun))

    plt.setp(ax.get_xticklabels(), rotation=30)

    plt.autoscale(enable=True, axis='x', tight=True)


def _date_tick(days, pos, calendar='gregorian', includeMonth=True):
    days = np.maximum(days, 0.)
    date = days_to_datetime(days, calendar)
    if includeMonth:
        return '{:04d}-{:02d}'.format(date.year, date.month)
    else:
        return '{:04d}'.format(date.year)


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
