# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function, \
    unicode_literals


def compute_moving_avg(ds, movingAveragePoints=12):  # {{{

    '''
    Compute the rolling mean of a data set

    Parameters
    ----------
    ds :  ``xarray.Dataset``
        a dataset to be averaged

    movingAveragePoints : int, optional
        The number of points (months) over which to perform the rolling average
        of the data set

    Returns
    -------
    ds : ``xarray.Dataset``
        The anomaly of the rolling time mean from the start of the simulation

    Authors
    -------
    Xylar Asay-Davis
    '''

    ds = ds.rolling(Time=movingAveragePoints,
                    center=True).mean().dropna('Time')

    return ds

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
