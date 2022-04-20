# -*- coding: utf-8 -*-
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
#


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
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    ds = ds.rolling(Time=movingAveragePoints,
                    center=True).mean().dropna('Time')

    return ds

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
