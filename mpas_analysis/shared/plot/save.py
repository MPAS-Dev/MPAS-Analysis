# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import matplotlib.pyplot as plt


def savefig(filename):
    """
    Waves the current plot to a file, then closes it.

    Parameters
    ----------
    filename : str
        the file name to be written
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    plt.savefig(filename, dpi='figure', bbox_inches='tight', pad_inches=0.1)

    plt.close()


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
