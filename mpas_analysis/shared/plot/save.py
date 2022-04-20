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

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os

import matplotlib.pyplot as plt


def savefig(filename, config, tight=True, pad_inches=0.1):
    """
    Saves the current plot to a file, then closes it.

    Parameters
    ----------
    filename : str
        the file name to be written

    config : mpas_tools.config.MpasConfigParser
        Configuration options

    tight : bool, optional
        whether to tightly crop the figure

    pad_inches : float, optional
        The boarder around the image
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    if tight:
        bbox_inches = 'tight'
    else:
        bbox_inches = None

    filenames = [filename]

    if config.getboolean('plot', 'pdf'):
        pdf_filename = '{}.pdf'.format(os.path.splitext(filename)[0])
        filenames.append(pdf_filename)

    for path in filenames:
        plt.savefig(path, dpi='figure', bbox_inches=bbox_inches,
                    pad_inches=pad_inches)

    plt.close()


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
