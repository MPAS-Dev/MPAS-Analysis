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
"""
Regression tests for helpers in ``mpas_analysis.__main__``.
"""

import os
from unittest.mock import Mock, patch

from mpas_analysis.test import TestCase


# Importing mpas_analysis.__main__ triggers matplotlib imports in some test
# environments, so use a writable cache directory.
os.environ.setdefault('MPLCONFIGDIR', '/tmp/matplotlib')

import mpas_analysis.__main__ as main


class TestMain(TestCase):
    def test_get_editable_install_dir_without_direct_url(self):
        distribution = Mock()
        distribution.read_text.return_value = None

        with patch.object(main.Distribution, 'from_name',
                          return_value=distribution):
            self.assertEqual(main.get_editable_install_dir('mpas_analysis'),
                             None)
