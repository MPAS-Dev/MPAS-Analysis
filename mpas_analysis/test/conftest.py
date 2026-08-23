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

import pytest


@pytest.fixture(autouse=True)
def run_in_tmp_dir(tmp_path, monkeypatch):
    """
    Run every test in its own temporary working directory.

    Some tests -- and the ``ESMF_RegridWeightGen`` subprocess that pyremap
    launches, which writes ``PET*.RegridWeightGen.Log`` -- create files
    relative to the current working directory.  Without this fixture, running
    the test suite leaves those files wherever pytest happened to be started,
    which is usually the root of the repository.

    All tests locate their inputs through absolute paths (``datadir`` and
    ``test_dir`` are both temporary directories), so changing the working
    directory is safe.
    """
    monkeypatch.chdir(tmp_path)
