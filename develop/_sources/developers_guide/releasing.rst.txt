.. _dev_releasing:

***********************
Releasing a New Version
***********************

This document describes the steps for maintainers to tag and release a new
version of ``MPAS-Analysis``, and to update the conda-forge feedstock.

Version Bump and Dependency Updates
===================================

1. **Update the Version Number**

   - Open a pull request (PR) to update the version number in the following
     two files:

     - ``mpas_analysis/version.py``

     - ``ci/recipe/meta.yaml``

   - Make sure the version follows `semantic versioning <https://semver.org/>`_.

2. **Check and Update Dependencies**

   - Ensure that dependencies and their constraints are up-to-date and
     consistent in:

     - ``ci/recipe/meta.yaml`` (dependencies for the conda-forge release)

     - ``pyproject.toml`` (dependencies for PyPI; used as a sanity check)

     - ``dev-spec.txt`` (development dependencies; should be a superset of
       those for the conda-forge release)

   - The dependencies in ``meta.yaml`` are the ones that will be used for the
     released package on conda-forge. The dependencies in ``pyproject.toml``
     are for PyPI and should be kept in sync as much as possible but are only
     there as a sanity check when we run ``pip check``. The ``dev-spec.txt``
     file should include all dependencies needed for development and testing.

   - Review and update dependency versions and constraints as needed.

3. **Make a PR and merge it**

Tagging and Publishing a Release
================================

4. **Tag the Release on GitHub**

   - Go to https://github.com/MPAS-Dev/MPAS-Analysis/releases and click on
     "Draft a new release".

   - Enter the appropriate tag for the release, following semantic versioning
     (e.g., ``1.13.0``; **do not** include a ``v`` in front).

   - Enter a release title (typically the release version **with** a ``v`` in
     front, e.g., ``v1.13.0``).

   - Write a description and/or use the "Generate release notes" button to
     auto-generate release notes.

   - If the release is ready, click "Publish release". Otherwise, save it as a
     draft.

Updating the conda-forge Feedstock
==================================

5. **Update the conda-forge Feedstock**

   - After the release is published, update and merge a PR for the new release
     at the conda-forge feedstock:
     https://github.com/conda-forge/mpas-analysis-feedstock

   - The conda-forge bot should automatically create a pull request for the
     new version within a few hours to a day after the release.

   - Compare the dependencies in the new release to those in the previous
     release and update the recipe as needed. To do this:

     - Find the most recent release at
       https://github.com/MPAS-Dev/MPAS-Analysis/releases

     - Use the "Compare" feature to select the previous release.

     - Under "changed files", locate ``ci/recipe/meta.yaml`` to see
       any dependency changes.

   - Review and update the feedstock PR as needed, then merge it.

   - If you are not already a maintainer of the feedstock, you can request to
     be added by creating a new issue at
     https://github.com/conda-forge/mpas-analysis-feedstock/issues, choosing
     "Bot command", and putting
     ``@conda-forge-admin, please add user @username`` as the subject.
