.. _dev_releasing:

***********************
Releasing a New Version
***********************

This document describes the steps for maintainers to tag and release a new
version of ``MPAS-Analysis``, and to update the conda-forge feedstock.

Version Bump and Dependency Updates
===================================

1. **Update the Version Number**

   - Manually update the version number in the following files:

     - ``mpas_analysis/version.py``
     - ``ci/recipe/meta.yaml``

   - Make sure the version follows `semantic versioning <https://semver.org/>`_.
     For release candidates, use versions like ``1.3.0rc1``.

2. **Check and Update Dependencies**

   - Ensure that dependencies and their constraints are up-to-date and
     consistent in:

     - ``ci/recipe/meta.yaml`` (for the conda-forge release)

     - ``pyproject.toml`` (for PyPI; used for sanity checks)

     - ``dev-spec.txt`` (development dependencies; should be a superset)

   - Use the GitHub "Compare" feature to check for dependency changes between releases:
     https://github.com/MPAS-Dev/MPAS-Analysis/compare

3. **Make a PR and merge it**

   - Open a PR for the version bump and dependency changes and merge once
     approved.

Tagging and Publishing a Release
================================

4. **Tagging a Release Candidate**

   - For release candidates, **do not create a GitHub release page**. Just
     create a tag from the command line:

     - Make sure your changes are merged into ``develop`` or your own update
       branch (e.g. ``update-to-1.3.0``) and your local repo is up to date.

     - Tag the release candidate (e.g., ``1.3.0rc1``):

       ::

           git checkout develop
           git fetch --all -p
           git reset --hard origin/develop
           git tag 1.3.0rc1
           git push origin 1.3.0rc1

       (Replace ``1.3.0rc1`` with your actual version, and ``develop`` with
       your branch if needed.)

     **Note:** This will only create a tag. No release page will be created
     on GitHub.

5. **Publishing a Stable Release**

   - For stable releases, create a GitHub release page as follows:

     - Go to https://github.com/MPAS-Dev/pyremap/releases

     - Click "Draft a new release"

     - Enter a tag (e.g., ``1.3.0``)

     - Set the release title to the version prefixed with ``v`` (e.g.,
       ``v1.3.0``)

     - Generate or manually write release notes

     - Click "Publish release"

Updating the conda-forge Feedstock
==================================

6. **Automatic Feedstock Update (Preferred Method)**

   - Wait for the ``regro-cf-autotick-bot`` to open a PR at:
     https://github.com/conda-forge/mpas-analysis-feedstock

   - This may take several hours to a day.

   - Review the PR:
     - Confirm the version bump and dependency changes
     - Merge once CI checks pass

   **Note:** If you are impatient, you can accellerate this process by creating
   a bot issue at: https://github.com/conda-forge/mpas-analysis-feedstock/issues
   with the subject ``@conda-forge-admin, please update version``.  This
   will open a new PR with the version within a few minutes.

7. **Manual Feedstock Update (Fallback Method)**

   If the bot PR does not appear or is too slow, update manually:

   - Download the release tarball:

     ::

         wget https://github.com/MPAS-Dev/MPAS-Analysis/archive/refs/tags/<version>.tar.gz

   - Compute the SHA256 checksum:

     ::

         shasum -a 256 <version>.tar.gz

   - In the ``meta.yaml`` of the feedstock recipe:
     - Set ``{% set version = "<version>" %}``
     - Set the new ``sha256`` value
     - Update dependencies if needed

   - Commit, push to a new branch, and open a PR against the feedstock
   - Follow any instructions in the PR template and merge once approved

Post Release Actions
====================

8. **Verify and Announce**

   - Install the package in a clean environment to test:

     ::

         conda create -n test-mpas -c conda-forge mpas-analysis=<version>

   - Optionally announce the release on relevant communication channels

   - Update any documentation or release notes as needed
