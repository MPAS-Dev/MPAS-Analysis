.. _dev_overview:

Overview
========

MPAS-Analysis is a `python package <https://docs.python.org/3/tutorial/modules.html#packages>`_.
All of the code in the package can be accessed in one of two ways.  The first
is the command-line interface with the ``mpas_analysis`` command, see
:ref:`dev_mpas_analysis_command`.  The second way is through import commands
like:

.. code-block:: python

    from mpas_analysis.shared.io import NameList


    namelist = NameList('namelist.ocean')

Before we dig into the details of how to develop new analysis tasks and other
infrastructure for MPAS-Analysis, we first give a little bit of background on
the design philosophy behind the package.

.. _dev_style:

Code Style
----------

As new code gets added, we are aiming to adhere fairly strictly to the
`PEP8 style guide <https://www.python.org/dev/peps/pep-0008/>`_.  Older code
was not written with this requirement so a lot of the existing code will not
follow PEP8.  As more code begins to use the PEP8 conventions, we may add a
bot to flag any PEP8 violations (similar to `compass <https://github.com/MPAS-Dev/compass>`_)
as part of each pull request.  Please consider using an editor that
automatically flags PEP8 violations during code development, such as
`pycharm <https://www.jetbrains.com/pycharm/>`_ or
`spyder <https://www.spyder-ide.org/>`_, or a linter, such as
`flake8 <https://flake8.pycqa.org/en/latest/>`_ or
`pep8 <https://pep8.readthedocs.io/>`_.  We discourage you from automatically
reformatting your code (e.g. with `autopep8 <https://github.com/hhatto/autopep8>`_)
because this can often produce undesirable and confusing results.

The `flake8 <https://flake8.pycqa.org/en/latest/>`_ utility for linting python
files to the PEP8 standard is included in the :ref:`dev_environment`. To use
flake8, just run ``flake8`` from any directory and it will return lint results
for all files recursively through all subdirectories.  You can also run it for a
single file or using wildcards (e.g., ``flake8 *.py``).  There also is a
`vim plugin <https://github.com/nvie/vim-flake8>`_ that runs the flake8 linter
from within vim.  If you are not using an IDE that lints automatically, it is
recommended you run flake8 from the command line or the vim plugin before
committing your code changes.

.. _dev_environment:

conda environment for development
---------------------------------

To develop the code, you will first need to clone the repo from
`https://github.com/MPAS-Dev/MPAS-Analysis <https://github.com/MPAS-Dev/MPAS-Analysis>`_
and/or add your own fork as a "remote".

Then, you will need to set up a conda environment from the MPAS-Analysis repo.
This environment will include the required dependencies for the development
branch defined in the file ``dev-spec.txt``.  The ``mpas_analysis`` package
will be into the conda environment in "edit mode", so that the the environment
uses the version directly out of the local branch and your changes are seen as
you make them.  First, we recommend some one-time setup to make sure you
have the ``conda-forge`` channel set as the highest priority and you have the
``mamba`` package:

.. code-block:: bash

   conda config --add channels conda-forge
   conda config --set channel_priority strict
   conda update -y --all
   conda install mamba

To create a development environment for the current branch of MPAS-Analysis,
run:

.. code-block:: bash

   mamba create -y -n mpas_dev --file dev-spec.txt
   conda activate mpas_dev
   python -m pip install -e .

.. _dev_environment_tricks:

tricks for developing 2 packages at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are developing another conda package at the same time (this is common
for `MPAS-Tools <https://github.com/MPAS-Dev/MPAS-Tools>`_ or
`geometric_features <https://github.com/MPAS-Dev/geometric_features>`_).
You can install both packages in "edit mode" in the same development
environment, e.g.:

.. code-block:: bash

   mamba create -y -n mpas_dev --file tools/MPAS-Tools/conda_package/dev-spec.txt \
       --file analysis/MPAS-Analysis/dev-spec.txt
   conda activate mpas_dev
   cd ~/code/tools/MPAS-Tools/conda_package
   python -m pip install -e .
   cd ~/code/analysis/MPAS-Analysis
   python -m pip install -e .

Obviously, the paths to the repos may be different in your local clones.  With
the ``mpas_dev`` environment as defined above, you can make changes to both
MPAS-Tools and MPAS-Analysis packages in their respective branches, and
these changes will be reflected when refer to the packages or call their
respective entry points (command-line tools).

One warning with MPAS-Tools is that c++ tools are not compiled or installed
using this approach.  Any modifications you make to those tools will go
unnoticed by this approach.