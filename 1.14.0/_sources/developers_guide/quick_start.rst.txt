Quick Start for Developers
==========================

This guide provides a condensed overview for developers to get started with
MPAS-Analysis development.

1. Fork and Clone the Repository
--------------------------------
- Fork `MPAS-Analysis <https://github.com/MPAS-Dev/MPAS-Analysis>`_ on GitHub.
- Clone the main repo and your fork locally:
  - Create a base directory (e.g., ``mpas-analysis``).
  - Clone the main repo:

    .. code-block:: bash

       git clone git@github.com:MPAS-Dev/MPAS-Analysis.git develop

  - Add your fork as a remote:

    .. code-block:: bash

       git remote add <username>/MPAS-Analysis git@github.com:<username>/MPAS-Analysis.git

2. Configure Git
----------------
- Set up your ``~/.gitconfig`` with your name and email (must match your
  GitHub account).
- Recommended: set editor, color, and useful aliases.

3. Set Up SSH Keys
------------------
- Add SSH keys to GitHub for push access.
- See: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account

4. Create a Development Worktree
--------------------------------
- Fetch latest changes:

    .. code-block:: bash

       git fetch --all -p

- Create a worktree for your feature branch:

    .. code-block:: bash

       git worktree add ../<feature_branch>

- Enter the worktree directory:

    .. code-block:: bash

       cd ../<feature_branch>

5. Set Up Conda Environment
---------------------------
- Install Miniforge3 (recommended) or Miniconda.
- For Miniconda, add ``conda-forge`` channel and set strict priority.
- Create environment:

    .. code-block:: bash

       conda create -y -n mpas_analysis_dev --file dev-spec.txt

- Activate:

    .. code-block:: bash

       conda activate mpas_analysis_dev

- Install MPAS-Analysis in edit mode:

    .. code-block:: bash

       python -m pip install --no-deps --no-build-isolation -e .

6. Activate Environment (each session)
--------------------------------------
- For bash:

    .. code-block:: bash

       source ~/miniforge3/etc/profile.d/conda.sh; conda activate mpas_analysis_dev

- For csh:

    .. code-block:: csh

       source ~/miniforge3/etc/profile.d/conda.csh; conda activate mpas_analysis_dev

7. Configure and Run MPAS-Analysis
----------------------------------
- Copy and edit a config file (e.g., ``example_e3sm.cfg``) for your run.
- Set required options: ``mainRunName``, ``baseDirectory``, ``mpasMeshName``, output paths.
- Set ``mapMpiTasks = 1`` and ``mapParallelExec = None`` for development environments.
- Export HDF5 file locking variable if needed:
  - Bash:

    .. code-block:: bash

       export HDF5_USE_FILE_LOCKING=FALSE

  - Csh:

    .. code-block:: csh

       setenv HDF5_USE_FILE_LOCKING FALSE

- Run analysis:

    .. code-block:: bash

       mpas_analysis -m <machine> <your_config>.cfg

8. View Results
---------------
- Output is a set of web pages in your specified output directory.
- On some systems, update permissions:

    .. code-block:: bash

       chmod -R ugo+rX <output_dir>

- See the main web page for links to results and provenance info.

Additional Recommendations
--------------------------
- Use VS Code for remote editing and linting (optional).

For more details, see the full :doc:`../tutorials/dev_getting_started`.
