Test Suite Infrastructure
=========================

The `suite` directory provides a comprehensive infrastructure for testing
MPAS-Analysis on supported machines (Anvil, Chrysalis, Perlmutter-CPU, and
Compy). The suite is designed to ensure code changes do not introduce
unexpected results and to validate MPAS-Analysis in various environments.

Overview of Test Scripts
------------------------

The main entry point is ``suite/run_suite.bash``.  It supports three modes:

1. **Developer Testing**: ``./suite/run_suite.bash --dev``

   - This is the recommended workflow for development in a Pixi environment.

   - Run it either from an active Pixi shell or with an explicit Pixi
     environment name:

     .. code-block:: bash

        $ pixi shell
        $ ./suite/run_suite.bash --dev

     or:

     .. code-block:: bash

        $ ./suite/run_suite.bash --dev --pixi-env py313

   - It builds the documentation, renders the suite configs, and submits the
     suite jobs using ``pixi run`` in the selected environment.

   - Each task produces a web page with results, accessible via the web portal.

   - After completion, check for successful web page generation, e.g.:

     .. code-block:: bash

        $ tail -n 3 chrysalis_test_suite/main_py3.13/mpas_analysis.o793058

   - To quickly identify unfinished or failed tasks:

     .. code-block:: bash

        $ grep -L "Web page:" chrysalis_test_suite/*/mpas_analysis.o*

   - Developers should run this suite manually on each pull request before
     merging and link the results in the PR.

2. **Package Build & Test**: ``./suite/run_suite.bash``

   - This mode builds the MPAS-Analysis conda package and tests it in fresh
     environments.

   - It creates conda environments for multiple Python versions, runs tests,
     builds documentation, and executes the analysis suite.

   - Recommended for more thorough validation, especially before releases.

3. **E3SM-Unified Deployment Testing**:
   ``./suite/run_suite.bash --e3sm-unified``

   - This mode is used during test deployments of E3SM-Unified to verify
     MPAS-Analysis works as expected within the deployment.

   - It is typically run by E3SM-Unified maintainers during deployment
     testing.

Supported Machines
------------------

The suite is designed to run only on supported machines:

- Anvil

- Chrysalis

- Perlmutter-CPU (`pm-cpu`)

- Compy

If you attempt to run the suite on an unsupported machine, you will receive an
error.

Modifying the Test Suite
------------------------

Developers may need to update the suite for new requirements:

- **Python Versions**:

  - The Python versions tested in package mode are defined at the top of
    ``suite/run_suite.bash`` (for example ``main_py=3.13`` and
    ``alt_py=3.12``).

  - To test additional versions, add them to the relevant script variables and
    loops.

- **Adding New Machines**:

  - Update the machine detection logic in `suite/setup.py` and add appropriate
    input/output paths for the new machine.

  - Ensure the new machine is supported in the scripts and the web portal
    configuration.

- **Adding/Modifying Tests**:

  - To add new tests, update the run lists in ``suite/run_suite.bash`` and
    provide corresponding config files in ``suite/configs``.

  - New tests could change which analysis tasks are run, the configuration for
    running tasks overall (e.g. how climatologies are computed), or how
    individual tasks are configured (e.g. focused on polar regions vs. global)

- **Changing Simulation Data**:

  - Update the simulation name and mesh in `suite/setup.py` if you wish to
    test on different output.

Best Practices
--------------

- Always run the test suite before merging a pull request.

- Link the results web page in your PR for reviewers.

- Use the quick check (`grep -L "Web page:" ...`) to ensure all tasks
  completed successfully.

- Update the suite scripts and configs as needed to keep pace with
  MPAS-Analysis development.

The suite templates live in ``suite/templates`` and the run-specific config
overrides live in ``suite/configs``.
