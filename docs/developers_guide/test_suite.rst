Test Suite Infrastructure
=========================

The `suite` directory provides a comprehensive infrastructure for testing
MPAS-Analysis on supported machines (Anvil, Chrysalis, Perlmutter-CPU, and
Compy). The suite is designed to ensure code changes do not introduce
unexpected results and to validate MPAS-Analysis in various environments.

Overview of Test Scripts
------------------------

There are three main scripts for running the test suite:

1. **run_dev_suite.bash** (Developer Testing)

   - Use this script after activating your development environment
     (must be named `mpas_analysis_dev`).

   - It builds the documentation and runs a series of analysis tasks on output
     from a low-resolution (QUwLI240) simulation.

   - Each task produces a web page with results, accessible via the web portal.

   - Example usage:

     .. code-block:: bash

        $ source ~/miniforge3/etc/profile.d/conda.sh
        $ conda activate mpas_analysis_dev
        $ ./suite/run_dev_suite.bash

   - After completion, check for successful web page generation, e.g.:

     .. code-block:: bash

        $ tail -n 3 chrysalis_test_suite/main_py3.11/mpas_analysis.o793058

     The last lines should include:

     .. code-block:: none

        Generating webpage for viewing results...
        Web page: https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/<username>/analysis_testing/chrysalis/<branch>/main_py3.11/

   - To quickly identify unfinished or failed tasks:

     .. code-block:: bash

        $ grep -L "Web page:" chrysalis_test_suite/*/mpas_analysis.o*

   - Developers should run this suite manually on each pull request before
     merging and link the results in the PR.

2. **run_suite.bash** (Package Build & Test)

   - Use this script to build the MPAS-Analysis conda package and test it in
     fresh environments.

   - It creates conda environments for multiple Python versions, runs tests,
     builds documentation, and executes the analysis suite.

   - Recommended for more thorough validation, especially before releases.

   - Example usage:

     .. code-block:: bash

        $ ./suite/run_suite.bash

3. **run_e3sm_unified_suite.bash** (E3SM-Unified Deployment Testing)

   - Used during test deployments of E3SM-Unified to verify MPAS-Analysis
     works as expected within the deployment.

   - Typically run by E3SM-Unified maintainers during deployment testing.

   - Example usage:

     .. code-block:: bash

        $ ./suite/run_e3sm_unified_suite.bash

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

  - The Python versions tested are defined in the scripts (e.g.,
    `main_py=3.11`, `alt_py=3.10`).

  - To test additional versions, add them to the relevant script variables and
    loops.

- **Adding New Machines**:

  - Update the machine detection logic in `suite/setup.py` and add appropriate
    input/output paths for the new machine.

  - Ensure the new machine is supported in the scripts and the web portal
    configuration.

- **Adding/Modifying Tests**:

  - To add new tests, update the list of runs in the scripts and
    provide corresponding config files in the `suite` directory.

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

For more details, see the comments and documentation within each script and
config file in the `suite` directory.
