Building the Documentation
==========================

With the ``mpas_analysis_dev`` environment activated, you can run:

.. code-block:: bash

    cd docs
    DOCS_VERSION=test make clean versioned-html

to build the docs locally in the ``_build/html`` subdirectory.

The docs should build cleanly.  If they don't, please attempt to fix the
errors and warnings even if they are not related to your changes.  We want
to keep the documentation in good shape.

Previewing the Documentation
----------------------------

When generating documentation on HPC machines, you will want to copy the html
output to the public web space to view it, or if the web portal is being
cranky, scp it to your local machine.

To preview the documentation locally, open the ``index.html`` file in the
``_build/html/test`` directory with your browser or try:

.. code-block:: bash

    cd _build/html
    python -m http.server 8000

Then, open http://0.0.0.0:8000/test/ in your browser.
