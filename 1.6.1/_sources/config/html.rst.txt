.. _config_html:

HTML
====

The ``[html]`` section simply specifies whether or not a webpage should be
generated for displaying the plots produced by the analysis::

  [html]
  ## options related to generating a webpage to display the analysis

  # generate the webpage?
  generate = True

The webpage is produced in the directory specified by ``htmlSubdirectory``
in the ``[output]`` section, see :ref:`config_output`.
