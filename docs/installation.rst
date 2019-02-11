Installation
============

Prerequisites
--------------------------

First, install the third-party packages listed below. Detailed installation instructions are available in `An Introduction to Whole-Cell Modeling <http://docs.karrlab.org/intro_to_wc_modeling/master/0.0.1/installation.html>`_.

* `ChemAxon Marvin <https://chemaxon.com/products/marvin>`_
* `Open Babel <http://openbabel.org>`_
* `Pip <https://pip.pypa.io>`_ >= 18.0
* `Python <https://www.python.org>`_ >= 3.6

Latest release From PyPI
---------------------------
Run the following command to install the latest release from PyPI. For most environments, the ``--process-dependency-links`` option is needed to install some of the dependencies from GitHub.::

    pip install --process-dependency-links bpforms

Latest revision from GitHub
---------------------------
Run the following command to install the latest version from GitHub. For most environments, the ``--process-dependency-links`` option is needed to install some of the dependencies from GitHub.::

    pip install --process-dependency-links git+git://github.com/KarrLab/bpforms.git#egg=bpforms
