Installation
============

Prerequisites
--------------------------

First, install the third-party packages listed below. Detailed installation instructions are available in `An Introduction to Whole-Cell Modeling <http://docs.karrlab.org/intro_to_wc_modeling/master/0.0.1/installation.html>`_.

* `ChemAxon Marvin <https://chemaxon.com/products/marvin>`_: optional to calculate major protonation and tautomerization states

    * `Java <https://www.java.com>`_ >= 1.8

* `Open Babel <http://openbabel.org>`_
* `Pip <https://pip.pypa.io>`_ >= 18.0
* `Python <https://www.python.org>`_ >= 3.6

To use ChemAxon Marvin to calculate major protonation and tautomerization states, set ``JAVA_HOME`` to the path to your Java virtual machine (JVM) and add Marvin to the Java class path::

   export JAVA_HOME=/usr/lib/jvm/default-java
   export CLASSPATH=$CLASSPATH:/opt/chemaxon/marvinsuite/lib/MarvinBeans.jar

Latest release From PyPI
---------------------------
Run the following command to install the latest release from PyPI. For most environments, the ``--process-dependency-links`` option is needed to install some of the dependencies from GitHub.::

    pip install --process-dependency-links bpforms[all]

Latest revision from GitHub
---------------------------
Run the following command to install the latest version from GitHub. For most environments, the ``--process-dependency-links`` option is needed to install some of the dependencies from GitHub.::

    pip install --process-dependency-links git+git://github.com/KarrLab/bpforms.git#egg=bpforms[all]
