Installation
============

The following is a brief guide to installing `BpForms`. The `Dockerfile <https://github.com/KarrLab/bpforms/blob/master/Dockerfile>`_ in the `BpForms` Git repository contains detailed instructions for how to install `BpForms` in Ubuntu Linux.

Prerequisites
--------------------------

First, install the third-party packages listed below.

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
Run the following command to install the latest release from PyPI.::

    pip install bpforms

Latest revision from GitHub
---------------------------
Run the following command to install the latest version from GitHub.::

    pip install git+https://github.com/KarrLab/pkg_utils.git#egg=pkg_utils
    pip install git+https://github.com/KarrLab/wc_utils.git#egg=wc_utils[chem]
    pip install git+https://github.com/KarrLab/bpforms.git#egg=bpforms


Installing the optional features
--------------------------------
To calculate major protonation and tautomerization states, `BpForms` must be installed with the `[protontation]` option::

    pip install bpforms[protontation]
    pip install git+https://github.com/KarrLab/bpforms.git#egg=bpforms[protontation]

To draw molecules, `BpForms` must be installed with the `[draw]` option::

    pip install bpforms[draw]
    pip install git+https://github.com/KarrLab/bpforms.git#egg=bpforms[draw]

To export the alphabets in OBO format, `BpForms` must be installed with the `[onto_export]` option::

    pip install bpforms[onto_export]
    pip install git+https://github.com/KarrLab/bpforms.git#egg=bpforms[onto_export]

To install the rest API, `BpForms` must be installed with the `[rest_api]` option::

    pip install bpforms[rest_api]
    pip install git+https://github.com/KarrLab/bpforms.git#egg=bpforms[rest_api]
