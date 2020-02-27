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

    pip install bpforms[all]

Latest revision from GitHub
---------------------------
Run the following command to install the latest version from GitHub.::

    pip install git+https://github.com/KarrLab/pkg_utils.git#egg=pkg_utils[all]
    pip install git+https://github.com/KarrLab/wc_utils.git#egg=wc_utils[all]
    pip install git+https://github.com/KarrLab/bpforms.git#egg=bpforms[all]
