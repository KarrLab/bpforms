[![PyPI package](https://img.shields.io/pypi/v/bpforms.svg)](https://pypi.python.org/pypi/bpforms)
[![Documentation](https://readthedocs.org/projects/bpforms/badge/?version=latest)](https://docs.karrlab.org/bpforms)
[![Test results](https://circleci.com/gh/KarrLab/bpforms.svg?style=shield)](https://circleci.com/gh/KarrLab/bpforms)
[![Test coverage](https://coveralls.io/repos/github/KarrLab/bpforms/badge.svg)](https://coveralls.io/github/KarrLab/bpforms)
[![Code analysis](https://api.codeclimate.com/v1/badges/e35081f676dfbb5ac46f/maintainability)](https://codeclimate.com/github/KarrLab/bpforms)
[![License](https://img.shields.io/github/license/KarrLab/bpforms.svg)](LICENSE)
![Analytics](https://ga-beacon.appspot.com/UA-86759801-1/bpforms/README.md?pixel)

# `BpForms`: toolkit for concretely describing non-canonical DNA, RNA, and proteins

`BpForms` is a set of tools for concretely representing the primary structures of non-canonical forms of biopolymers, such as oxidized DNA, methylated RNA, and acetylated proteins, and calculating properties of non-canonical biopolymers.

`BpForms` encompasses five tools:

* A grammar for concretely describing the primary structures of non-canonical biopolymers. See the [documentation](https://docs.karrlab.org/bpforms/) for more information. For example, the following text represents a modified DNA molecule that contains a deoxyinosine monomeric form at the fourth position.
  ```
  ACG[id: "dI"
       | structure: "[H][C@]1(O)C[C@@]([H])(O[C@]1([H])CO)N1C=NC2=C1N=CN=C2O"]T
  ```

  This concrete representation enables the `BpForms` software tools to calculate properties of non-canonical biopolymers.

* Tools for calculating properties of non-canonical biopolymers including their chemical formulae, molecular weights, charges, and major protonation and tautomerization states.
  * A web app: [https://bpforms.org](https://bpforms.org)
  * A JSON REST API: [https://bpforms.org/api](https://bpforms.org/api)
  * A command line interface. See the [documentation](https://docs.karrlab.org/bpforms/master/0.0.1/cli.html) for more information.
  * A Python API. See the [documentation](https://docs.karrlab.org/bpforms/master/0.0.1/python_api.html) for more information.

`BpForms` was motivated by the need to concretely represent the biochemistry of DNA modification, DNA repair, post-transcriptional processing, and post-translational processing in [whole-cell computational models](https://www.wholecell.org). `BpForms` is also a valuable tool for experimental proteomics and synthetic biology. In particular, we developed `BpForms` because there were no notations, schemas, data models, or file formats for concretely representing non-canonical forms of biopolymers, despite the existence of several databases and ontologies of DNA, RNA, and protein modifications, the [ProForma Proteoform Notation](https://www.topdownproteomics.org/resources/proforma/), and the [MOMODICS](http://modomics.genesilico.pl/) codes for modified RNA bases.

*BpForms* can be combined with [*BcForms*](https://www.bcforms.org) to concretely describe the primary structure of complexes.

## Installation
1. Install the third-party dependencies listed below. Detailed installation instructions are available in [An Introduction to Whole-Cell Modeling](http://docs.karrlab.org/intro_to_wc_modeling/master/0.0.1/installation.html).

    * [ChemAxon Marvin](https://chemaxon.com/products/marvin): optional to calculate major protonation and tautomerization states
      * [Java](https://www.java.com) >= 1.8
    * [Open Babel](http://openbabel.org)
    * [Pip](https://pip.pypa.io) >= 19.0
    * [Python](https://www.python.org) >= 3.6

2. To use Marvin to calculate major protonation and tautomerization states, set ``JAVA_HOME`` to the path to your Java virtual machine (JVM)
   ```
   export JAVA_HOME=/usr/lib/jvm/default-java
   ```

3. To use Marvin to calculate major protonation and tautomerization states, add Marvin to the Java class path
   ```
   export CLASSPATH=$CLASSPATH:/opt/chemaxon/marvinsuite/lib/MarvinBeans.jar
   ```

4. Install this package

    * Install the latest release from PyPI:
      ```
      pip install bpforms
      ```

    * Install the latest revision from GitHub:
      ```
      pip install git+https://github.com/KarrLab/wc_utils.git#egg=wc_utils[all]
      pip install git+https://github.com/KarrLab/bpforms.git#egg=bpforms
      ```

    * To install the rest API, `BpForms` must be installed with the `[all]` option:
      ```
      pip install bpforms[all]
      pip install git+https://github.com/KarrLab/bpforms.git#egg=bpforms[all]
      ```

## Examples, tutorial, and documentation
Please see the [documentation](https://docs.karrlab.org/bpforms). An [interactive tutorial](https://sandbox.karrlab.org/notebooks/bpforms/Tutorial.ipynb) is also available in the whole-cell modeling sandbox.

## License
The package is released under the [MIT license](LICENSE).

## Citing `BpForms`
Lang PF, Chebaro Y & Jonathan R. Karr. BpForms: a toolkit for concretely describing modified DNA, RNA and proteins. arXiv:1903.10042. [:link:](https://arxiv.org/abs/1903.10042)

## Development team
This package was developed by the [Karr Lab](https://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai in New York, USA.

* [Jonathan Karr](https://www.karrlab.org)
* [Yassmine Chebaro](https://www.linkedin.com/in/yassmine-chebaro-6bb8a05/)
* [Paul Lang](http://www.dtc.ox.ac.uk/people/17/langp/)

## Questions and comments
Please contact the [Karr Lab](mailto:info@karrlab.org) with any questions or comments.
