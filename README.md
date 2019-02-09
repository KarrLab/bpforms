[![PyPI package](https://img.shields.io/pypi/v/bpforms.svg)](https://pypi.python.org/pypi/bpforms)
[![Documentation](https://readthedocs.org/projects/bpforms/badge/?version=latest)](https://docs.karrlab.org/bpforms)
[![Test results](https://circleci.com/gh/KarrLab/bpforms.svg?style=shield)](https://circleci.com/gh/KarrLab/bpforms)
[![Test coverage](https://coveralls.io/repos/github/KarrLab/bpforms/badge.svg)](https://coveralls.io/github/KarrLab/bpforms)
[![Code analysis](https://api.codeclimate.com/v1/badges/e35081f676dfbb5ac46f/maintainability)](https://codeclimate.com/github/KarrLab/bpforms)
[![License](https://img.shields.io/github/license/KarrLab/bpforms.svg)](LICENSE)
![Analytics](https://ga-beacon.appspot.com/UA-86759801-1/bpforms/README.md?pixel)

# BpForms: unambiguous representation of modified DNA, RNA, and proteins

BpForms is a set of tools for unambiguously representing the structures of modified forms of biopolymers such as DNA, RNA, and protein. 

* The BpForms notation can unambiguously represent the structure of modified forms of biopolymers. For example, the following represents a modified DNA molecule that contains a deoxyinosine residue at the fourth position.
  ```
  ACG[id: "dI" 
       | structure: InChI=1S
          /C10H12N4O4
          /c15-2-6-5(16)1-7(18-6)14-4-13-8-9(14)11-3-12-10(8)17
          /h3-7,15-16H,1-2H2,(H,11,12,17)
          /t5-,6+,7+
          /m0
          /s1
          ]T
  ```
* This concrete representation of modified biopolymers enables the BpForms software tools to calculate the formulae, molecular weights, and charges of biopolymers, as well as automatically protonate biopolymers for specific pHs.

BpForms emcompasses five tools:

* [Notation for describing biopolymers](https://docs.karrlab.org/bpforms/)
* Web-based graphical interface: [https://bpforms.org](https://bpforms.org)
* [REST JSON API](https://docs.karrlab.org/bpforms/)
* [Command line interface](https://docs.karrlab.org/bpforms/)
* [Python API](https://docs.karrlab.org/bpforms/)

BpForms was motivated by the need to concretely represent the biochemistry of DNA modification, DNA repair, post-transcriptional processing, and post-translational processing in [whole-cell computational models](https://www.wholecell.org). In addition, BpForms are a valuable tool for experimental proteomics. In particular, we developed BpForms because there were no notations, schemas, data models, or file formats for concretely representing modified forms of biopolymers, despite the existence of several databases and ontologies of DNA, RNA, and protein modifications and the [ProForma Proteoform Notation](https://www.topdownproteomics.org/resources/proforma/).

The BpForms syntax was inspired by the ProForma Proteoform Notation. BpForms improves upon this syntax in several ways:
 
* BpForms separates the representation of modified biopolymers from the chemical processes which generate them. 
* BpForms clarifies the representation of multiply modified monomers. This is necessary to represent the combinatorial complexity of modified DNA, RNA, and proteins.
* BpForms can represent any modification and, therefore, is not limited to previously enumerated modifications. This is also necessary to represent the combinatorial complexity of modified DNA, RNA, and proteins.
* BpForms supports two additional types of uncertainty in the structures of biopolyers: uncertainty in the positions of modifications and uncertainty in the charges of modifications.
* BpForms has a concrete grammar. This enables error checking, as well the calculation of formulae, masses, and charges which is essential for modeling.

## Installation
1. Install the third-party dependencies listed below. Detailed installation instructions are available in [An Introduction to Whole-Cell Modeling](http://docs.karrlab.org/intro_to_wc_modeling/master/0.0.1/installation.html).
        
    * [ChemAxon Marvin](https://chemaxon.com/products/marvin)
    * [Open Babel](http://openbabel.org)
    * [Pip](https://pip.pypa.io) >= 18.0
    * [Python](https://www.python.org) >= 3.6

2. Install this package 

    * Install the latest release from PyPI
      ```
      pip install bpforms
      ```

    * Install the latest revision from GitHub
      ```
      pip install git+git://github.com/KarrLab/bpforms#egg=bpforms
      ```

## Examples, tutorial, and documentation
Please see the [documentation](https://docs.karrlab.org/bpforms).

## License
The package is released under the [MIT license](LICENSE).

## Development team
This package was developed by the [Karr Lab](https://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai in New York, USA.

* [Jonathan Karr](https://www.karrlab.org)
* [Yassmine Chebaro](https://www.linkedin.com/in/yassmine-chebaro-6bb8a05/)
* [Paul Lang](http://www.dtc.ox.ac.uk/people/17/langp/)

## Questions and comments
Please contact the [Karr Lab](https://www.karrlab.org) with any questions or comments.
