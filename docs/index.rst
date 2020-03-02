`BpForms` documentation
=======================

`BpForms` is a set of tools for concretely representing the primary structures of non-canonical forms of biopolymers, such as oxidized DNA, methylated RNA, and acetylated proteins, and calculating properties of non-canonical biopolymers.

`BpForms` emcompasses five tools:

* A grammar for concretely describing the primary structures of non-canonical biopolymers. See :numref:`grammar` for detailed information. For example, the following represents a modified DNA molecule that contains a deoxyinosine monomeric form at the fourth position:
  ::

    ACG{dI}

  This concrete representation enables the `BpForms` software tools to calculate properties of non-canonical biopolymers.

* Tools for calculating properties of non-canonical biopolymers including their chemical formulae, molecular weights, charges, and major protonation and tautomerization states.

  * A web-based graphical interface: See `https://bpforms.org <https://bpforms.org>`_ and :numref:`graphical_web_interface`.
  * A JSON REST API: See `https://bpforms.org/api <https://bpforms.org/api>`_ and :numref:`rest_api`.
  * A command line interface: See :numref:`cli`.
  * A Python API: See :numref:`python_api`.

`BpForms` was motivated by the need to concretely represent the biochemistry of DNA modification, DNA repair, post-transcriptional processing, and post-translational processing in `whole-cell computational models <https://www.wholecell.org>`_. `BpForms` is also a valuable tool for experimental proteomics and synthetic biology. In particular, we developed `BpForms` because there were no notations, schemas, data models, or file formats for concretely representing modified forms of biopolymers, despite the existence of several databases and ontologies of DNA, RNA, and protein modifications, the `ProForma Proteoform Notation <https://www.topdownproteomics.org/resources/proforma/>`_, and the `MOMODICS <http://modomics.genesilico.pl>`_ codes for modified RNA bases.

`BpForms` can be combined with `BcForms <https://www.bcforms.org>`_ to concretely describe the primary structure of complexes.

Contents
--------

.. toctree::
   :maxdepth: 3
   :numbered:

   use_cases.rst
   installation.rst
   tutorial.rst
   alphabets.rst
   crosslinks.rst
   integrations.rst
   resources.rst
   comparison.rst
   future.rst
   contributing.rst
   API documentation <source/modules.rst>
   about.rst
