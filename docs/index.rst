`BpForms` documentation
=======================


`BpForms` is a set of tools for unambiguously representing the structures of modified forms of biopolymers such as DNA, RNA, and protein. 

* The `BpForms` notation can unambiguously represent the structure of modified forms of biopolymers. For example, the following represents a modified DNA molecule that contains a deoxyinosine base at the fourth position.::
  
    ACG[
        id: "dI" | structure: InChI=1S
            /C10H12N4O4
            /c15-2-6-5(16)1-7(18-6)14-4-13-8-9(14)11-3-12-10(8)17
            /h3-7,15-16H,1-2H2,(H,11,12,17)
            /t5-,6+,7+
            /m0
            /s1
        ]T

* This concrete representation of modified biopolymers enables the `BpForms` software tools to calculate the chemical formulae, molecular weights, and charges of biopolymers, as well as automatically protonate biopolymers for specific pHs.

`BpForms` emcompasses five tools:

* Notation for describing biopolymers: See :numref:`notation`.
* Web-based graphical interface: See `https://bpforms.org <https://bpforms.org>`_ and :numref:`graphical_web_interface`.
* REST JSON API: See :numref:`rest_api`.
* Command line interface: See :numref:`cli`.
* Python API: See :numref:`python_api`.

`BpForms` was motivated by the need to concretely represent the biochemistry of DNA modification, DNA repair, post-transcriptional processing, and post-translational processing in `whole-cell computational models <https://www.wholecell.org>`_. In addition, `BpForms` are a valuable tool for experimental proteomics. In particular, we developed `BpForms` because there were no notations, schemas, data models, or file formats for concretely representing modified forms of biopolymers, despite the existence of several databases and ontologies of DNA, RNA, and protein modifications and the `ProForma Proteoform Notation <https://www.topdownproteomics.org/resources/proforma/>`_.

The `BpForms` syntax was inspired by the ProForma Proteoform Notation. `BpForms` improves upon this syntax in several ways:
 
* `BpForms` separates the representation of modified biopolymers from the chemical processes which generate them. 
* `BpForms` clarifies the representation of multiply modified monomers. This is necessary to represent the combinatorial complexity of modified DNA, RNA, and proteins.
* `BpForms` can be customized to represent any modification and, therefore, is not limited to previously enumerated modifications. This is also necessary to represent the combinatorial complexity of modified DNA, RNA, and proteins.
* `BpForms` supports two additional types of uncertainty in the structures of biopolyers: uncertainty in the position of a modified nucleotide/amino acid and uncertainty in its charge.
* `BpForms` has a concrete grammar. This enables error checking, as well the calculation of chemical formulae, masses, and charges which is essential for modeling.

Contents
--------

.. toctree::
   :maxdepth: 3
   :numbered:

   installation.rst
   tutorial.rst
   alphabets.rst
   resources.rst
   API documentation <source/modules.rst>
   about.rst
