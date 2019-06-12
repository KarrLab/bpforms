Advantages of `BpForms` over previous formats, databases, and software tools
----------------------------------------------------------------------------

Advantages over BioPAX
^^^^^^^^^^^^^^^^^^^^^^

* `BpForms` can represent non-canonical DNA, RNA, and proteins
* `BpForms` can represent any modification and, therefore, is not limited to modifications that have been previously enumerated in the `PSI Molecular Interation ontology <https://www.ebi.ac.uk/ols/ontologies/mi>`. This is necessary to represent the combinatorial complexity of non-canonical DNA, RNA, and proteins.
* `BpForms` concretely represents the primary structure of biopolymers. In constrast, the chemical structure implied by BioPAX features and chemical reactions is unclear.
* `BpForms` is easier to embed into other files such as SBML-encoded models.
* `BpForms` includes software for interpreting descriptions of non-canonical biopolymers. This enables calculations of properties such as chemical formulae, masses, and charges, which is essential for modeling and other applications.


Advantages over DNAmod
^^^^^^^^^^^^^^^^^^^^^^

* The `BpForms` DNA alphabet is internally consistent. Each monomeric form represents a nucleobase. This ensures the monomeric forms can be composed into polymers, as the semantic meaning of the polymers is well-defined.


Advantages over MODOMICS
^^^^^^^^^^^^^^^^^^^^^^^^

* `BpForms` can represent DNA, RNA, and proteins including the diversity of ids of monomeric forms used by MODOMICS, RESID, and other databases.
* `BpForms` can represent any modification, and is not limited to the modifications catalogued in MODOMICS.
* `BpForms` can represent bonds between non-adjacent monomeric forms, such as crosslinks.
* `BpForms` can represent circular biopolymers.
* `BpForms` can capture uncertainity in the structures of biopolymers. This is essential for proteomics.
* `BpForms` has a concrete grammar.
* All of the monomeric forms in the `BpForms` RNA alphabet have defined structures, unlike several MODOMICS entries which have undefined `BASE` groups. This gaurantees that `BpForms` polymers specify concrete structures.
* The `BpForms` RNA alphabet consistenty represents only nucleosides. This gaurantees that the monomeric forms are composable. In constract, MODOMICS includes both nucleosides and bases, which makes the composition of the MODOMICS monomeric forms ill-defined.
* The `BpForms` RNA alphabet encompasses monomeric forms from both MODOMICS and the RNA Modification Database. This adds two additional monomeric forms that are not present in MODOMICS.
* `BpForms` includes software for error checking descriptions of non-canonical biopolymers. This includes verifying that monomeric forms that only have left bonding sites (e.g. 5' caps) only appear at the first position and that monomeric forms that only have right bonding sites (e.g. 3' caps) only appear at the last position.
* `BpForms` includes software for interpreting descriptions of non-canonical biopolymers. This enables calculations of properties such as chemical formulae, masses, and charges, which is essential for modeling and other applications.


Advantages over `ProForma Proteoform Notation <http://www.topdownproteomics.org/resources/proforma/>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `BpForms` can represent DNA, RNA, and proteins.
* `BpForms` can represent any modification and, therefore, is not limited to modifications that have been previously enumerated in databases and ontologies. This is necessary to represent the combinatorial complexity of non-canonical DNA, RNA, and proteins.
* `BpForms` can represent monomeric forms which can only bind to the right and left or which don't have backbones such as 3' and 5' caps.
* `BpForms` can represent bonds between non-adjacent monomeric forms, such as disulfide bonds.
* `BpForms` can represent circular biopolymers.
* `BpForms` separates the representation of non-canonical biopolymers from the chemical processes which generate them.
* `BpForms` can capture additional uncertainty in the structures of biopolymers: uncertainty in the position of a non-canonical monomeric form within a sequence, and uncertainty in the chemical identity of a non-canonical monomeric form (e.g., deviation from its expected mass or charge).
* `BpForms` has a concrete grammar.
* `BpForms` includes software for error checking descriptions of non-canonical biopolymers.
* `BpForms` includes software for interpreting descriptions of non-canonical biopolymers. This enables calculations of properties such as chemical formulae, masses, and charges, which is essential for modeling and other applications.


Advantages over RESID
^^^^^^^^^^^^^^^^^^^^^

* Each monomeric form in the `BpForms` protein alphabet has a defined structure. This gaurantees that polymers have well-defined structures. In constrast, RESID has numerous entires without defined structures.
* The composability of the monomeric forms in the `BpForms` protein alphabet is well-defined. Each form has at most one left-binding-terminus (C) and at most one right-binding-terminus (N). This eliminates confusion about the meaning of composition monomeric forms with multiple N and C-termini. In contrast, RESID has numerous entries with multiple N or C-termini whose composition into polymers is ill-defined.
* The `BpForms` protein alphabet encompasses entries from additional databases.


Advantages over the RNA Modification Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Each monomeric form in the `BpForms` protein alphabet has a machine-readable structure. This gaurantees that polymers have well-defined structures. In constrast, RESID has numerous entires without defined structures. In contrast, the RNA Modification Database only provides images and CAS ids, neither or which can easily be converted into SMILES.


Advantages over the Synthetic Biology Open Language (SBOL)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `BpForms` concretely captures the primary structure of non-canonical biopolymers. In particular, it concretely captures the covalent bonds between monomeric forms. In contrast, SBOL's sequence annotations capture insufficient information to define the primary structure of a non-canonical biopolymer. The chemical meaning of these sequence annotations are unclear.
* `BpForms` directly captures the primary structure of biopolymers. In constract, SBOL indirectly captures structures via the reactions that produce them via operations such as cutting.
* `BpForms` can capture uncertainity in the structures of biopolymers. This is essential for proteomics.
* `BpForms` is easier to embed into other files such as SBML-encoded models.
* `BpForms` includes software for interpreting descriptions of non-canonical biopolymers. This enables calculations of properties such as chemical formulae, masses, and charges, which is essential for modeling and other applications.
