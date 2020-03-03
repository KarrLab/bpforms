Comparison of `BpForms` with other formats, databases, and software
-------------------------------------------------------------------

Comparison with `BigSMILES <https://doi.org/10.1021/acscentsci.9b00476>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `BpForms` is backward compatible with the IUPAC/IUBMB format
* `BpForms` supports high-level naming of residues, crosslinks, nicks, and circularity
* With `BcForms`, `BpForms` can be compared into yet more abstract descriptions of complexes
* `BpForms` can capture missing or uncertain knowledge about molecules
* There are software tools for working with `BpForms`


Comparison with BioPAX
^^^^^^^^^^^^^^^^^^^^^^

* `BpForms` can represent non-canonical DNA, RNA, and proteins
* `BpForms` can represent any modification and, therefore, is not limited to modifications that have been previously enumerated in the `PSI Molecular Interation ontology <https://www.ebi.ac.uk/ols/ontologies/mi>`. This is necessary to represent the combinatorial complexity of non-canonical DNA, RNA, and proteins.
* `BpForms` concretely represents the primary structure of biopolymers. In constrast, the chemical structure implied by BioPAX features and chemical reactions is unclear.
* `BpForms` is easier to embed into other files such as SBML-encoded models.
* `BpForms` includes software for interpreting descriptions of non-canonical biopolymers. This enables calculations of properties such as chemical formulae, masses, and charges, which is essential for modeling and other applications.


Comparison with DNAmod
^^^^^^^^^^^^^^^^^^^^^^

* The `BpForms` DNA alphabet is internally consistent. Each monomeric form represents a nucleotide monophosphate. This ensures the monomeric forms can be composed into polymers, as the semantic meaning of the polymers is well-defined.

Comparison with `HELM <https://www.pistoiaalliance.org/helm-project/>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* The `BpForms` grammar is fully backward compatible with the IUPAC/IUBMB format. While the HELM format is similar to the IUPAC/IUBMB format, the HELM representation of an unmodified sequence is not equivalent to its IUPAC/IUBMB representation.
* `BpForms` supports inline definitions of residues. This enables users to describe molecules that involve new residues without having to create their own alphabet or edit the public alphabets. In contrast, HELM has a more hierarchical design that requires users to describe all new residues in a custom alphabet. We believe this is more cumbersome, particular for use cases where users want to embed descriptions of molecules into other files (e.g., Excel workbooks or PDF documents) without having to deal with references to files that define alphabets.
* `BpForms` provides a high-level representation of crosslinks. HELM does not support this.
* `BpForms` can capture missing or uncertain information about molecules, which is essential for biological research. HELM does not support this.
* With `BcForms`, `BpForms` can be composed into yet higher-level descriptions of complexes. In contrast, HELM can only capture complexes whose subunits are linked by crosslinks. `BcForms` is more flexible, and can represent a complex as a bag of subunits. This flexibility is essential for representing complexes that do not involve crosslinks or for which crosslink information is not available.

Comparison with MODOMICS
^^^^^^^^^^^^^^^^^^^^^^^^

* `BpForms` can represent DNA, RNA, and proteins including the diversity of ids of monomeric forms used by MODOMICS, RESID, and other databases.
* `BpForms` can represent any modification, and is not limited to the modifications catalogued in MODOMICS.
* `BpForms` concretely captures the bonds between adjacent monomeric forms, avoiding ambiguity about how monomeric forms are composed into polymers. This is particularly important for monomeric forms that have multiple 3' and 5' sites.
* `BpForms` can represent bonds between non-adjacent monomeric forms, such as crosslinks.
* `BpForms` can represent circular biopolymers.
* `BpForms` can capture uncertainity in the structures of biopolymers. This is essential for proteomics.
* `BpForms` has a concrete grammar.
* All of the monomeric forms in the `BpForms` RNA alphabet have defined structures, unlike several MODOMICS entries which have undefined `BASE` groups. This gaurantees that `BpForms` polymers specify concrete structures.
* The `BpForms` RNA alphabet consistenty represents only nucleotide monophosphates. This gaurantees that the monomeric forms are composable. In constract, MODOMICS includes both nucleosides and bases, which makes the composition of the MODOMICS monomeric forms ill-defined.
* The `BpForms` RNA alphabet encompasses monomeric forms from both MODOMICS and the RNA Modification Database. This adds two additional monomeric forms that are not present in MODOMICS.
* `BpForms` includes software for error checking descriptions of non-canonical biopolymers. This includes verifying that monomeric forms that only have left bonding sites (e.g. 5' caps) only appear at the first position and that monomeric forms that only have right bonding sites (e.g. 3' caps) only appear at the last position.
* `BpForms` includes software for interpreting descriptions of non-canonical biopolymers. This enables calculations of properties such as chemical formulae, masses, and charges, which is essential for modeling and other applications.


Comparison with the `PDB Chemical Component Dictionary/Ligand <http://www.rcsb.org/pdb/ligand/chemAdvSearch.do>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* The `BpForms` protein alphabet includes numerous additional monomeric forms from RESID.
* The `BpForms` protein alphabet consistently represents the C termini as COH. In constrast, the PDB CCD represents C termini as both COOH and COH.
* The composition of `BpForms` monomeric forms into sequences is well defined. The `BpForms` protein alphabet explicitly describes the location of the C and N termini of each monomeric form. This enables the `BpForms` software to verify that BpForms describe valid atomic structures. The semantics for combining PDB CCD monomeric forms into sequences is unclear. Consequently, sequences of PDB CCD monomeric forms cannot be validated. In particular, monomeric forms with multiple C and N termini have ambiguous bonding becuase the PDB CCD monomeric forms do not capture the left and right bonding sites.


Comparison with `PDB format <http://www.wwpdb.org/documentation/file-format>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `BpForms` is a more compact, human readable description of the primary structure of biopolymers.
* It is easier to compose monomeric forms into `BpForms` than the PDB format.
* `BpForms` is easier to embed into other standards and formats.


Comparison with `PDB SEQRES annotations <http://www.wwpdb.org/documentation/file-format>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `BpForms` can capture any monomeric form, including monomeric forms which are not part of the PDB Chemical Component Dictionary.
* `BpForms` concretely captures the bonds between successive monomeric forms. PDB SEQRES annotations are ambiguous for monomeric forms with multiple C and N termini.
* `BpForms` can capture circularity.
* `BpForms` can capture crosslinks.
* There is no defined semantics for generating atomic structures from PDB SEQRES annotations.
* The `BpForms` software can verify that `BpForms` describe valid atomic structures. PDB SEQRES annotations cannot be verified because there is no defined semantics for generating atomic structures from these annotations.

Comparison with `ProForma Proteoform Notation <http://www.topdownproteomics.org/resources/proforma/>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `BpForms` can represent DNA, RNA, and proteins.
* `BpForms` can represent any modification and, therefore, is not limited to modifications that have been previously enumerated in databases and ontologies. This is necessary to represent the combinatorial complexity of non-canonical DNA, RNA, and proteins.
* `BpForms` concretely captures the bonds between adjacent monomeric forms, avoiding ambiguity about how monomeric forms are composed into polymers. This is particularly important for monomeric forms that have multiple C and N termini, which affects numerous entries in RESID.
* `BpForms` can represent monomeric forms which can only bond to the right and left or which don't have backbones such as 3' and 5' caps.
* `BpForms` can represent bonds between non-adjacent monomeric forms, such as disulfide bonds.
* `BpForms` can represent circular biopolymers.
* `BpForms` separates the representation of non-canonical biopolymers from the chemical processes which generate them.
* `BpForms` can capture additional uncertainty in the structures of biopolymers: uncertainty in the position of a non-canonical monomeric form within a sequence, and uncertainty in the chemical identity of a non-canonical monomeric form (e.g., deviation from its expected mass or charge).
* `BpForms` has a concrete grammar.
* `BpForms` includes software for error checking descriptions of non-canonical biopolymers.
* `BpForms` includes software for interpreting descriptions of non-canonical biopolymers. This enables calculations of properties such as chemical formulae, masses, and charges, which is essential for modeling and other applications.


Comparison with RESID
^^^^^^^^^^^^^^^^^^^^^

* Each monomeric form in the `BpForms` protein alphabet has a defined structure. This gaurantees that polymers have well-defined structures. In constrast, RESID has numerous entires without defined structures.
* The composability of the monomeric forms in the `BpForms` protein alphabet is well-defined. Each form has at most one left-bonding-terminus (C) and at most one right-bonding-terminus (N). This eliminates confusion about the meaning of composition monomeric forms with multiple N and C-termini. In contrast, RESID has numerous entries with multiple N or C-termini whose composition into polymers is ill-defined.
* The `BpForms` protein alphabet encompasses entries from additional databases.


Comparison with the RNA Modification Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Each monomeric form in the `BpForms` protein alphabet has a machine-readable structure. This gaurantees that polymers have well-defined structures. In constrast, RESID has numerous entires without defined structures. In contrast, the RNA Modification Database only provides images and CAS ids, neither or which can easily be converted into SMILES.


Comparison with the Synthetic Biology Open Language (SBOL)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `BpForms` concretely captures the primary structure of non-canonical biopolymers. In particular, it concretely captures the covalent bonds between monomeric forms. In contrast, SBOL's sequence annotations capture insufficient information to define the primary structure of a non-canonical biopolymer. The chemical meaning of these sequence annotations are unclear.
* `BpForms` directly captures the primary structure of biopolymers. In constract, SBOL indirectly captures structures via the reactions that produce them via operations such as cutting.
* `BpForms` can capture uncertainity in the structures of biopolymers. This is essential for proteomics.
* `BpForms` is easier to embed into other files such as SBML-encoded models.
* `BpForms` includes software for interpreting descriptions of non-canonical biopolymers. This enables calculations of properties such as chemical formulae, masses, and charges, which is essential for modeling and other applications.


Comparison with the `World-wide Monomer Reference Database <http://www.monomer.org/#/main>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* The `BpForms` alphabets include many more residues.
* The residues in the `BpForms` alphabets include synonyms, comments, and unification links with other databases.
