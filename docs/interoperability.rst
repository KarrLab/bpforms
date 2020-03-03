Interoperating between ontologies and with other formats
--------------------------------------------------------

Interoperating between user alphabets of residues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The goal of `BpForms` is to faciliate precise communication about polymers. From a technical perspective, the easiest way to achieve this would be to require all users to use the same alphabets, require each character in each alphabet to represent a unique molecular structure, and not allow users to define residues inline or to define their own alphabets. With only one set of alphabets, it would be easy to compare descriptions of molecules and generate high-level descriptions of any differences. However, we believe that `BpForms` should enable users to define their own residues and even entire alphabets because (a) even the simplest process for contributing residues to the public alphabets could be a barrier for some users, and we prefer users to use `BpForms` rather than not describe polymers precisely at all and (b) some fields may prefer to manage their own alphabets (without a large number of users yet, we think it may be difficult to influence every field -- structural biology, transcriptomics, proteomics, systems biology, synthetic biology -- to converage on one set of alphabets).

Consequently, there will likely be a need to interoperate between alphabets. Because `BpForms` describes concrete molecular structures (i.e. atoms and bonds), `BpForms` can compare molecules described with different alphabets and convert descriptions of molecules between alphabets using the molecular meanings of the descriptions of molecules. Below are our plans for facilitating such interoperability:

* Expand the grammar and data structures to capture the alphabet (including user alphabets) used to describe each molecule.
* Implement a method for comparing molecules described with different alphabets by generating the molecular structures for the molecules and comparing these structures.
* Implement a method for migrating descriptions of molecules between alphabets by (i) using the molecular structures of the residues in the alphabets to build a mapping from the residues of the source alphabet to the residues of the target alphabet, (ii) using this mapping to migrate all mapped residues, and (iii) replacing all unmapped residues with user-defined residues.


Interoperating between user ontologies of crosslinks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Interoperability between ontologies of crosslinks can be managed similar to the interoperating between alphabets of residues as discussed above.

Interoperating with other formats
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Similarly, `BpForms` can interoperate (e.g., compare, calculate differences) with any other molecularly-precise format, such as the InChI, SMILES, and BigSMILES, through converting each format to a common molecular representation.
