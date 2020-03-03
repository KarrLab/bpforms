Managing revisions to the grammar and the ontologies and user ontologies
------------------------------------------------------------------------

Over time, we will likely introduce new capabilties into the grammar (e.g., to represent additional types of uncertainty) and expand (and, if necessary, correct) the ontologies. In addition, we anticipate that users will develop their own ontologies of residues and crosslinks. Consequently, there will likely be needs to (a) compare molecules that are described with different versions of the grammar, different versions of the ontologies, and different ontologies and (b) convert descriptions of molecules between versions of the grammar, versions ontologies, and different ontologies.

Comparing descriptions of molecules described with different versions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Because each combination of (a) a version of the grammar and (b) a set of versions of the ontologies will be capable of representing the molecular structures (atoms and bonds) of polymers, it will be possible to compare molecules described with different versions of the grammar, different versions of the ontologies, and different ontologies at the molecular level.

Below are our plans for faciliating such comparisons:

* Separate the alphabets of residues and ontology of crosslinks into their own Git repositories so their revisions can be separately tracked from that of the `BpForms` software.
* Create a repository for user ontologies. This could be implemented as a Git repository.
* Expand the grammar and data structures to capture the ontologies (including user ontologies) used to describe each molecule.
* Expand the grammar and data structures to capture the versions of the grammar and the versions of the ontologies used to describe each molecule.
* Implement a method for comparing molecules describe with different versions of the grammar, different ontologies, or different versions of the ontologies. This method will use the Git repositories to resolve the specific versions of the grammar and ontologies.

Migrating descriptions of molecules between different ontologies and versions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To make it easier for researchers to work with revisions to the grammar and ontologies, it will also be useful to develop a utility for migrating descriptions of molecules between different versions of the grammar and/or different ontologies and versions of the ontologies. Below are our plans for facilitating such migrations.

* Implement data structures for capturing changes to the grammar and ontologies.
* Use these data structures to implement a method for migrating descriptions of molecules between successive versions of the grammar or ontologies.
* Implement a second method for migrating descriptions of molecules between arbitrary versions of the grammar and ontologies by chaining together multiple executions of the first method.
