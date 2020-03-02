Ontology of crosslinks
----------------------
As of March 2020, the ontology of crosslinks contains 36 crosslinks between protein residues.

We developed the ontology based on entries in `RESID <https://pir.georgetown.edu/resid/>`_ which represent crosslinked dipeptides:

#. We searched RESID for entries that represent crosslinked dipeptides. 
#. We identified the individual residues which participate in each dimer. 
#. We used `ChemAxon Marvin <https://chemaxon.com/products/marvin>`_ to identify the atoms involved in each crosslink. 
#. We used `Open Babel <https://pir.georgetown.edu/resid/>`_ to determine the indices of these atoms in the canonical SMILES ordering of the atoms. 
#. We manually assigned an id and name to each crosslink. 

We invite the community to help us curate additional crosslinks, including DNA-DNA, RNA-RNA, DNA-protein, and RNA-protein crosslinks. Please contact us at `info@karrlab.org <mailto:info@karrlab.org>`_ to get involved, or please submit additional crosslinks via GitHub pull requests or issues.
