Alphabets of residues
---------------------

`BpForms` comes with six alphabets:

* Canonical DNA: The canonical DNA nucleotide monophosphates.
* Canonical RNA: The canonical RNA nucleotide monophosphates.
* Canonical protein: The 20 canonical protein residues.
* DNA: The canonical DNA nucleotide monophosphates, plus non-canonical DNA nucleotide monophosphates based on `DNAmod <https://dnamod.hoffmanlab.org>`_.
* RNA: The canonical RNA nucleotide monophosphates, plus non-canonical RNA nucleotide monophosphates based on `MODOMICS <http://modomics.genesilico.pl/modifications/>`_ and the `RNA Modification Database <https://mods.rna.albany.edu/mods/>`_.
* Protein: The 20 canonical protein residues, plus the non-canonical protein residues in the `PDB Chemical Component Dictionary <http://www.wwpdb.org/data/ccd>`_ and `RESID <https://pir.georgetown.edu/resid/>`_.


Construction of the public alphabets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The alphabets of DNA, RNA, and protein residues were developed by merging residues from multiple databases:

* The DNA alphabet was developed by combining the deoxyribose nucleotide monophosphates and 3' and 5' DNA ends from the `PDB Chemical Component Dictionary <http://www.wwpdb.org/data/ccd>`_ (PDB CCD) with the verified DNA nucleobases from `DNAmod <https://dnamod.hoffmanlab.org/>`_ and the deoxyribose nucleosides from `REPAIRtoire <http://repairtoire.genesilico.pl/damage/>`_ that had concrete structures. 
* The RNA alphabet was developed by combining the ribose nucleotide monophosphates and 3' and 5' RNA ends from the PDB CCD with the ribose nucleosides from `MODOMICS <http://modomics.genesilico.pl/modifications/>`_ and the `RNA Modification Database <https://mods.rna.albany.edu/>`_ that had concrete structures. 
* The protein alphabet was developed by merging residues and ends from the PDB CCD and `RESID <https://pir.georgetown.edu/resid/>`_. 

Specifically, the alphabets were constructed as follows:

#. We downloaded, scraped, and manually extracted residues from DNAmod, MODOMICS, the PDB CCD, REPAIRtoire, RESID, the RNA Modification Database. 
#. We parsed each database into a list of residues. 
#. We rejected residues with incompletely defined structures, as well as inconsistent residues such as nucleotides from DNAmod. 
#. We normalized the DNA and RNA residues to nucleotide monophosphates and normalized the protein residues to amino acids. For example, we transformed the DNAmod entries to nucleotides by adding deoxyribose monophosphate to each nucleobase. 
#. We merged the repeated residues. This included residues that had the same molecular structure, that the upstream sources annotated were equivalent, or that had similar names. 
#. We identified the atom indices of the left/preceding and right/following bonding sites in each residue. For DNA and RNA, the left and right bonding sites comprised the phosphorus atom in the phosphate group bonded to the 5' carbon and the oxygen atom bonded to the 3' carbon, respectively. For protein residues, these comprised the nitrogen atom in the amino group and the acidic oxygen atom in the carboxyl group.

The above steps are implemented by :py:func:`bpforms.util.build_alphabets`.


Updating the public alphabets
"""""""""""""""""""""""""""""
:py:func:`bpforms.util.build_alphabets` can be used to pull updates to DNAmod, MODOMICS, the PDB CCD, REPAIRtoire, RESID, and the RNA Modification Database into the alphabets. We plan to use this function to update the alphabets annually or as needed.


Building additional alphabets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Users can also create additional alphabets by creating additional instances of ``bpforms.Alphabet``. Users can either add monomeric forms programmatically or load them from YAML files. Users can add monomeric forms programmatically by creating instances of ``bpforms.Monomer`` and adding them to the ``monomers`` attribute of ``bpforms.Alphabet``, which is a dictionary that maps the character codes of monomeric forms to monomeric forms. Users can load monomeric forms from YAML files by using the ``from_yaml`` method of ``bpforms.Alphabet``. Please see `bpforms/alphabet/dna.yml <https://github.com/KarrLab/bpforms/blob/master/bpforms/alphabet/dna.yml>`_ for an example of the YAML alphabet format.


Contributing to the alphabets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We welcome contributions to the alphabets, as well as new alphabets. See :numref:`contributing` for information about how to contribute to `BpForms`.
