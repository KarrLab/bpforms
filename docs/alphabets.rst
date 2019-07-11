Alphabets
------------------

`BpForms` comes with six alphabets:

* Canonical DNA: The canonical DNA nucleotide monophosphates.
* Canonical RNA: The canonical RNA nucleotide monophosphates.
* Canonical protein: The 20 canonical protein residues.
* DNA: The canonical DNA nucleotide monophosphates, plus non-canonical DNA nucleotide monophosphates based on `DNAmod <https://dnamod.hoffmanlab.org>`_.
* RNA: The canonical RNA nucleotide monophosphates, plus non-canonical RNA nucleotide monophosphates based on `MODOMICS <http://modomics.genesilico.pl/modifications/>`_ and the `RNA Modification Database <https://mods.rna.albany.edu/mods/>`_.
* Protein: The 20 canonical protein residues, plus the non-canonical protein residues in the `PDB Chemical Component Dictionary <http://www.wwpdb.org/data/ccd>`_ and `RESID <https://pir.georgetown.edu/resid/>`_.

Building additional alphabets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Users can also create additional alphabets by creating additional instances of ``bpforms.Alphabet``. Users can either add monomeric forms programmatically or load them from YAML files. Users can add monomeric forms programmatically by creating instances of ``bpforms.Monomer`` and adding them to the ``monomers`` attribute of ``bpforms.Alphabet``, which is a dictionary that maps the character codes of monomeric forms to monomeric forms. Users can load monomeric forms from YAML files by using the ``from_yaml`` method of ``bpforms.Alphabet``. Please see `bpforms/alphabet/dna.yml <https://github.com/KarrLab/bpforms/blob/master/bpforms/alphabet/dna.yml>`_ for an example of the YAML alphabet format.


Contributing to the alphabets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We welcome contributions to the alphabets, as well as new alphabets. See :numref:`contributing` for information about how to contribute to `BpForms`.
