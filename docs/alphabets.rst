Alphabets
------------------

`BpForms` comes with six alphabets:

* Canonical DNA: The four canonical DNA nucleobases.
* Canonical RNA: The four canonical RNA nucleosides.
* Canonical protein: The 20 canonical protein residues.
* DNA: The four canonical DNA nucleobases, plus the non-canonical DNA nucleobases in `DNAmod <https://dnamod.hoffmanlab.org>`_.
* RNA: The four canonical RNA nucleosides, plus the non-canonical RNA nucleosides in `MODOMICS <http://modomics.genesilico.pl/modifications/>`_.
* Protein: The 20 canonical protein residues, plus the non-canonical protein residues in `RESID <https://pir.georgetown.edu/resid/>`_.

To support compatibility with the entries in these and other databases, `BpForms` represents biopolymers as sequence of ``monomeric forms`` linked together via a ``backbone``. For compatibility with DNAmod, the monomeric forms form the DNA alphabets are nucleobases and the backbones for the DNA alphabets is deoxyribose 5-phosphate. For compatibility with MODOMICS, the monomeric forms of the RNA alphabets are nucleosides and the backbones for the RNA alphabets is hydrogen phosphate. For compatibility with RESID, the monomeric forms of the protein alphabets are amino acid stubs and the backbones of the protein alphabets is hydroxide.

Building additional alphabets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Users can also create additional alphabets by creating additional instances of ``bpforms.Alphabet``. Users can either add monomeric forms programmatically or load them from YAML files. Users can add monomeric forms programmatically by creating instances of ``bpforms.Monomer`` and adding them to the ``monomers`` attribute of ``bpforms.Alphabet``, which is a dictionary that maps the character codes of monomeric forms to monomeric forms. Users can load monomeric forms from YAML files by using the ``from_yaml`` method of ``bpforms.Alphabet``. Please see `bpforms/alphabet/dna.yml <https://github.com/KarrLab/bpforms/blob/master/bpforms/alphabet/dna.yml>`_ for an example of the YAML alphabet format.


Contributing to the alphabets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We welcome contributions to the alphabets, as well as new alphabets. See :numref:`contributing` for information about how to contribute to `BpForms`.
