Alphabets
------------------

`BpForms` comes with six alphabets:

* Canonical DNA: The four canonical DNA nucleobases.
* Canonical RNA: The four canonical RNA nucleosides.
* Canonical protein: The 20 canonical protein residues.
* DNA: The four canonical DNA nucleobases, plus the non-canonical DNA nucleobases in `DNAmod <https://dnamod.hoffmanlab.org>`_.
* RNA: The four canonical RNA nucleosides, plus the non-canonical RNA nucleosides in `MODOMICS <http://modomics.genesilico.pl/modifications/>`_.
* Protein: The 20 canonical protein residues, plus the non-canonical protein residues in `RESID <https://pir.georgetown.edu/resid/>`_.

To support compatibility with the entries in these and other databases, `BpForms` represents biopolymers as so called ``monomers`` linked together via a ``Backbone``. For DNAmod based alphabets, ``monomers`` are nucleobases, and the ``Backbone`` is deoxyribose 5-phosphate. For MODOMICS based alphabets ``monomers`` are nucleosides and the ``Backbone`` is hydrogen phosphate. For RESID based alphabets, ``monomers`` are amino acids, and the ``Backbone`` is hydroxide.


Building additional alphabets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Users can also create additional alphabets by creating additional instances of ``bpforms.Alphabet``. Users can either add monomers programmatically or load them from YAML files. Users can add monomers programmatically by creating instances of ``bpforms.Monomer`` and adding them to the ``monomers`` attribute of ``bpforms.Alphabet``, which is a dictionary that maps the character codes of monomers to monomers. Users can load monomers from YAML files by using the ``from_yaml`` method of ``bpforms.Alphabet``. Please see `bpforms/alphabet/dna.yml <https://github.com/KarrLab/bpforms/blob/master/bpforms/alphabet/dna.yml>`_ for an example of the YAML alphabet format.


Contributing to the alphabets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We welcome contributions to the alphabets, as well as new alphabets. See :numref:`contributing` for information about how to contribute to `BpForms`.
