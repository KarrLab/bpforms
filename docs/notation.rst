.. _notation:

`BpForms` notation
------------------

The `BpForms` notation unambiguously represents the primary structure of biopolymer forms that contain canonical and non-canonical monomers using (a) a syntax similar to FASTA and (b) extended alphabets for DNA, RNA, and proteins based on DNAmod, MODOMICS, and RESID. This enables `BpForms` to calculate the formula, molecular, weight, and charge of biopolymer forms.

* Canonical monomers are indicated by their single upper case character codes.
* Non-canonical monomers defined in the alphabets that have single-character codes are indicated by these codes.
* Non-canonical monomers defined in the alphabets that have multiple-character codes are indicated by these codes delimited by curly brackets.
* Non-canonical monomers that are not defined in the alphabet can be defined "inline" with multiple attributes separated by vertical pipes ("|") enclosed inside square brackets. The structures of these monomers are defined in SMILES format using the ``structure`` attribute. Additional attributes can provide metadata about monomers such as their ids and names.

`BpForms` contains three pre-built canonical DNA, RNA and protein alphabets and three extended DNA, RNA, and protein alphabets based on DNAmod, MODOMICS, and RESID. Users can also create additional custom alphabets.

Examples:

* ``[id: "dI" | name: "deoxyinosine"]ACGC``: represents deoxyinosine at the first position
* ``AC[id: "dI" | name: "deoxyinosine"]GC``: represents deoxyinosine at the third position
* ``ACGC{6A}``: represents methyl-6-adenosine at the last position


Structure
^^^^^^^^^

The ``structure`` attribute describes the chemical structure of the monomer as a SMILES-encoded string::

    [id: "dI" |
        structure: "O=C1NC=NC2=C1N=CN2"
        ]

The ``monomer-bond-atom``, ``monomer-displaced-atom``, ``left-bond-atom``, ``left-displaced-atom``, ``right-bond-atom``, and ``right-displaced-atom`` attributes describe the linkages between the monomer and the backbone and between successive monomers::
    
    [id: "dI" |
        structure: "O=C1NC=NC2=C1N=CN2"
        ]

We recommend defining these attributes for each monomer. These attributes must be defined to calculate the structure, formula, molecular weight, and charge of the biopolymer.


Uncertainty
^^^^^^^^^^^

`BpForms` can also represent two types of uncertainty in the structures of biopolymer forms.

* The ``delta-mass`` and ``delta-charge`` attributes can describe uncertainty in the chemical identities of monomers. For example, ``[id: "dAMP" | delta-mass: 1 | delta-charge: 1]`` indicates the presence of an additional proton exact location is not known.
* The ``position`` attribute can describe uncertainty in the position of a monomer. For example, ``[id: "5mC" | position: 2-3]`` indicates that 5mC may occur anywhere between the second and third position.


Metadata
^^^^^^^^

`BpForms` can also represent several types of metadata:

* The ``id`` and ``name`` attributes provide human-readable labels for monomers. Only one id and one name is allowed per monomer::

    [id: "dI"
        | name: "deoxyinosine"
        ]

* The ``synonym`` attribute provides additional human-readable labels. Each monomer can have multiple synonyms::

    [id: "dI"
        | synonym: "2'-deoxyinosine"
        | synonym: "2'-deoxyinosine, 9-[(2R,4S,5R)-4-hydroxy-5-(hydroxymethyl)tetrahydrofuran-2-yl]-9H-purin-6-ol"
        ]

* The ``identifier`` attribute describes references to entries in external databases. Each monomer can have multiple identifiers. The namespaces and ids of identifers must be separated by "/"::

    [id: "dI"
        | identifier: "biocyc.compound" / "DEOXYINOSINE"
        | identifier: "chebi" / "CHEBI:28997"
        | identifier: "pubchem.compound" / "65058"
        ]

* The ``base-monomer`` attribute describes other monomer(s) which the monomer is generated from. The value of this attribute must be the code of a monomer in the alphabet. Each monomer can have one or more bases. This annotation is needed to generate more informative FASTA sequences for `BpForms`::

    [id: "m2A"
        | base-monomer: "A"
        ]

* The ``comments`` attribute describes additional information about each monomer. Each monomer can only have one comment::

    [id: "dI"
        | comments: "A purine 2'-deoxyribonucleoside that is inosine in which the
                     hydroxy group at position 2' is replaced by a hydrogen."
        ]


Syntax
^^^^^^

* Monomers that are in the alphabet are indicated by a single character or multiple characters delimiated by curly brackets.
* Monomers that are not in the alphabet are defined "inline" with one or more attributes separated by vertical pipes ("|") inside square brackets.

  * All of the attributes are optional. However, the `structure` attribute is required to compute the formula, molecular weight, and charge of the biopolymer.
  * Attributes are separated by vertical pipes ("|").
  * Attributes and their values are separated by colons (":").
  * White spaces are ignored.
  * The values of the `id`, `name`, `synonym`, and `comments` attributes must be enclosed in quotes ('"').
  * The namespace and id of each identifer must be separated by "/".

* The positions of the monomers in the string indicates in their location in the sequence.


Grammar
^^^^^^^

The following is the definition of the `BpForms` grammar. The grammar is defined in `Lark syntax <https://lark-parser.readthedocs.io/en/latest/grammar/>`_ which is based on `EBNF syntax <https://en.wikipedia.org/wiki/Extended_Backus%E2%80%93Naur_form>`_.

.. literalinclude:: ../bpforms/grammar.lark
    :language: text


Examples
^^^^^^^^

* DNA::
    
    ACGT[id: "dI" | structure: "O=C1NC=NC2=C1N=CN2" | monomer-bond-atom: Monomer / N / 10 / 0 | monomer-displaced-atom: Monomer / H / 10 / 0]AG{m2A}

* RNA:: 

    {01G}CGU[id: "01A" | structure: "COC1C(O)C(OC1n1cnc2c1ncn(c2=N)C)CO"]AG[id: "019A" | structure: "COC1C(O)C(OC1n1cnc2c1ncn(c2=O)C)CO"]

* Protein::

    ARGKL[id: "AA0318" | structure: "COC(=O)[C@@H]([NH3+])CCCC[NH3+]"]YRCG[id: "AA0567" | structure: "CC=CC(=O)NCCCC[C@@H](C=O)[NH3+]"]
