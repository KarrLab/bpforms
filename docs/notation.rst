.. _notation:

`BpForms` notation
------------------

The `BpForms` notation unambiguously represents the primary structure of biopolymer forms that contain canonical and non-canonical monomers using (a) a syntax similar to FASTA and (b) extended alphabets for DNA, RNA, and proteins to describe monomers. `BpForms` contains three pre-built canonical DNA, RNA and protein alphabets and three extended DNA, RNA, and protein alphabets curated from the DNAmod, MODOMICS, and RESID databases, respectively. Users can also create additional custom alphabets. These alphabets are associated with their corresponding biopolyer form, which allows `BpForms` to calculate the chemical structure (e.g. in SMILES format), chemical formula, molecular weight, and charge of a biopolymer.

* Monomers that are present in the alphabet are indicated by a single character or multiple characters delimited by curly brackets.
* Monomers that are not in the alphabet are defined "inline" with one or more attributes separated by vertical pipes ("|") inside square brackets.

  * All of the attributes are optional. However, the `structure` attribute is required to compute the formula, molecular weight, and charge of the biopolymer.
  * Attributes are separated by vertical pipes ("|").
  * Attributes and their values are separated by colons (":").
  * White spaces are ignored.
  * The values of the `id`, `name`, `synonym`, and `comments` attributes must be enclosed in quotes ('"').
  * The namespace and id of each identifer must be separated by "/".

* The positions of the monomers in the string indicate their location in the sequence, as illustrated by the following DNA polymers:

  * ``ACTGCC``: represents alphabet defined thymine at the third position
  * ``aCTGCC``: represents alphabet defined 6-methyladenine at the first position
  * ``ACGC{dI}``: represents alphabet defined hypoxanthine at the last position
  * ``AC[id: "c1" | name: "custom_1"]GC``: represents inline defined custom_1 at the third position

Sections 2.1.1 - 2.1.3 describe the attributes of monomers.


Structure
^^^^^^^^^

The ``structure`` attribute describes the chemical structure of the monomer as a SMILES-encoded string::

    [id: "dI" |
        structure: "O=C1NC=NC2=C1N=CN2"
        ]

We recommend defining this attribute for each monomer. Theis attribute must be defined to calculate the structure, formula, molecular weight, and charge of the biopolymer.


Linkages
^^^^^^^^

The ``backbone-bond-atom``, ``backbone-displaced-atom``, ``left-bond-atom``, ``left-displaced-atom``, ``right-bond-atom``, and ``right-displaced-atom`` attributes describe the linkages between the monomer and the backbone and between successive monomers::

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


Grammar
^^^^^^^

The following is the definition of the `BpForms` grammar. The grammar is defined in `Lark syntax <https://lark-parser.readthedocs.io/en/latest/grammar/>`_ which is based on `EBNF syntax <https://en.wikipedia.org/wiki/Extended_Backus%E2%80%93Naur_form>`_.

.. literalinclude:: ../bpforms/grammar.lark
    :language: text


Examples
^^^^^^^^

* DNA::

    ACGT[id: "dI" | structure: "O=C1NC=NC2=C1N=CN2" | backbone-bond-atom: Monomer / N / 10 / 0 | backbone-displaced-atom: Monomer / H / 10 / 0]AG{m2A}

* RNA::

    {01G}CGU[id: "01A" | structure: "COC1C(O)C(OC1n1cnc2c1ncn(c2=N)C)CO"]AG[id: "019A" | structure: "COC1C(O)C(OC1n1cnc2c1ncn(c2=O)C)CO"]

* Protein::

    ARGKL[id: "AA0318" | structure: "COC(=O)[C@@H]([NH3+])CCCC[NH3+]"]YRCG[id: "AA0567" | structure: "CC=CC(=O)NCCCC[C@@H](C=O)[NH3+]"]


Comparison to ProForma Proteoform Notation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The BpForms syntax was inspired by the `ProForma Proteoform Notation <http://www.topdownproteomics.org/resources/proforma/>`_. BpForms improves upon this syntax in several ways:

* BpForms separates the representation of non-canonical biopolymers from the chemical processes which generate them.
* BpForms can represent any modification and, therefore, is not limited to modifications that have been previously enumerated in databases and ontologies. This is necessary to represent the combinatorial complexity of non-canonical DNA, RNA, and proteins.
* BpForms can capture additional uncertainty in the structures of biopolymers: uncertainty in the position of a non-canonical monomeric form within a sequence, and uncertainty in the chemical identity of a non-canonical monomeric form (e.g., deviation from its expected mass or charge).
* BpForms has a concrete grammar. This enables error checking, as well as the calculation of chemical formulae, masses, and charges, which is essential for modeling and other applications.
* We have written software tools for verifying descriptions of non-canonical biopolymers and calculating their properties
