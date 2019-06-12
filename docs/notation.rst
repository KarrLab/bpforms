.. _notation:

`BpForms` notation
------------------

The `BpForms` notation unambiguously represents the primary structure of biopolymer forms that contain canonical and non-canonical monomeric forms using (a) a syntax similar to FASTA and (b) extended alphabets for DNA, RNA, and proteins to describe monomeric forms. `BpForms` contains three pre-built canonical DNA, RNA and protein alphabets and three extended DNA, RNA, and protein alphabets curated from the DNAmod, MODOMICS, and RESID databases, respectively. Users can also create additional custom alphabets. These alphabets are associated with their corresponding biopolyer form, which allows `BpForms` to calculate the chemical structure (e.g. in SMILES format), chemical formula, molecular weight, and charge of a biopolymer.

* Monomeric forms that are present in the alphabet are indicated by a single character or multiple characters delimited by curly brackets.
* Monomeric forms that are not in the alphabet are defined "inline" with one or more attributes separated by vertical pipes ("|") inside square brackets.

  * All of the attributes are optional. However, the `structure` attribute is required to compute the formula, molecular weight, and charge of the biopolymer.
  * Attributes are separated by vertical pipes ("|").
  * Attributes and their values are separated by colons (":").
  * White spaces are ignored.
  * The values of the `id`, `name`, `synonym`, and `comments` attributes must be enclosed in quotes ('"').
  * The namespace and id of each identifer must be separated by "/".

* The positions of the monomeric forms in the string indicate their location in the sequence, as illustrated by the following DNA polymers:

  * ``ACTGCC``: represents alphabet defined thymine at the third position
  * ``aCTGCC``: represents alphabet defined 6-methyladenine at the first position
  * ``ACGC{dI}``: represents alphabet defined hypoxanthine at the last position
  * ``AC[id: "c1" | name: "custom_1"]GC``: represents inline defined custom_1 at the third position

Sections 2.1.1 - 2.1.3 describe the attributes of monomeric forms.


Structure
^^^^^^^^^

The ``structure`` attribute describes the chemical structure of the monomeric form as a SMILES-encoded string::

    [id: "dI" |
        structure: "O=C1NC=NC2=C1N=CN2"
        ]

We recommend defining this attribute for each monomeric form. Theis attribute must be defined to calculate the structure, formula, molecular weight, and charge of the biopolymer.


Linkages
^^^^^^^^

The ``backbone-bond-atom``, ``backbone-displaced-atom``, ``right-bond-atom``, ``right-displaced-atom``, ``left-bond-atom``, and ``left-displaced-atom`` attributes describe the linkages between the monomeric form and the backbone and between successive monomeric forms. The values of these attributes represent (i) the element in the monomeric form involved in the bond, (ii) the atom index of the element in the monomeric form involved in the bond, and (iii) the charge of atom involved in the bond.::

    [id: "dI"
        | structure: "O=C1NC=NC2=C1N=CN2"
        | backbone-bond-atom: N10
        | backbone-displaced-atom: H10
        ]

We recommend defining these attributes for each monomeric form. These attributes must be defined to calculate the structure, formula, molecular weight, and charge of the biopolymer.


Circular topology
^^^^^^^^^^^^^^^^^

The ``circular`` attribute indicates that a biopolymer has a cicular topology::

    AAA | circular


Intrachain crosslinks
^^^^^^^^^^^^^^^^^^^^^

The ``crosslink``,  ``right-bond-atom``, ``right-displaced-atom``, ``left-bond-atom``, and ``left-displaced-atom`` attributes describe additional covalent bonds between non-adjacent monomeric forms, such as DNA crosslinks caused chemotherapeutics and disulfide bonds between cysteines in proteins. The values of the ``right-bond-atom``, ``right-displaced-atom``, ``left-bond-atom``, and ``left-displaced-atom`` attributes indicate the atoms involved in or displaced by the formation of each covalent bond. The values of these attributes represent (i) the index of the monomeric form involved in the bond, (ii) the element in the monomeric form involved in the bond, (iii) the atom index of the element in the monomeric form involved in the bond, and (iv) the charge of atom involved in the bond.::
  
  AC | crosslink: [
    right-bond-atom: 1N4 |
    left-bond-atom: 2C8 |
    right-displaced-atom: 1H4+2 |
    left-displaced-atom: 2H8+1
    ]

  AC | crosslink: [
    left-bond-atom: 1C2 |
    right-bond-atom: 2N1 |
    left-displaced-atom: 1N1 |
    left-displaced-atom: 1H1 |
    left-displaced-atom: 1H1 |
    right-displaced-atom: 2H1
    ]


Uncertainty
^^^^^^^^^^^

`BpForms` can also represent two types of uncertainty in the structures of biopolymer forms.

* The ``delta-mass`` and ``delta-charge`` attributes can describe uncertainty in the chemical identities of monomeric forms. For example, ``[id: "dAMP" | delta-mass: 1 | delta-charge: 1]`` indicates the presence of an additional proton exact location is not known.
* The ``position`` attribute can describe uncertainty in the position of a monomeric form. For example, ``[id: "5mC" | position: 2-3]`` indicates that 5mC may occur anywhere between the second and third position.


Metadata
^^^^^^^^

`BpForms` can also represent several types of metadata:

* The ``id`` and ``name`` attributes provide human-readable labels for monomeric forms. Only one id and one name is allowed per monomeric form::

    [id: "dI"
        | name: "deoxyinosine"
        ]

* The ``synonym`` attribute provides additional human-readable labels. Each monomeric form can have multiple synonyms::

    [id: "dI"
        | synonym: "2'-deoxyinosine"
        | synonym: "2'-deoxyinosine, 9-[(2R,4S,5R)-4-hydroxy-5-(hydroxymethyl)tetrahydrofuran-2-yl]-9H-purin-6-ol"
        ]

* The ``identifier`` attribute describes references to entries in external databases. Each monomeric form can have multiple identifiers. The namespaces and ids of identifers must be separated by "/"::

    [id: "dI"
        | identifier: "DEOXYINOSINE" @ "biocyc.compound"
        | identifier: "CHEBI:28997" @ "chebi"
        | identifier: "65058" @ "pubchem.compound"
        ]

* The ``base-monomer`` attribute describes other monomer form(s) which the monomeric form is generated from. The value of this attribute must be the code of a monomeric form in the alphabet. Each monomeric form can have one or more bases. This annotation is needed to generate more informative FASTA sequences for `BpForms`::

    [id: "m2A"
        | base-monomer: "A"
        ]

* The ``comments`` attribute describes additional information about each monomeric form. Each monomeric form can only have one comment::

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

    ACGT[id: "dI" | structure: "O=C1NC=NC2=C1N=CN2" | backbone-bond-atom: N10 | backbone-displaced-atom: H10]AG{m2A}

* RNA::

    {01G}CGU[id: "01A" | structure: "COC1C(O)C(OC1n1cnc2c1ncn(c2=N)C)CO"]AG[id: "019A" | structure: "COC1C(O)C(OC1n1cnc2c1ncn(c2=O)C)CO"]

* Protein::

    ARGKL[id: "AA0318" | structure: "COC(=O)[C@@H]([NH3+])CCCC[NH3+]"]YRCG[id: "AA0567" | structure: "CC=CC(=O)NCCCC[C@@H](C=O)[NH3+]"]

* Cicular DNA::

    ACGT | circular
