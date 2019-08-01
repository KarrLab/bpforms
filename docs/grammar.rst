.. _grammar:

`BpForms` grammar
------------------

The `BpForms` grammar unambiguously represents the primary structure of biopolymer forms that contain canonical and non-canonical monomeric forms using (a) a syntax similar to IUPAC/IUBMB and (b) extended alphabets for DNA, RNA, and proteins to describe monomeric forms. `BpForms` contains three pre-built canonical DNA, RNA and protein alphabets and three extended DNA, RNA, and protein alphabets curated from DNAmod, MODOMICS, the PDB Chemical Component Dictionary, RESID, and the RNA Modification Database. Users can also create additional custom alphabets. These alphabets are associated with their corresponding biopolyer form, which allows `BpForms` to calculate the chemical structure (e.g., in SMILES format), chemical formula, molecular weight, and charge of a biopolymer.

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

Sections 2.1.1 - 2.1.6 describe the attributes of monomeric forms. Please also see `BpForms.org <https://www.bpforms.org>`_ for more information and examples.


Structure
^^^^^^^^^

The ``structure`` attribute describes the chemical structure of the monomeric form as a SMILES-encoded string, preferably with the atoms canonicaly ordered (e.g., in Open Babel canonical SMILES format).::

    [id: "dI" |
        structure: "OC[C@H]1O[C@H](C[C@@H]1O)[N+]1(C=Nc2c1nc[nH]c2=O)C1CC(C(O1)COP(=O)([O-])[O-])O"
        ]

We recommend defining this attribute for each monomeric form. Theis attribute must be defined to calculate the structure, formula, molecular weight, and charge of the biopolymer.


Bonds with adjacent monomeric forms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``r-bond-atom``, ``r-displaced-atom``, ``l-bond-atom``, and ``l-displaced-atom`` attributes describe the bonds between successive monomeric forms. The values of these attributes represent 

* The element in the monomeric form involved in the bond, 
* The atom index of the element in the monomeric form involved in the bond, according to the atom ordering in the structure of the monomer which, preferably, should be ordered according to the canonical SMILES algorithm, and 
* The charge of atom involved in the bond.

::

    [id: "dI"
        | structure: "OC[C@H]1O[C@H](C[C@@H]1O)[N+]1(C=Nc2c1nc[nH]c2=O)C1CC(C(O1)COP(=O)([O-])[O-])O"
        ]

We recommend defining these attributes for each monomeric form. These attributes must be defined to calculate the structure, formula, molecular weight, and charge of the biopolymer.


Circular topology
^^^^^^^^^^^^^^^^^

The ``circular`` attribute indicates that a biopolymer has a cicular topology::

    AAA | circular


Intrachain crosslinks
^^^^^^^^^^^^^^^^^^^^^

The ``x-link``,  ``r-bond-atom``, ``r-displaced-atom``, ``l-bond-atom``, and ``l-displaced-atom`` attributes describe additional covalent bonds between non-adjacent monomeric forms, such as DNA crosslinks caused chemotherapeutics and disulfide bonds between cysteines in proteins. The values of the ``r-bond-atom``, ``r-displaced-atom``, ``l-bond-atom``, and ``l-displaced-atom`` attributes indicate the atoms involved in or displaced by the formation of each covalent bond. The values of these attributes represent 

  * The index of the monomeric form involved in the bond,
  * The element in the monomeric form involved in the bond,
  * The atom index of the element in the monomeric form involved in the bond, according to the atom ordering in the structure of the monomer which, preferably, should be ordered according to the canonical SMILES algorithm, and 
  * The charge of atom involved in the bond.

::

  AC | x-link: [
    r-bond-atom: 2O1 |
    l-bond-atom: 1P9 |
    r-displaced-atom: 2H1 |
    l-displaced-atom: 1O12-1
    ]

  CRC | x-link: [
    l-bond-atom: 1S11 |
    l-displaced-atom: 1H11 |
    r-bond-atom: 3S11 |
    r-displaced-atom: 3H11
    ]


Uncertainty
^^^^^^^^^^^

`BpForms` can also represent two types of uncertainty in the structures of biopolymer forms.

* The ``delta-mass`` and ``delta-charge`` attributes can describe uncertainty in the chemical identities of monomeric forms. For example, ``[id: "dAMP" | delta-mass: 1 | delta-charge: 1]`` indicates the presence of an additional proton exact location is not known.
* The ``position`` attribute can describe uncertainty in the position of a monomeric form. For example, ``[id: "5mC" | position: 2-3]`` indicates that 5mC may occur anywhere between the second and third position; ``[id: "5mC" | position: 4-8 [A | C]]`` indicates that 5mC may occur at any A or C between the fourth and eighth positions.


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

* The ``base-monomer`` attribute describes other monomer form(s) which the monomeric form is generated from. The value of this attribute must be the code of a monomeric form in the alphabet. Each monomeric form can have one or more bases. This annotation can be used to generate canonical IUPAC/IUBMB sequences for `BpForms`::

    [id: "m2A"
        | base-monomer: "A"
        ]

* The ``comments`` attribute describes additional information about each monomeric form. Each monomeric form can only have one comment::

    [id: "dI"
        | comments: "A purine 2'-deoxyribonucleotide monophosphate that is inosine in which the
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

    ACGT[id: "dI" 
    | structure: "OC[C@H]1O[C@H](C[C@@H]1O)[N+]1(C=Nc2c1nc[nH]c2=O)C1CC(C(O1)COP(=O)([O-])[O-])O"
    | r-bond-atom: O34
    | l-bond-atom: P30
    | r-displaced-atom: H32
    | l-displaced-atom: O33-1
    ]AG{m2A}

* RNA::

    {01G}CGU[id: "01A" 
      | structure: "COC1C(O)C(OC1n1cnc2c1ncn(c2=N)C)COP(=O)([O-])[O-]"
      | r-bond-atom: O5
      | l-bond-atom: P22
      | r-displaced-atom: H5
      | l-displaced-atom: O25-1
    ]
    AG[id: "019A" 
      | structure: "COC1C(O)C(OC1n1cnc2c1ncn(c2=O)C)COP(=O)([O-])[O-]"
      | r-bond-atom: O5
      | l-bond-atom: P22
      | r-displaced-atom: H5
      | l-displaced-atom: O25-1
    ]

* Protein::

    ARGK[id: "AA0567" 
      | structure: "C/C=C/C(=O)NCCCC[C@@H](C(=O)O)[NH3+]"
      | l-bond-atom: N16-1
      | l-displaced-atom: H16
      | l-displaced-atom: H16+1
      | r-bond-atom: C13
      | r-displaced-atom: O15
      | r-displaced-atom: H15
    ]LYRCG[id: "AA0318" 
      | structure: "COC(=O)[C@@H]([NH3+])CCCC[NH3+]"
      | l-bond-atom: N7-1
      | l-displaced-atom: H7
      | l-displaced-atom: H7+1
    ]

* Cicular DNA::

    ACGT | circular
