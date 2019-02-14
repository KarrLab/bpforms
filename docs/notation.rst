.. _notation:

`BpForms` notation
------------------

The `BpForms` notation unambiguously represents the primary structure of biopolymer forms that contain canonical and non-canonical monomers using (a) a syntax similar to FASTA and (b) extended alphabets for DNA, RNA, and proteins based on DNAmod, MODOMICS, and RESID. This enables `BpForms` to calculate the formula, molecular, weight, and charge of biopolymer forms.

* Canonical monomers are indicated by their single upper case character codes.
* Non-canonical monomers defined in the alphabets that have single-character codes are indicated by these codes.
* Non-canonical monomers defined in the alphabets that have multiple-character codes are indicated by these codes delimited by curly brackets.
* Non-canonical monomers that are not defined in the alphabet can be defined "inline" with multiple attributes separated by vertical pipes ("|") enclosed inside square brackets. The structures of these monomers are defined in InChI format using the ``structure`` attribute. Additional attributes can provide metadata about monomers such as their ids and names.

`BpForms` contains three pre-built canonical DNA, RNA and protein alphabets and three extended DNA, RNA, and protein alphabets based on DNAmod, MODOMICS, and RESID. Users can also create additional custom alphabets.

Examples:

* ``[id: "dI" | name: "deoxyinosine"]ACGC``: represents deoxyinosine at the first position
* ``AC[id: "dI" | name: "deoxyinosine"]GC``: represents deoxyinosine at the third position
* ``ACGC{m6A}``: represents methyl-6-adenosine at the last position


Structure
^^^^^^^^^

The ``structure`` attribute describes the chemical structure of the monomer as an InChI-encoded string. We recommend that each monomer have a structure. This attribute must be defined to calculate the formula, molecular weight, and charge of the biopolymer.::

    [id: "dI" |
        structure: InChI=1S
            /C10H12N4O4
            /c15-2-6-5(16)1-7(18-6)14-4-13-8-9(14)11-3-12-10(8)17
            /h3-7,15-16H,1-2H2,(H,11,12,17)
            /t5-,6+,7+
            /m0
            /s1
        ]


Uncertainty
^^^^^^^^^^^

`BpForms` can also represent two types of uncertainty in the structures of biopolymer forms.

* The ``delta-mass`` and ``delta-charge`` attributes can describe uncertainty in the chemical identities of monomers. For example, ``[id: "dAMP" | delta-mass: 1 | delta-charge: 1]`` indicates the presence of an additional proton exact location is not known.
* The ``position`` attribute can describe uncertainty in the position of a monomer. For example, ``[id: "5mC" | position: 2-3]`` indicates that 5mC may occur anywhere between the second and third position.


Metadata
^^^^^^^^

`BpForms` can also represent several types of metadata:

* The ``id`` and ``name`` attributes provide human-readable labels for monomers. Only one id and one name is allowed per monomer.::

    [id: "dI"
        | name: "deoxyinosine"
        ]

* The ``synonym`` attribute provides additional human-readable labels. Each monomer can have multiple synonyms.::

    [id: "dI"
        | synonym: "2'-deoxyinosine"
        | synonym: "2'-deoxyinosine, 9-[(2R,4S,5R)-4-hydroxy-5-(hydroxymethyl)tetrahydrofuran-2-yl]-9H-purin-6-ol"
        ]

* The ``identifier`` attribute describes references to entries in external databases. Each monomer can have multiple identifiers. The namespaces and ids of identifers must be separated by "/".::

    [id: "dI"
        | identifier: "biocyc.compound" / "DEOXYINOSINE"
        | identifier: "chebi" / "CHEBI:28997"
        | identifier: "pubchem.compound" / "65058"
        ]

* The ``base-monomer`` attribute describes other monomer(s) which the monomer is generated from. The value of this attribute must be the code of a monomer in the 
alphabet. Each monomer can have one or more bases. This annotation is needed to generate more informative FASTA sequences for ``BpForm``s::

    [id: "m2A"
        | base-monomer: "A"
        ]

* The ``comments`` attribute describes additional information about each monomer. Each monomer can only have one comment.::

    [id: "dI"
        | comments: "A purine 2'-deoxyribonucleoside that is inosine in which the
                     hydroxy group at position 2' is replaced by a hydrogen."
        ]


Syntax
^^^^^^

* Monomers that are in the alphabet are indicated by a single character or multiple characters delimiated by curly brackets.
* Monomers that are not in the alphabet are defined "inline" with one or more attributes separated by verticle pipes ("|") inside square brackets.

  * All of the attributes can optional. However, the `structure` attribute is required to compute the formula, molecular weight, and charge of the biopolymer..
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
    
    ACGT[id: "dI" | structure: InChI=1S/C10H12N4O4/c15-2-6-5(16)1-7(18-6)14-4-13-8-9(14)11-3-12-10(8)17
    /h3-7,15-16H,1-2H2,(H,11,12,17)/t5-,6+,7+/m0/s1]AG[id: "m6A" | structure: InChI=1S/C6H7N5
    /c1-7-5-4-6(10-2-8-4)11-3-9-5/h2-3H,1H3,(H2,7,8,9,10,11)]{m2A}

* RNA:: 

    {m6A}CGU[id: "m1G" | structure: InChI=1S/C11H15N5O5/c1-15-9(20)5-8(14-11(15)12)16(3-13-5)10-7(19)6(18)4(2-17)21-10
    /h3-4,6-7,10,17-19H,2H2,1H3,(H2,12,14)/t4-,6-,7-,10-/m1/s1]AG[id: "m1A" | structure: InChI=1S/C11H15N5O4
    /c1-15-3-14-10-6(9(15)12)13-4-16(10)11-8(19)7(18)5(2-17)20-11/h3-5,7-8,11-12,17-19H,2H2,1H3/t5-,7-,8-,11-/m1/s1]

* Protein::

    ARGKL[id: "m3Arg" | structure: InChI=1S/C7H16N4O2/c1-4(5(8)6(12)13)2-3-11-7(9)10
    /h4-5H,2-3,8H2,1H3,(H,12,13)(H4,9,10,11)/t4?,5-/m0/s1]YRCG[id: "lysidine" | structure: InChI=1S/C4H8N2
    /c1-4-5-2-3-6-4/h2-3H2,1H3,(H,5,6)]
