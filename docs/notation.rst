.. _notation:

`BpForms` notation
------------------

The `BpForms` notation represents biopolymer forms as FASTA sequences with (a) bases denoted by multiple characters delimited by curly brackets and (b) modified residues described by multiple attributes separated by "|" inside square brackets. The locations of modified residues are encoded by their positions within the sequence.

    * ``[id: "dI" | name: "deoxyinosine"]ACGC``: represents deoxyinosine at the first position
    * ``AC[id: "dI" | name: "deoxyinosine"]GC``: represents deoxyinosine at the third position
    * ``ACGC{m6A}``: represents methyl-6-adenosine at the last position


Structure
^^^^^^^^^

The ``structure`` attribute can describe the chemical structure of modified residues. This should be an InChI-encoded string. Each modified residue can have one structure.::

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

* The ``delta-mass`` ``delta-charge`` attributes can describe uncertainty in the chemical identities of modified residues. For example, ``[id: "dAMP" | delta-mass: 1 | delta-charge: 1]`` indicates that the presence of an additional hydrogen atom whose exact location is not known.
* The ``position`` attribute can describe uncertainty in the positions of modified residues. For example, ``[id: "dI" | position: 2-3]`` indicates that deoxyinosine may occur anywhere between the second and third position.


Metadata
^^^^^^^^

`BpForms` can also represent several types of metadata:

* The ``id`` and ``name`` attributes can be used to provide human-readable labels for modified residues. Only one id and one name is allowed per residue::

    [id: "dI"
        | name: "deoxyinosine"
        ]

* The ``synonym`` attribute can be used to provide additional human-readable labels. Each residue can have multiple synonyms.::

    [id: "dI"
        | synonym: "2'-deoxyinosine"
        | synonym: "2'-deoxyinosine, 9-[(2R,4S,5R)-4-hydroxy-5-(hydroxymethyl)tetrahydrofuran-2-yl]-9H-purin-6-ol"
        ]

* The ``identifier`` attribute can be used to provide references to entries in external databases. Each residue can have multiple identifiers. The namespaces and ids of identifers must be separated by "/".::

    [id: "dI"
        | identifier: biocyc.compound/DEOXYINOSINE
        | identifier: chebi/CHEBI:28997
        | identifier: pubchem.compound/65058
        ]

* The ``comments`` attribute can be used to describe additional information about each residue. Each residue can only have one comment.::

    [id: "dI"
        | comments: "A purine 2'-deoxyribonucleoside that is inosine in which the
                     hydroxy group at position 2' is replaced by a hydrogen."
        ]


Syntax
^^^^^^

* Bases in the alphabet that are denoted by multiple characters are enclosed in curly brackets ("{" and "}")
* Modified residues are described by square brackets ("[" and "]")
* The position of the square brackets describes their location within the sequence
* Attributes of modified residues are separated by vertical pipes ("|")
* Attributes and their values are separated by colons (":")
* White spaces are ignored around the attribute and field separators
* All of the attributes can optional. However, the structure attribute is necessary to compute the formula, molecular weight, and charge.
* String-valued attribute values (id, name, synonym, comments) must be enclosed in quotes ('"')
* The namespaces and ids of identifers must be separated by "/".


Grammar
^^^^^^^

The following is the definition of the `BpForms` grammar.

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
