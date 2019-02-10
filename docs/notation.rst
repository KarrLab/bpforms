.. _notation:

`BpForms` notation
------------------

The `BpForms` notation represents biopolymer forms like FASTA sequences, using an extended alphabet to describe canonical and non-canonical nucleotides/amino acids. Each nucleotide/amino acid is unambiguously defined by its ``structure`` attribute in InChI format. Optional attributes (see below) may be added, using “|” as attribute separator. `BpForms` contains an in built-in alphabet consisting of all canonical and hundreds of non-canonical nucleotides/amino acids, and supports user defined custom alphabets. In the FASTA like biopolymer sequence the nucleotides/amino acids are represented using either of three possible notations. (1) Any arbitrary nucleotide/amino acid (even if not defined in the used alphabet) can be defined from whithin the sequence by a ``structure`` attribute, enclosed in square brackets (inline_base notation). (2) Non-canonical nucleotides/amino acids of the alphabet can be denoted by their ``id`` attribute, enclosed in parenthesis. (3) Canonical nucleoteds/amino acids can be denoted using their single letter code.

    * ``[id: "dI" | name: "deoxyinosine"]ACGC``: represents deoxyinosine at the first position
    * ``AC[id: "dI" | name: "deoxyinosine"]GC``: represents deoxyinosine at the third position
    * ``ACGC(m6A)``: represents methyl-6-adenosine at the last position


Structure
^^^^^^^^^

The ``structure`` attribute describes the chemical structure of modified residues as an InChI-encoded string. Each modified residue must have exactly one structure.::

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

* The ``delta-mass`` and ``delta-charge`` attributes can describe uncertainty in the chemical identities of modified residues. For example, ``[id: "dAMP" | delta-mass: 1 | delta-charge: 1]`` indicates the presence of an additional proton.
* The ``position`` attribute can describe uncertainty in the position of a modified nucleotide/amino acid within the polymer sequence. For example, ``[id: "dI" | position: 2-3]`` indicates that deoxyinosine may occur anywhere between the second and third position.


Metadata
^^^^^^^^

`BpForms` can also represent several types of metadata:

* The ``id`` and ``name`` attributes can be used to provide human-readable labels for nucleotides/amino acids. Only one id and one name is allowed per residue::

    [id: "dI"
        | name: "deoxyinosine"
        ]

* The ``synonym`` attribute can be used to provide additional human-readable labels. Each residue can have multiple synonyms::

    [id: "dI"
        | synonym: "2'-deoxyinosine"
        | synonym: "2'-deoxyinosine, 9-[(2R,4S,5R)-4-hydroxy-5-(hydroxymethyl)tetrahydrofuran-2-yl]-9H-purin-6-ol"
        ]

* The ``identifier`` attribute can be used to provide references to entries in external databases. Each residue can have multiple identifiers. The namespaces and ids of identifers must be separated by "/"::

    [id: "dI"
        | identifier: biocyc.compound/DEOXYINOSINE
        | identifier: chebi/CHEBI:28997
        | identifier: pubchem.compound/65058
        ]

* The ``comments`` attribute can be used to describe additional information about each residue. Each residue can only have one comment::

    [id: "dI"
        | comments: "A purine 2'-deoxyribonucleoside that is inosine in which the
                     hydroxy group at position 2' is replaced by a hydrogen."
        ]


Syntax
^^^^^^

* The position of a residue describes their location within the sequence.
* Canonical residues can be denoted with their single letter code.
* Residues in the alphabet that can be denoted by parentheses ("(" and ")").
* Any residue can be denoted by square brackets ("[" and "]").
  * Attributes of modified residues are separated by vertical pipes ("|").
  * Attributes and their values are separated by colons (":").
  * White spaces are ignored around the attribute and field separators.
  * All of the attributes can optional. However, the structure attribute is necessary to compute the formula, molecular weight, and charge.
  * String-valued attribute values (id, name, synonym, comments) must be enclosed in quotes ('"').
  * The namespaces and ids of identifers must be separated by "/".


Grammar
^^^^^^^

The following is the definition of the `BpForms` grammar.::

    ?start: seq
    seq: base+
    ?base: alphabet_base | inline_base
    alphabet_base: CHAR | DELIMITED_CHARS
    inline_base: "[" WS* inline_base_attr (ATTR_SEP inline_base_attr)* WS* "]"
    ?inline_base_attr: id | name | synonym | identifier | structure | delta_mass | delta_charge | position | comments
    ?id: "id" FIELD_SEP ESCAPED_STRING
    ?name: "name" FIELD_SEP ESCAPED_STRING
    ?synonym: "synonym" FIELD_SEP ESCAPED_STRING
    ?identifier: "identifier" FIELD_SEP identifier_ns IDENTIFIER_SEP identifier_id
    ?identifier_ns: ESCAPED_STRING
    ?identifier_id: ESCAPED_STRING
    ?structure: "structure" FIELD_SEP INCHI
    ?delta_mass: "delta-mass" FIELD_SEP DALTON
    ?delta_charge: "delta-charge" FIELD_SEP CHARGE
    ?position: "position" FIELD_SEP START_POSITION? "-" END_POSITION?
    ?comments: "comments" FIELD_SEP ESCAPED_STRING
    ATTR_SEP: WS* "|" WS*
    FIELD_SEP: WS* ":" WS*
    IDENTIFIER_SEP: WS* "/" WS*
    CHAR: /[A-Z]/
    DELIMITED_CHARS: "(" /[^\(\) ]*[A-Z][^\(\) ]*/ ")"
    INCHI: /InChI=1S\/[A-Za-z0-9\(\)\-\+,\/]+/
    DALTON: /[\-\+]?[0-9]+(\.[0-9]*)?/
    CHARGE: /[\-\+]?[0-9]+/
    START_POSITION: INT
    END_POSITION: INT
    WS: /[ \t\f\r\n]/+
    ESCAPED_STRING : "\"" _STRING_ESC_INNER "\""
    _STRING_ESC_INNER: _STRING_INNER /(?<!\\)(\\\\)*?/
    _STRING_INNER: /.*?/
    INT: DIGIT+
    DIGIT: "0".."9"


Examples
^^^^^^^^

* DNA::
    
    ACGT[id: "dI" | structure: InChI=1S/C10H12N4O4/c15-2-6-5(16)1-7(18-6)14-4-13-8-9(14)11-3-12-10(8)17
    /h3-7,15-16H,1-2H2,(H,11,12,17)/t5-,6+,7+/m0/s1]AG[id: "m6A" | structure: InChI=1S/C6H7N5
    /c1-7-5-4-6(10-2-8-4)11-3-9-5/h2-3H,1H3,(H2,7,8,9,10,11)](m2A)

* RNA:: 

    (m6A)CGU[id: "m1G" | structure: InChI=1S/C11H15N5O5/c1-15-9(20)5-8(14-11(15)12)16(3-13-5)10-7(19)6(18)4(2-17)21-10
    /h3-4,6-7,10,17-19H,2H2,1H3,(H2,12,14)/t4-,6-,7-,10-/m1/s1]AG[id: "m1A" | structure: InChI=1S/C11H15N5O4
    /c1-15-3-14-10-6(9(15)12)13-4-16(10)11-8(19)7(18)5(2-17)20-11/h3-5,7-8,11-12,17-19H,2H2,1H3/t5-,7-,8-,11-/m1/s1]

* Protein::

    ARGKL[id: "m3Arg" | structure: InChI=1S/C7H16N4O2/c1-4(5(8)6(12)13)2-3-11-7(9)10
    /h4-5H,2-3,8H2,1H3,(H,12,13)(H4,9,10,11)/t4?,5-/m0/s1]YRCG[id: "lysidine" | structure: InChI=1S/C4H8N2
    /c1-4-5-2-3-6-4/h2-3H2,1H3,(H,5,6)]
