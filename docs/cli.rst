.. _cli:

Command line interface
----------------------

The command line interface provides four functions to easily manipulate `BpForms`-encoded descriptions of biopolymers.

* **Get help with the `BpForms` command line interface.** The following commands return inline help information about the command line interface::

    bpforms
    bpforms -h
    bpforms --help

* **Validate a `BpForms`-encoded description of a form of a biopolymer.** The following command can be used to verify if a sequence of monomeric forms is syntactically and semantically valid. The command line interface will print any errors to the standard error::

    bpforms validate --help

    bpforms validate <alphabet_name> <bpform_sequence>

    bpforms validate dna 'ACGT | circular'
    # Form is valid

    bpforms validate protein 'CRATUG'
    # Form is valid

* **Calculate the protonation of a form of a biopolymer.** The following command can be used to calculate the major protonation and tautomerization state of each monomeric form in a biopolymer form. This command will print the sequence of the major protonation/tautomerization states. Note, this function requires a structure for each monomeric form::

    bpforms get-major-micro-species --help

    bpforms get-major-micro-species <alphabet_name> <bpform_sequence> <ph_value>

    bpforms get-major-micro-species dna 'ACGT | circular' 7
    # Cc1cn(C2CC(O)C(COP(=O)([O-])OC3CC(OC3COP(=O)([O-])OC3CC(OC3COP(=O)([O-])OC3CC(OC3COP(=O)([O-])[O-])n3cnc4c(N)ncnc34)n3ccc(N)nc3=O)n3cnc4c3nc(N)[nH]c4=O)O2)c(=O)[nH]c1=O

    bpforms get-major-micro-species protein 'CRATUG' 7
    # C[C@@H](O)[C@H](NC(=O)[C@H](C)NC(=O)[C@H](CCCNC(=[NH2+])N)NC(=O)[C@@H]([NH3+])CS)C(=O)N[C@@H](C[SeH])C(=O)NCC(=O)[O-]

* **Calculate physical properties of a form of a biopolymer.** The following command can be used to calculate the length, formula, mass, and charge of a bipolymer form. The optional ``ph`` argument can be used to calculate the major protonation and tautomerization state of the biopolymer at a specific pH. Note, this function requires a structure for each monomeric form::

    bpforms get-properties --help

    bpforms get-properties <alphabet_name> <bpform_sequence> [--ph ph]

    bpforms get-properties dna 'ACGT | circular'
    # Length: 4
    # Structure: O(C1CC(OC1COP(=O)([O-])[O-])n1cnc2c1ncnc2N)P(=O)(OCC1C(OP(=O)(OCC2C(OP(=O)(OCC3C(O)CC(O3)n3cc(C)c(=O)[nH]c3=O)[O-])CC(O2)n2cnc3c2nc(N)[nH]c3=O)[O-])CC(O1)n1ccc(nc1=O)N)[O-]
    # Formula: C39H46N15O25P4
    # Molecular weight: 1248.772047992
    # Charge: -5

    bpforms get-properties protein 'CRATUG'
    # Length: 6
    # Structure: C(=O)([C@@H]([NH3+])CS)N[C@H](C(=O)N[C@@H](C)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@H](C(=O)NCC(=O)O)C[SeH])CCCNC(=[NH2+])N
    # Formula: C21H41N9O8SSe
    # Molecular weight: 658.645
    # Charge: 2
