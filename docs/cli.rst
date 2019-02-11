.. _cli:

Command line interface
----------------------

The command line interface provides four functions to easily manipulate BpForms.

* **Help.** The following commands return inline help information about the command line interface::

    bpforms
    bpforms -h
    bpforms --help

* **BpForms validation.** The following command can be used to verify if a base sequence is syntactically valid. The command line interface will print any syntax errors to the standard error::

    bpforms validate <alphabet_name> <bpform_sequence>

* **Protation of biopolymer forms.** The following command can be used to calculate the major protonation state of each base in a biopolymer form. This command will print the sequence of the major protonation states.  Note, this function requires a structure for each modified base.::

    bpforms protonate <alphabet_name> <bpform_sequence> <ph_value>

* **Calculation of the physical properties of biopolymer forms.** The following command can be used to calculate the length, formula, mass, and charge of a bipolymer form. The optional ``ph`` argument can be used to protonate the biopolymer to a specific pH. Note, this function requires a structure for each modified base.::

    bpforms get-properties <alphabet_name> <bpform_sequence> [--ph ph]
