Limitations, alternatives, and future directions
------------------------------------------------

`BpForms` has several known limitations. Here we describe these limitations and, where possible suggest alternative solutions to these limitations.

* `BpForms` has limited ability to represent non-canonical backbones. However, non-canonical backbones can be encoded within monomeric forms.

    #. Create another alphabet.
    #. Encode the canonical backbone within the monomeric forms of the alphabet (e.g., build an alphabet of nucleotides rather than an alphabet of nucleosides or nucleobases).
    #. Add additional monomeric form to the alphabet that encode non-canonical backbones.

* `BpForms` has limited ability to represent non-canonical backbone-backbone or backbone-monomer bonds. Users can represent non-canonical bonds by encoding non-canonical bonds within instances of ``bpforms.Monomer`` that represent multiple monomeric forms. However, the `BpForms` notation cannot directly describe non-canonical bonds. The `BpForms` notation could be extended to clearly and compactly describe non-canonical bonds by adding an optional syntax for describing the atoms that bond and the atoms that are displaced by bonds. This syntax could use angle brackets to delimit bonds, a vertical pipe to delimit bonded atoms from displaced atoms, and the IUPAC numbering system to indicate the atoms that bonds and the atoms that are displaced (e.g., ``A<O4, P6 | H4, O6>C``).
* ``bpforms.BpForm.circular`` represents the circularity of biopolymers. This attribute must be set programmatically. The `BpForms` notation cannot not represent the circularity biopolymers.
* `BpForms` cannot represent the secondary structure (e.g., base pairing) of biopolymers. Extending `BpForms` to represent secondary structure would enable `BpForms` to represent single strand breaks, crosslinks, and other deviations from the canonical secondary structures of biopolymers.
* `BpForms` cannot represent branched biopolymers such as lipids.

Most importantly, `BpForms` should be extended to represent secondary structure. This would enable `BpForms` to describe DNA forms involved in DNA damage and repair.
