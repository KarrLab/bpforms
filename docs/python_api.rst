.. _python_api:

Python API
----------

The following tutorial illustrates how to use the `BpForms` Python API. An `interactive version of this tutorial <https://sandbox.karrlab.org/notebooks/bpforms/Tutorial.ipynb>`_ is also available in the whole-cell modeling sandbox.

Importing `BpForms`
^^^^^^^^^^^^^^^^^^^

Run this command to import `BpForms`::

    import bpforms


Creating biopolymer forms
^^^^^^^^^^^^^^^^^^^^^^^^^

Use the `BpForms` grammar and the ``bpforms.BpForm.from_str`` method to create an instance of ``bpforms.BpForm`` that represents a form of a biopolymer::

    dna_form = bpforms.DnaForm().from_str('ACG{m2C}AC')


Getting and setting monomeric forms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Individual monomeric forms and slices of monomeric forms can be get and set similar to lists::

    dna_form[0]
        => <bpforms.core.Monomer at 0x7fb365341240>

    dna_form[1] = bpforms.dna_alphabet.monomers.A

    dna_form[1:3]
        => [<bpforms.core.Monomer at 0x7fb365341240>, <bpforms.core.Monomer at 0x7fb365330cf8>]

    dna_form[1:3] = bpforms.DnaForm().from_str('TA')


Getting and setting the base of a monomeric form
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Optionally, `BpForms` can track the monomeric forms that are generated from a monomeric form (e.g. m2A is generated from A). This can be get and set using the ``bpforms.Monomer.base_monomers`` attribute. This attribute is a ``set`` of ``bpforms.Monomer``::

    di_monomer = dna_form[3]
    di_monomer.base_monomers
        => set(<bpforms.core.Monomer at 0x7fb365341240>)


Protonation and tautomerization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculate the major protonation and tautomerization state of each monomeric form in the biopolymer form::

    dna_form.get_major_micro_species(8., major_tautomer=True)


Calculation of physical properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use these commands to calculate the length, formula, molecular weight, and charge of the biopolymer form::

    len(dna_form)
        => 6

    dna_form.get_formula()
        => AttrDefault(<class 'float'>, False, {'C': 59.0, 'N': 23.0, 'O': 35.0, 'P': 6.0, 'H': 72.0})

    dna_form.get_mol_wt()
        => 1849.193571988

    dna_form.get_charge()
        => -7


Generating IUPAC/IUBMB sequences for `BpForms`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``get_canonical_seq`` method generates IUPAC/IUBMB representations of `BpForms`. Where annotated, this method uses the ``base_monomers`` attribute to represent modified monomeric forms using the code for their root (e.g. m2A is represented as "A"). Monomeric forms that don't have their base annotated are represented as "N" and "X" for nucleic acids and proteins, respectively::

    dna_form.get_canonical_seq()
        => ATANAC


Determine if two biopolymers describe the same structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use the following command to determine if two instances of :obj:`BpForm` describe the same biopolymer::

    dna_form_1 = bpforms.DnaForm().from_str('ACGT')
    dna_form_2 = bpforms.DnaForm().from_str('ACGT')
    dna_form_3 = bpforms.DnaForm().from_str('GCTC')

    dna_form_1.is_equal(dna_form_2)
        => True

    dna_form_1.is_equal(dna_form_3)
        => False
