.. _python_api:

Python API
----------

The following tutorial illustrates how to use the `BpForms` Python API. An `interactive version of this tutorial <https://sandbox.karrlab.org/notebooks/bpforms/Tutorial.ipynb>`_ is also available in the whole-cell modeling sandox.

Importing `BpForms`
^^^^^^^^^^^^^^^^^^^

Run this command to import `BpForms`.::

    import bpforms


Creating biopolymer forms
^^^^^^^^^^^^^^^^^^^^^^^^^

Use the `BpForms` notation and the ``bpforms.BpForm.from_str`` method to create an instance of ``bpforms.BpForm`` that represents a form of a biopolymer.::

    dna_form = bpforms.DnaForm().from_str('''ACG[
        id: "dI" 
        | structure: "O=C1NC=NC2=C1N=CN2"
        | base-monomer: "A"
        ]AC'''.replace('\n', '').replace(' ', ''))


Getting and setting monomers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Individual monomers and slices of monomers can be get and set similar to lists.::

    dna_form[0]
        => <bpforms.core.Monomer at 0x7fb365341240>
    
    dna_form[1] = bpforms.dna_alphabet.monomers.A
    
    dna_form[1:3] 
        => [<bpforms.core.Monomer at 0x7fb365341240>, <bpforms.core.Monomer at 0x7fb365330cf8>]
    
    dna_form[1:3] = bpforms.DnaForm().from_str('TA')


Getting and setting the base of a monomer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Optionally, `BpForms` can track the monomers that are generated from a monomer (e.g. m2A is generated from A). This can be get and set using the ``bpforms.Monomer.base_monomers`` attribute. This attribute is a ``set`` of ``bpforms.Monomer``.::

    di_monomer = dna_form[3]
    di_monomer.base_monomers
        => set(<bpforms.core.Monomer at 0x7fb365341240>)
    di_monomer.base_monomers.add(bpforms.Monomer())


Protonation and tautomerization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculate the major protation and tautomerization state of each monomer in the biopolymer form.::

    dna_form.get_major_micro_species(8., major_tautomer=True)


Calculation of physical properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use these commands to calculate the length, formula, molecular weight, and charge of the biopolymer form.::

    len(dna_form)
        => 6
    
    dna_form.get_formula()
        => AttrDefault(<class 'float'>, False, {'C': 59.0, 'N': 24.0, 'O': 37.0, 'P': 5.0, 'H': 66.0})
    
    dna_form.get_mol_wt()
        => 1858.17680999
    
    dna_form.get_charge()
        => -7


Generating FASTA sequences for `BpForms`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``get_fasta`` method generates FASTA representations of `BpForms`. Where annotated, this method uses the ``base_monomers`` attribute to represent modified monomers using the code for their root (e.g. m2A is represented as "A"). Monomers that don't have their base annotated are represented as "N".::

    dna_form.get_fasta()
        => ACGAAC


Determine if two biopolymers describe the same structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use the following command to determine if two instances of :obj:`BpForm` describe the same biopolymer.::

    dna_form_1 = bpforms.DnaForm().from_str('ACGT')
    dna_form_2 = bpforms.DnaForm().from_str('ACGT')
    dna_form_3 = bpforms.DnaForm().from_str('GCTC')

    dna_form_1.is_equal(dna_form_2)
        => True
    
    dna_form_1.is_equal(dna_form_3)
        => False
