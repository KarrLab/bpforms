.. _python_api:

Python API
----------


Importing `BpForms`
^^^^^^^^^^^^^^^^^^^

First, run this command to import ``BpForms``.::

    import bpforms


Creating biopolymer forms
^^^^^^^^^^^^^^^^^^^^^^^^^

Second, use the `BpForms` notation and the `BpForm.from_str` method to create an instance of `BpForm` that represents a form of a biopolymer.::

    dna_form = bpforms.DnaForm().from_str('''ACG[
        id: "dI" | 
        structure: InChI=1S
            /C10H12N4O4
            /c15-2-6-5(16)1-7(18-6)14-4-13-8-9(14)11-3-12-10(8)17
            /h3-7,15-16H,1-2H2,(H,11,12,17)
            /t5-,6+,7+
            /m0
            /s1
        ]AC'''.replace('\n', '').replace(' ', ''))


Getting a setting residues
^^^^^^^^^^^^^^^^^^^^^^^^^^
Third, access individual residues and slices of residues similar to lists.::

    dna_form[0]
        => <bpforms.core.Base at 0x7fb365341240>
    
    dna_form[1] = bpforms.dna_alphabet.bases.A
    
    dna_form[1:3] 
        => [<bpforms.core.Base at 0x7fb365341240>, <bpforms.core.Base at 0x7fb365330cf8>]
    
    dna_form[1:3] = bpforms.DnaForm().from_str('TA')


Protonation
^^^^^^^^^^^
Fourth, calculate the major protation state of each residue in the biopolymer form.::

    dna_form.protonate(8.)


Calculation of physical properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Fifth, use these commands to calculate the length, formula, molecular weight, and charge of the biopolymer form.::

    len(dna_form)
        => 6
    
    dna_form.get_formula()
        => AttrDefault(<class 'float'>, False, {'C': 59.0, 'N': 24.0, 'O': 37.0, 'P': 5.0, 'H': 66.0})
    
    dna_form.get_mol_wt()
        => 1858.17680999
    
    dna_form.get_charge()
        => -7


Determine if two biopolymers describe the same structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Six, use the following command to determine if two instances of `BpForm` describe the same biopolymer.::

    dna_form_1 = bpforms.DnaForm().from_str('ACGT')
    dna_form_2 = bpforms.DnaForm().from_str('ACGT')
    dna_form_3 = bpforms.DnaForm().from_str('GCTC')

    dna_form_1.is_equal(dna_form_2)
        => True
    
    dna_form_1.is_equal(dna_form_3)
        => False
