.. _rest_api:

REST API
--------

The REST JSON API enables researchers to programmatically validate and calculate the properties of biopolymer forms.

The root URL for the REST API is `https://bpforms.org/api/ <https://bpforms.org/api/>`_. The REST API provides two endpoints.

* ``/alphabet``: Returns a list of available alphabets.::

    [
        {
            'id': 'dna', 
            'name': 'DNA',
            'description': '...'
        },
        ...
    ]

* ``/alphabet/{alphabet: string}`` returns an associative array of monomers in the alphabet and their properties::

    {
        A: {
            id: 'dAMP',
            structure: 'InChI=1s/...',
            ...
        },
        ...
    }

* ``/bpform/{alphabet: rna, dna, or  protein}/{monomer_seq: string}(/{ph: float})?``: optionally, protonates the biopolymer form to the specified pH and returns its length, chemical formula, mass and charge. If the form is invalid, this returns an error message.::

    {
        alphabet: <string>,
        monomer_seq: <string>,
        length: <integer>,
        formula: <associative array>,
        mol_wt: <float>,
        charge: <integer>
    }

    OR 

    {
        message: <string>,
        details: <string>
    }
