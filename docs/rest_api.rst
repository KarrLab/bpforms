.. _rest_api:

REST API
--------

The REST JSON API enables researchers to programmatically validate and calculate the properties of biopolymer forms.

The base URL for the REST API is `https://bpforms.org/api/ <https://bpforms.org/api/>`_. The REST API provides two endpoints.

* ``/``: Returns metadata about the API including a description, the version, and a list of endpoints.::

    {
        description: <string>,
        endpoints: <array>,
        version: <string>
    }

* ``/alphabet``: Returns a list of available alphabets.::

    ['dna', 'rna', 'protein']

* ``/alphabet/{alphabet: string}``: Returns an associative array of bases in the alphabet and their properties.::

    {
        A: {
            id: 'dAMP',
            structure: 'InChI=1s/...',
            ...
        },
        ...
    }

* ``/bpform/properties/{alphabet: rna, dna, or  protein}/{base_seq: string}(/{ph: float})?``: If the form is valid, this, optionally, protonates the form to the specified pH and returns the length, formula, mass, and charge of the form. If the form is invalid, this returns an error message.::

    {
        alphabet: <string>,
        base_seq: <string>,
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
