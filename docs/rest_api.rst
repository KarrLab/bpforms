.. _rest_api:

REST API
--------

The REST JSON API enables researchers to programmatically validate and calculate the properties of biopolymer forms.

The base URL for the REST API is `https://bpforms.org/api/ <https://bpforms.org/api/>`_. The REST API provides two endpoints.

* ``/version``: Returns the version of the API.::

    0.0.1

* ``/get-properties?base-seq={}``: If the form is valid, this returns the length, formula, mass, and charge of the form. If the form is invalid, this returns an error message.::

    {
        length: <int>,
        formula: <associative array>,
        mass: <float>,
        charge: <int>,
    }

    OR 

    {
        error: <str>,
    }

* ``/get-properties?base-seq={}&ph={}``: If the form is valid, this protonates the form to the specified pH and returns the length, formula, mass, and charge of the form. If the form is invalid, this returns an error message.::

    {
        base-seq: <str>,
        length: <int>,
        formula: <associative array>,
        mass: <float>,
        charge: <int>,
    }

    OR 

    {
        error: <str>,
    }
