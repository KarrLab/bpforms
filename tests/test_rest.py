""" Test of bpforms.core

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms import rest
from wc_utils.util.chem import EmpiricalFormula
import bpforms
import flask_api.exceptions
import mock
import unittest


class RestTestCase(unittest.TestCase):
    def test_ApiException_handle_api_exception(self):
        client = rest.app.test_client()

        rv = client.get('/undefined-endpoint')
        self.assertEqual(rv.status_code, 404)

        @rest.app.route("/error_1/")
        def error_1():
            raise rest.ApiException('the first message', details='the details')

        @rest.app.route("/error_2/")
        def error_2():
            raise flask_api.exceptions.APIException('the second message')

        @rest.app.route("/error_3/")
        def error_3():
            raise Exception('the third message')

        rv = client.get('/error_1/')
        self.assertEqual(rv.status_code, 400)
        self.assertEqual(rv.get_json(), {'message': 'the first message', 'details': 'the details'})

        rv = client.get('/error_2/')
        self.assertEqual(rv.status_code, 500)
        self.assertEqual(rv.get_json(), {'message': 'the second message'})

        rv = client.get('/error_3/')
        self.assertEqual(rv.status_code, 500)
        self.assertEqual(rv.get_json(), None)

    def test_ApiException(self):
        exc = rest.ApiException('the first message', details='the details')
        self.assertEqual(exc.to_dict(), {'message': 'the first message', 'details': 'the details'})
        self.assertEqual(str(exc), 'the first message')

    def test_default(self):
        client = rest.app.test_client()
        rv = client.get('/api/')
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json()['version'], bpforms.__version__)

    def test_get_properties(self):
        client = rest.app.test_client()
        rv = client.get('/api/get-properties/dna/ACGT/')
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            'alphabet': 'dna',
            'base_seq': 'ACGT',
            'length': 4,
            'formula': dict(EmpiricalFormula('C39H46N15O28P4')),
            'mol_wt': 1296.769047992,
            'charge': -5,
        })

    def test_get_properties_with_ph(self):
        client = rest.app.test_client()
        rv = client.get('/api/get-properties/dna/ACGT/7./')
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json()['alphabet'], 'dna')
        self.assertEqual(rv.get_json()['base_seq'], 'ACGT')
        self.assertEqual(rv.get_json()['length'], 4)

    def test_get_properties_errors(self):
        client = rest.app.test_client()

        rv = client.get('/api/get-properties/lipid/ACGT/7./')
        self.assertEqual(rv.status_code, 400)

        rv = client.get('/api/get-properties/dna/ACG[T/7./')
        self.assertEqual(rv.status_code, 400)

        rv = client.get('/api/get-properties/dna/ACGT/a/')
        self.assertEqual(rv.status_code, 400)
