""" Test of bpforms.core

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms import core
from bpforms import rest
from bpforms.alphabet import dna
from wc_utils.util.chem import EmpiricalFormula
import bpforms
import mock
import unittest


class RestTestCase(unittest.TestCase):
    def tearDown(self):
        dna.dna_alphabet.from_yaml(dna.filename)

    def test_PrefixMiddleware(self):
        rest.PrefixMiddleware(rest.app).__call__({'PATH_INFO': 'x'}, lambda x, y: None)

    def test_get_bpform_properties(self):
        client = rest.app.test_client()
        rv = client.post('/api/bpform/', json=dict(alphabet='dna', monomer_seq='ACGT'))
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            'alphabet': 'dna',
            'monomer_seq': 'ACGT',
            'length': 4,
            'formula': dict(EmpiricalFormula('C39H46O25N15P4')),
            'mol_wt': 1248.772047992,
            'charge': -5,
        })

    def test_get_bpform_properties_with_ph(self):
        client = rest.app.test_client()
        rv = client.post('/api/bpform/', json=dict(
            alphabet='dna',
            monomer_seq='ACGT',
            ph=7.,
            major_tautomer=True))
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json()['alphabet'], 'dna')
        self.assertEqual(rv.get_json()['monomer_seq'], 'ACGT')
        self.assertEqual(rv.get_json()['length'], 4)

    def test_get_bpform_properties_errors(self):
        client = rest.app.test_client()

        rv = client.post('/api/bpform/', json=dict(alphabet='lipid', monomer_seq='ACGT'))
        self.assertEqual(rv.status_code, 400)

        rv = client.post('/api/bpform/', json=dict(alphabet='dna', monomer_seq='ACG[T'))
        self.assertEqual(rv.status_code, 400)

        rv = client.post('/api/bpform/', json=dict(alphabet='dna', monomer_seq='ACGT', ph='a'))
        self.assertEqual(rv.status_code, 400)

        rv = client.post('/api/bpform/', json=dict(alphabet='dna', monomer_seq='ACGT', ph=7., major_tautomer='a'))
        self.assertEqual(rv.status_code, 400)

    def test_get_alphabets(self):
        client = rest.app.test_client()
        rv = client.get('/api/alphabet/')
        self.assertEqual(rv.status_code, 200)
        alphabets = rv.get_json()
        self.assertEqual(alphabets['dna']['name'], 'DNAmod DNA nucleobases')
        self.assertEqual(alphabets['canonical_dna']['name'], 'Canonical DNA nucleobases')

    def test_get_alphabet(self):
        client = rest.app.test_client()

        rv = client.get('/api/alphabet/dna/')
        self.assertEqual(rv.status_code, 200)
        alphabet = core.Alphabet().from_dict(rv.get_json())
        self.assertTrue(dna.dna_alphabet.is_equal(alphabet))

        rv = client.get('/api/alphabet/lipid/')
        self.assertEqual(rv.status_code, 400)
