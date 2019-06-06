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
        rv = client.post('/api/bpform/', json=dict(alphabet='dna', seq='ACGT'))
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            'alphabet': 'dna',
            'seq': 'ACGT',
            'length': 4,
            'structure': ('Nc1c2ncn(c2ncn1)C1CC(OP(=O)(OCC2C(OP(=O)(OCC3C(OP(=O)(OCC4C(O)CC(n5cc(C)c(=O)'
                          '[nH]c5=O)O4)[O-])CC(n4c5nc(N)[nH]c(=O)c5nc4)O3)[O-])CC(n3c(=O)nc(N)cc3)O2)[O-])'
                          'C(O1)COP(=O)([O-])[O-]'),
            'formula': dict(EmpiricalFormula('C39H46O25N15P4')),
            'mol_wt': 1248.772047992,
            'charge': -5,
            'warnings': None,
        })

    def test_get_bpform_properties_with_ph(self):
        client = rest.app.test_client()
        rv = client.post('/api/bpform/', json=dict(
            alphabet='dna',
            seq='ACGT',
            ph=7.,
            major_tautomer=True))
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json()['structure'], (
            'Cc1cn(C2CC(O)C(COP(=O)([O-])OC3CC(OC3COP(=O)([O-])OC3CC(OC3COP(=O)([O-])'
            'OC3CC(OC3COP(=O)([O-])[O-])n3cnc4c(N)ncnc34)n3ccc(N)nc3=O)n3cnc4c3nc(N)[nH]'
            'c4=O)O2)c(=O)[nH]c1=O'))
        self.assertEqual(rv.get_json()['alphabet'], 'dna')
        self.assertEqual(rv.get_json()['seq'], 'ACGT')
        self.assertEqual(rv.get_json()['length'], 4)

    def test_get_bpform_properties_circular(self):
        client = rest.app.test_client()
        rv = client.post('/api/bpform/', json=dict(alphabet='dna', seq='ACGT', circular=True))
        self.assertEqual(rv.status_code, 200)
        print(rv.get_json())
        self.assertEqual(rv.get_json(), {
            'alphabet': 'dna',
            'seq': 'ACGT',
            'length': 4,
            'structure': ('Nc1c2ncn(c2ncn1)C1CC2OP(=O)(OCC3C(OP(=O)(OCC4C(OP(=O)(OCC5C(OP(=O)'
                          '(OCC2O1)[O-])CC(n1cc(C)c(=O)[nH]c1=O)O5)[O-])CC(n1c2nc(N)[nH]c(=O)'
                          'c2nc1)O4)[O-])CC(n1c(=O)nc(N)cc1)O3)[O-]'),
            'formula': dict(EmpiricalFormula('C39H45O24N15P4')),
            'mol_wt': 1231.7650479919998,
            'charge': -4,
            'warnings': None,
        })

    def test_get_bpform_properties_errors(self):
        client = rest.app.test_client()

        rv = client.post('/api/bpform/', json=dict(alphabet='lipid', seq='ACGT'))
        self.assertEqual(rv.status_code, 400)

        rv = client.post('/api/bpform/', json=dict(alphabet='dna', seq='ACG[T'))
        self.assertEqual(rv.status_code, 400)

        rv = client.post('/api/bpform/', json=dict(alphabet='dna', seq='ACGT', ph='a'))
        self.assertEqual(rv.status_code, 400)

        rv = client.post('/api/bpform/', json=dict(alphabet='dna', seq='ACGT', ph=7., major_tautomer='a'))
        self.assertEqual(rv.status_code, 400)

        with mock.patch.object(core.BpForm, 'export', side_effect=Exception('unable to generate SMILES')):
            rv = client.post('/api/bpform/', json=dict(alphabet='dna', seq='ACGT'))
        self.assertEqual(rv.status_code, 200)

        with mock.patch.object(core.BpForm, 'get_major_micro_species', side_effect=Exception('unable to generate SMILES')):
            rv = client.post('/api/bpform/', json=dict(alphabet='dna', seq='ACGT', ph=7.))
        self.assertEqual(rv.status_code, 200)

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
