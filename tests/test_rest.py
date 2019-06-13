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
import imghdr
import importlib
import mock
import os
import shutil
import tempfile
import unittest
import warnings


class RestTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        shutil.rmtree(rest.cache_dir)
        importlib.reload(rest)

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
            'seq': 'ACGT | circular',
            'length': 4,
            'structure': ('Nc1c2ncn(c2ncn1)C1CC2OP(=O)(OCC3C(OP(=O)(OCC4C(OP(=O)(OCC5C(OP(=O)'
                          '(OCC2O1)[O-])CC(n1cc(C)c(=O)[nH]c1=O)O5)[O-])CC(n1c2nc(N)[nH]c(=O)'
                          'c2nc1)O4)[O-])CC(n1c(=O)nc(N)cc1)O3)[O-]'),
            'formula': dict(EmpiricalFormula('C39H45O24N15P4')),
            'mol_wt': 1231.7650479919998,
            'charge': -4,
            'warnings': None,
        })

    def test_get_bpform_properties_no_structure(self):
        client = rest.app.test_client()
        with mock.patch.object(core.BpForm, 'get_structure', return_value=None):
            rv = client.post('/api/bpform/', json=dict(alphabet='dna', seq='ACGT'))
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            'alphabet': 'dna',
            'seq': 'ACGT',
            'length': 4,
            'structure': None,
            'formula': dict(EmpiricalFormula('C39H46O25N15P4')),
            'mol_wt': 1248.772047992,
            'charge': -5,
            'warnings': None,
        })

    def test_get_bpform_properties_warnings(self):
        client = rest.app.test_client()

        def side_effect(*args):
            warnings.warn('My warning', UserWarning)
            return None
        with mock.patch.object(core.BpForm, 'get_structure', side_effect=side_effect):
            rv = client.post('/api/bpform/', json=dict(alphabet='dna', seq='ACGT'))
        self.assertEqual(rv.status_code, 200)
        self.assertIn('My warning', rv.get_json()['warnings'])

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

        rv = client.post('/api/bpform/', json=dict(alphabet='dna',
                                                   seq=('ACGT[id: "dI"'
                                                        ' | structure: "{}"'
                                                        ' | backbone-bond-atom: N10'
                                                        ' | backbone-displaced-atom: H10 ]').format(
                                                       'O=C1NC=NC2=C1N=CN2'
                                                   )))
        self.assertEqual(rv.status_code, 200)

        rv = client.post('/api/bpform/', json=dict(alphabet='dna',
                                                   seq=('ACGT[id: "dI"'
                                                        ' | structure: "{}"'
                                                        ' | backbone-displaced-atom: H10 ]').format(
                                                       'O=C1NC=NC2=C1N=CN2'
                                                   )))
        self.assertEqual(rv.status_code, 400)

    def test_get_alphabets(self):
        client = rest.app.test_client()
        rv = client.get('/api/alphabet/')
        self.assertEqual(rv.status_code, 200)
        alphabets = rv.get_json()
        self.assertEqual(alphabets['dna']['name'], 'DNA nucleobases')
        self.assertEqual(alphabets['canonical_dna']['name'], 'Canonical DNA nucleobases')

    def test_get_alphabet(self):
        client = rest.app.test_client()

        rv = client.get('/api/alphabet/dna/')
        self.assertEqual(rv.status_code, 200)
        rv_json = rv.get_json()
        self.assertTrue(rv_json['monomers']['A']['binds_backbone'])
        self.assertTrue(rv_json['monomers']['A']['binds_left'])
        self.assertTrue(rv_json['monomers']['A']['binds_right'])
        self.assertEqual(rv_json['monomers']['A']['mol_wt'], 135.13)
        self.assertEqual(rv_json['monomers']['A']['formula'], {'C': 5.0, 'H': 5.0, 'N': 5.0})
        self.assertEqual(rv_json['monomers']['A']['charge'], 0)
        alphabet = core.Alphabet().from_dict(rv_json)
        self.assertTrue(dna.dna_alphabet.is_equal(alphabet))

        rv = client.get('/api/alphabet/lipid/')
        self.assertEqual(rv.status_code, 400)

    @unittest.skipIf(os.getenv('CIRCLECI', '0') in ['1', 'true'], 'Skip long test in CircleCI')
    def test_get_alphabet_caching(self):
        client = rest.app.test_client()

        rv = client.get('/api/alphabet/protein/')
        self.assertEqual(rv.status_code, 200)

        rv = client.get('/api/alphabet/protein/')
        self.assertEqual(rv.status_code, 200)

    def test_get_monomer(self):
        client = rest.app.test_client()

        rv = client.get('/api/alphabet/dna/A/')
        self.assertEqual(rv.status_code, 200)
        rv_json = rv.get_json()
        self.assertTrue(rv_json['binds_backbone'])
        self.assertTrue(rv_json['binds_left'])
        self.assertTrue(rv_json['binds_right'])
        self.assertEqual(rv_json['mol_wt'], 135.13)
        self.assertEqual(rv_json['formula'], {'C': 5.0, 'H': 5.0, 'N': 5.0})
        self.assertEqual(rv_json['charge'], 0)

        rv = client.get('/api/alphabet/dna/A/svg/')
        self.assertEqual(rv.status_code, 200)
        self.assertTrue(rv.data.decode('utf-8').startswith('<?xml'))

        rv = client.get('/api/alphabet/dna/A/png/')
        self.assertEqual(rv.status_code, 200)
        png_file, png_filename = tempfile.mkstemp()
        os.close(png_file)
        with open(png_filename, 'wb') as file:
            file.write(rv.data)
        self.assertEqual(imghdr.what(png_filename), 'png')
        os.remove(png_filename)

        rv = client.get('/api/alphabet/lipid/A/')
        self.assertEqual(rv.status_code, 400)

        rv = client.get('/api/alphabet/dna/XXXX/')
        self.assertEqual(rv.status_code, 400)

        rv = client.get('/api/alphabet/dna/A/unknown/')
        self.assertEqual(rv.status_code, 400)
