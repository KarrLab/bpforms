""" Test of bpforms.alphabet.rna

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import IdentifierSet
from bpforms.alphabet import rna
from wc_utils.util.chem import EmpiricalFormula
import mock
import os.path
import requests
import shutil
import tempfile
import unittest


class RnaTestCase(unittest.TestCase):
    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_rna_alphabet(self):
        self.assertEqual(rna.rna_alphabet.monomers.A.get_formula(), EmpiricalFormula('C10O4N5H13'))
        self.assertEqual(rna.rna_alphabet.monomers.C.get_formula(), EmpiricalFormula('C9O5N3H13'))
        self.assertEqual(rna.rna_alphabet.monomers.G.get_formula(), EmpiricalFormula('C10O5N5H13'))
        self.assertEqual(rna.rna_alphabet.monomers.U.get_formula(), EmpiricalFormula('C9O6N2H12'))

        self.assertEqual(rna.canonical_rna_alphabet.monomers.A.get_formula(), EmpiricalFormula('C10O4N5H13'))
        self.assertEqual(rna.canonical_rna_alphabet.monomers.C.get_formula(), EmpiricalFormula('C9O5N3H13'))
        self.assertEqual(rna.canonical_rna_alphabet.monomers.G.get_formula(), EmpiricalFormula('C10O5N5H13'))
        self.assertEqual(rna.canonical_rna_alphabet.monomers.U.get_formula(), EmpiricalFormula('C9O6N2H12'))

    def test_RnaForm_init(self):
        rna.RnaForm()
        rna.CanonicalRnaForm()

    def test_RnaForm_properties(self):
        monomers = rna.canonical_rna_alphabet.monomers
        form = rna.CanonicalRnaForm().from_str('ACGU')
        self.assertEqual(form.get_formula(), EmpiricalFormula('H-1O3P') * 4 + EmpiricalFormula('OH') * -3
                         + monomers.A.get_formula()
                         + monomers.C.get_formula()
                         + monomers.G.get_formula()
                         + monomers.U.get_formula())
        self.assertEqual(form.get_formula(), EmpiricalFormula('C38H44O29N15P4'))
        self.assertEqual(form.get_charge(), -5)

    def test_RnaAlphabetBuilder_list_monomers(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        with mock.patch.object(rna.RnaAlphabetBuilder, 'get_monomer_details', return_value=(None, IdentifierSet())):
            alphabet = rna.RnaAlphabetBuilder().run(path=path)
        self.assertEqual(alphabet.monomers.A.get_formula(), EmpiricalFormula('C10O4N5H13'))
        self.assertTrue(os.path.isfile(path))

    def test_RnaAlphabetBuilder_get_monomer_details(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        session = requests.Session()
        self.assertEqual(rna.RnaAlphabetBuilder().get_monomer_details('(pN)2′3′>p', session), (None, IdentifierSet()))
        self.assertEqual(rna.RnaAlphabetBuilder().get_monomer_details('Xm', session), (None, IdentifierSet()))

    def test_RnaAlphabetBuilder(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        alphabet = rna.RnaAlphabetBuilder(_max_monomers=3).run(path=path)
        self.assertEqual(alphabet.monomers.A.get_formula(), EmpiricalFormula('C10O4N5H13'))
        self.assertTrue(os.path.isfile(path))
