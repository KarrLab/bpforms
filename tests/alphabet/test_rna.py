""" Test of bpforms.alphabet.rna

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.alphabet import rna
from wc_utils.util.chem import EmpiricalFormula
import os.path
import shutil
import tempfile
import unittest


class RnaTestCase(unittest.TestCase):
    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_rna_alphabet(self):
        self.assertEqual(rna.rna_alphabet.A.get_formula(), EmpiricalFormula('C10H12N5O7P'))
        self.assertEqual(rna.rna_alphabet.C.get_formula(), EmpiricalFormula('C9H12N3O8P'))
        self.assertEqual(rna.rna_alphabet.G.get_formula(), EmpiricalFormula('C10H12N5O8P'))
        self.assertEqual(rna.rna_alphabet.U.get_formula(), EmpiricalFormula('C9H11N2O9P'))

    def test_RnaForm_init(self):
        rna.RnaForm()

    def test_RnaAlphabetBuilder(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        alphabet = rna.RnaAlphabetBuilder().run(path=path)
        self.assertEqual(alphabet.A.get_formula(), EmpiricalFormula('C10H12N5O7P'))
        self.assertTrue(os.path.isfile(path))
