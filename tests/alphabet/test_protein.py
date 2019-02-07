""" Test of bpforms.alphabet.protein

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.alphabet import protein
from wc_utils.util.chem import EmpiricalFormula
import os.path
import shutil
import tempfile
import unittest


class ProteinTestCase(unittest.TestCase):
    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    @unittest.skip('Todo')
    def test_protein_alphabet(self):
        pass

    def test_ProteinForm_init(self):
        protein.ProteinForm()
        protein.CanonicalProteinForm()

    def test_ProteinAlphabetBuilder(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        alphabet = protein.ProteinAlphabetBuilder().run(path=path)        
        self.assertEqual(alphabet.bases.A.get_formula(), EmpiricalFormula('C3H7NO2'))
        self.assertTrue(os.path.isfile(path))
