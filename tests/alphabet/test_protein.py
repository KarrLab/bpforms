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

    def test_protein_alphabet(self):
        self.assertEqual(protein.protein_alphabet.bases.F.get_formula(), EmpiricalFormula('C9N1O2H11'))
        self.assertEqual(protein.canonical_protein_alphabet.bases.F.get_formula(), EmpiricalFormula('C9N1O2H11'))

    def test_ProteinForm_init(self):
        protein.ProteinForm()
        protein.CanonicalProteinForm()

    def test_ProteinForm_properties(self):
        bases = protein.canonical_protein_alphabet.bases
        form = protein.CanonicalProteinForm().from_str('AAA')
        self.assertEqual(form.get_formula(), EmpiricalFormula('H2') * -2
                         + bases.A.get_formula()
                         + bases.A.get_formula()
                         + bases.A.get_formula())
        self.assertEqual(form.get_charge(), 3 * bases.A.get_charge())

    def test_ProteinAlphabetBuilder(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        alphabet = protein.ProteinAlphabetBuilder().run(path=path)
        self.assertEqual(alphabet.bases.F.get_formula(), EmpiricalFormula('C9N1O2H11'))
        self.assertTrue(os.path.isfile(path))
