""" Test of bpforms.alphabet.dna

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.alphabet import dna
from wc_utils.util.chem import EmpiricalFormula
import os.path
import shutil
import tempfile
import unittest


class DnaTestCase(unittest.TestCase):
    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_dna_alphabet(self):
        self.assertEqual(dna.dna_alphabet.monomers.A.get_formula(), EmpiricalFormula('C10H12N5O6P'))
        self.assertEqual(dna.dna_alphabet.monomers.C.get_formula(), EmpiricalFormula('C9H12N3O7P'))
        self.assertEqual(dna.dna_alphabet.monomers.G.get_formula(), EmpiricalFormula('C10H12N5O7P'))
        self.assertEqual(dna.dna_alphabet.monomers.T.get_formula(), EmpiricalFormula('C10H13N2O8P'))

        self.assertEqual(dna.canonical_dna_alphabet.monomers.A.get_formula(), EmpiricalFormula('C10H12N5O6P'))
        self.assertEqual(dna.canonical_dna_alphabet.monomers.C.get_formula(), EmpiricalFormula('C9H12N3O7P'))
        self.assertEqual(dna.canonical_dna_alphabet.monomers.G.get_formula(), EmpiricalFormula('C10H12N5O7P'))
        self.assertEqual(dna.canonical_dna_alphabet.monomers.T.get_formula(), EmpiricalFormula('C10H13N2O8P'))

    def test_DnaForm_init(self):
        dna.DnaForm()
        dna.CanonicalDnaForm()

    def test_DnaForm_properties(self):
        monomers = dna.canonical_dna_alphabet.monomers
        form = dna.CanonicalDnaForm().from_str('ACGT')
        self.assertEqual(form.get_formula(), EmpiricalFormula('H') * -3
                         + monomers.A.get_formula()
                         + monomers.C.get_formula()
                         + monomers.G.get_formula()
                         + monomers.T.get_formula())
        self.assertEqual(form.get_formula(), EmpiricalFormula('C39H46O28N15P4'))
        self.assertEqual(form.get_charge(), -5)

    def test_DnaAlphabetBuilder(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        alphabet = dna.DnaAlphabetBuilder(_max_monomers=3).run(path=path)
        self.assertEqual(alphabet.monomers.A.get_formula(), EmpiricalFormula('C10H12N5O6P'))
        self.assertTrue(os.path.isfile(path))
