""" Test of bpforms.rna

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms import rna
from wc_utils.util.chem import EmpiricalFormula
import unittest


class RnaTestCase(unittest.TestCase):
    def test_rna_alphabet(self):
        self.assertEqual(rna.rna_alphabet.A.get_formula(), EmpiricalFormula('C10H12N5O7P'))
        self.assertEqual(rna.rna_alphabet.C.get_formula(), EmpiricalFormula('C9H12N3O8P'))
        self.assertEqual(rna.rna_alphabet.G.get_formula(), EmpiricalFormula('C10H12N5O8P'))
        self.assertEqual(rna.rna_alphabet.U.get_formula(), EmpiricalFormula('C9H11N2O9P'))

    def test_RnaForm_init(self):
        rna.RnaForm()
