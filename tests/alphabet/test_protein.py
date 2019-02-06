""" Test of bpforms.alphabet.protein

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.alphabet import protein
from wc_utils.util.chem import EmpiricalFormula
import unittest


class ProteinTestCase(unittest.TestCase):
    @unittest.skip('Todo')
    def test_protein_alphabet(self):
        pass

    def test_ProteinForm_init(self):
        protein.ProteinForm()
