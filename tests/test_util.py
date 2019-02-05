""" Test of bpforms.util

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms import dna
from bpforms import rna
from bpforms import protein
from bpforms import util
import unittest


class UtilTestCase(unittest.TestCase):
    def test_get_form(self):
        self.assertEqual(util.get_form('dna'), dna.DnaForm)
        self.assertEqual(util.get_form('rna'), rna.RnaForm)
        self.assertEqual(util.get_form('protein'), protein.ProteinForm)
        with self.assertRaises(ValueError):
            util.get_form('lipid')
