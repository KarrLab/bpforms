""" Test of bpforms.alphabet.protein

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Alphabet, Identifier, IdentifierSet, Monomer
from bpforms.alphabet import protein
from wc_utils.util.chem import EmpiricalFormula
import os.path
import shutil
import requests
import tempfile
import unittest


class ProteinTestCase(unittest.TestCase):
    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_protein_alphabet(self):
        self.assertEqual(protein.protein_alphabet.monomers.F.get_formula(), EmpiricalFormula('C9H12NO'))
        self.assertEqual(protein.canonical_protein_alphabet.monomers.F.get_formula(), EmpiricalFormula('C9H11NO'))

    def test_ProteinForm_init(self):
        protein.ProteinForm()
        protein.CanonicalProteinForm()

    def test_ProteinForm_properties(self):
        monomers = protein.canonical_protein_alphabet.monomers
        form = protein.CanonicalProteinForm().from_str('AAA')
        self.assertEqual(form.get_formula(), EmpiricalFormula('O') * 3 + EmpiricalFormula('H2O') * -2
                         + monomers.A.get_formula()
                         + monomers.A.get_formula()
                         + monomers.A.get_formula())
        self.assertEqual(form.get_formula(), EmpiricalFormula('C9H17N3O4'))
        self.assertEqual(form.get_charge(), 3 * monomers.A.get_charge())

    def test_ProteinForm_get_monomer_details(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        session = requests.Session()
        self.assertEqual(protein.ProteinAlphabetBuilder().get_monomer_details('AA0037', session)[0], 'S')
        self.assertEqual(protein.ProteinAlphabetBuilder().get_monomer_details('AA0037', session)[3], {'AA0016'})
        # add more for the identifiers
        # self.assertEqual(protein.ProteinAlphabetBuilder().get_monomer_details('AA0037', session)[2],
                         # IdentifierSet([Identifier('mod', 'MOD:00046')]))
 
    def test_ProteinAlphabetBuilder(self):
        path = os.path.join(self.dirname, 'alphabet.yml')

        alphabet = protein.ProteinAlphabetBuilder().run(path=path)
        self.assertEqual(alphabet.monomers.F.get_formula(), EmpiricalFormula('C9H11NO'))

        self.assertTrue(os.path.isfile(path))
        alphabet = Alphabet().from_yaml(path)
        self.assertEqual(alphabet.monomers.F.get_formula(), EmpiricalFormula('C9H11NO'))
