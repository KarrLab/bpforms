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

ALA_smiles = 'CC([NH3+])C([O-])=O'
di_ALA_smiles = 'CC([NH3+])C(=O)[NH2+]C(C)C([O-])=O'
tri_ALA_smiles = 'CC([NH3+])C(=O)[NH2+]C(C)C(=O)[NH2+]C(C)C([O-])=O'

ALA_inchi = 'InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)'
di_ALA_inchi = 'InChI=1S/C6H12N2O3/c1-3(7)5(9)8-4(2)6(10)11/h3-4H,7H2,1-2H3,(H,8,9)(H,10,11)/p+1'
tri_ALA_inchi = 'InChI=1S/C9H17N3O4/c1-4(10)7(13)11-5(2)8(14)12-6(3)9(15)16/h4-6H,10H2,1-3H3,(H,11,13)(H,12,14)(H,15,16)/p+2'


class ProteinTestCase(unittest.TestCase):
    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_protein_alphabet(self):
        self.assertEqual(protein.protein_alphabet.monomers.F.get_formula(), EmpiricalFormula('C9H12NO'))
        self.assertEqual(protein.canonical_protein_alphabet.monomers.F.get_formula(), EmpiricalFormula('C9H12NO'))

    def test_ProteinForm_init(self):
        protein.ProteinForm()
        protein.CanonicalProteinForm()

    def test_ProteinForm_properties(self):
        monomers = protein.canonical_protein_alphabet.monomers

        form = protein.CanonicalProteinForm().from_str('A')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C3H7NO2'))
        self.assertEqual(form.get_charge(), 0)        
        self.assertEqual(form.export('inchi'), ALA_inchi)

        form = protein.CanonicalProteinForm().from_str('AA')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C6H12N2O3'))
        self.assertEqual(form.get_charge(), 0)
        self.assertEqual(form.export('inchi'), di_ALA_inchi)

        form = protein.CanonicalProteinForm().from_str('AAA')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C9H17N3O4'))
        self.assertEqual(form.get_charge(), 0)
        self.assertEqual(form.export('inchi'), tri_ALA_inchi)

        form = protein.ProteinForm().from_str('AAA')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C9H17N3O4'))
        self.assertEqual(form.get_charge(), 0)
        self.assertEqual(form.export('inchi'), tri_ALA_inchi)

    def test_ProteinForm_get_monomer_details(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        session = requests.Session()
        self.assertEqual(protein.ProteinAlphabetBuilder().get_monomer_details('AA0037', session)[0], 'S')
        self.assertEqual(protein.ProteinAlphabetBuilder().get_monomer_details('AA0037', session)[3], {'AA0016'})

        identifiers = protein.ProteinAlphabetBuilder().get_monomer_details('AA0037', session)[2]
        identifiers_test = set([
            Identifier('pdb.ligand', 'SEP'),
            Identifier('mod', 'MOD:00046'),
            Identifier('go', 'GO:0018105'),
            Identifier('cas', '407-41-0'),
            Identifier('chebi', 'CHEBI:45522')
        ])
        self.assertEqual(identifiers, identifiers_test)

    def test_ProteinAlphabetBuilder(self):
        path = os.path.join(self.dirname, 'alphabet.yml')

        alphabet = protein.ProteinAlphabetBuilder().run(path=path)
        self.assertEqual(alphabet.monomers.F.get_formula(), EmpiricalFormula('C9H11NO'))

        self.assertTrue(os.path.isfile(path))
        alphabet = Alphabet().from_yaml(path)
        self.assertEqual(alphabet.monomers.F.get_formula(), EmpiricalFormula('C9H11NO'))
