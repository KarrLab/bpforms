""" Test of bpforms.alphabet.protein

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Alphabet, Identifier, IdentifierSet, Monomer
from bpforms.alphabet import protein
from bpforms.util import validate_bpform_linkages
from wc_utils.util.chem import EmpiricalFormula, OpenBabelUtils
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
        self.tmp_pdbfile = os.path.join(self.dirname, 'AA0005.PDB')

        pdb = ''
        pdb += 'ATOM      1  N   Cys     1       0.000   0.000   0.000'+'\n'
        pdb += 'ATOM      2  HN  Cys     1      -0.560  -0.366   0.772'+'\n'
        pdb += 'ATOM      3  CA  Cys     1       1.454   0.000   0.000'+'\n'
        pdb += 'ATOM      4  HA  Cys     1       1.815  -0.552  -0.872'+'\n'
        pdb += 'ATOM      5  C   Cys     1       1.798  -0.764   1.252'+'\n'
        pdb += 'ATOM      6  O   Cys     1       0.970  -0.863   2.161'+'\n'
        pdb += 'ATOM      7  CB  Cys     1       1.979   1.438   0.000'+'\n'
        pdb += 'ATOM      8  HB1 Cys     1       1.580   1.957  -0.874'+'\n'
        pdb += 'ATOM      9  HB2 Cys     1       1.629   1.971   0.887'+'\n'
        pdb += 'ATOM     10  SG  Cys     1       3.784   1.560  -0.073'+'\n'
        pdb += 'ATOM     11  HG  Cys     1       4.169   1.730  -1.356'+'\n'

        with open(self.tmp_pdbfile, 'w') as file:
            file.write(pdb)

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

    def test_ProteinForm_get_monomer_isotope_structure(self):
        path = os.path.join(self.dirname, 'alphabet.yml')

        structure_isotope, indexn, indexc = protein.ProteinAlphabetBuilder().get_monomer_isotope_structure('AA0005', self.tmp_pdbfile)
        structure = protein.ProteinAlphabetBuilder().get_monomer_structure('AA0005', self.tmp_pdbfile)

        # check if correct isotopic species are present
        self.assertEqual(OpenBabelUtils.get_inchi(structure_isotope), 'InChI=1S/C3H7NOS/c4-3(1-5)2-6/h1,3,6H,2,4H2/t3-/m1/s1/i1+1,4+1')
        # check if correct index for N and C atoms
        self.assertEqual(indexn, 4)
        self.assertEqual(indexc, 1)
        # just in case check that original structure has not been modified
        self.assertEqual(OpenBabelUtils.get_inchi(structure), 'InChI=1S/C3H7NOS/c4-3(1-5)2-6/h1,3,6H,2,4H2/t3-/m1/s1')

    def test_ProteinAlphabetBuilder(self):
        path = os.path.join(self.dirname, 'alphabet.yml')

        alphabet = protein.ProteinAlphabetBuilder().run(path=path)
        self.assertEqual(alphabet.monomers.F.get_formula(), EmpiricalFormula('C9H11NO'))

        self.assertTrue(os.path.isfile(path))
        alphabet = Alphabet().from_yaml(path)
        self.assertEqual(alphabet.monomers.F.get_formula(), EmpiricalFormula('C9H11NO'))

    def test_validate_linkages(self):
        validate_bpform_linkages(protein.CanonicalProteinForm)
        validate_bpform_linkages(protein.ProteinForm)
