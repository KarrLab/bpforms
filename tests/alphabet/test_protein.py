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
import openbabel
import os.path
import pkg_resources
import shutil
import requests
import tempfile
import unittest

ALA_smiles = 'CC([NH3+])C([O-])=O'
di_ALA_smiles = 'CC([NH3+])C(=O)[NH2+]C(C)C([O-])=O'
tri_ALA_smiles = 'CC([NH3+])C(=O)[NH2+]C(C)C(=O)[NH2+]C(C)C([O-])=O'

ALA_inchi = 'InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1'
di_ALA_inchi = 'InChI=1S/C6H12N2O3/c1-3(7)5(9)8-4(2)6(10)11/h3-4H,7H2,1-2H3,(H,8,9)(H,10,11)/t3-,4?/m0/s1'
tri_ALA_inchi = 'InChI=1S/C9H17N3O4/c1-4(10)7(13)11-5(2)8(14)12-6(3)9(15)16/h4-6H,10H2,1-3H3,(H,11,13)(H,12,14)(H,15,16)/t4-,5?,6?/m0/s1'


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

    def test_ProteinAlphabetBuilder_is_termini(self):
        builder = protein.ProteinAlphabetBuilder()

        mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        conv.SetInFormat('smiles')
        conv.ReadString(mol, 'C[C@H]([NH3+])C=O')

        for i_atom in range(mol.NumAtoms()):
            if i_atom + 1 == 4:
                self.assertTrue(builder.is_n_terminus(mol.GetAtom(i_atom + 1)))
            else:
                self.assertFalse(builder.is_n_terminus(mol.GetAtom(i_atom + 1)))

            if i_atom + 1 == 8:
                self.assertTrue(builder.is_c_terminus(mol.GetAtom(i_atom + 1)))
            else:
                self.assertFalse(builder.is_c_terminus(mol.GetAtom(i_atom + 1)))

        self.assertFalse(builder.is_n_terminus(None))

        conv.ReadString(mol, 'N')
        atom = mol.GetAtom(1)
        self.assertFalse(builder.is_n_terminus(atom))

        conv.ReadString(mol, '[NH3+]O')
        atom = mol.GetAtom(1)
        self.assertFalse(builder.is_n_terminus(atom))

        conv.ReadString(mol, '[NH4+]')
        atom = mol.GetAtom(1)
        self.assertFalse(builder.is_n_terminus(atom))

        conv.ReadString(mol, 'C[NH+](C)C')
        atom = mol.GetAtom(2)
        self.assertFalse(builder.is_n_terminus(atom))

        self.assertFalse(builder.is_c_terminus(None))

        conv.ReadString(mol, 'C[C@H]([NH3+])C=O')
        atom = mol.GetAtom(8)
        atom.SetFormalCharge(1)
        self.assertFalse(builder.is_c_terminus(atom))

        conv.ReadString(mol, 'CC(C)O')
        atom = mol.GetAtom(2)
        self.assertFalse(builder.is_c_terminus(atom))

        conv.ReadString(mol, 'CCO')
        atom = mol.GetAtom(2)
        self.assertFalse(builder.is_c_terminus(atom))

        conv.ReadString(mol, 'C=C=O')
        atom = mol.GetAtom(2)
        self.assertFalse(builder.is_c_terminus(atom))

    def test_ProteinAlphabetBuilder_is_termini_2(self):
        builder = protein.ProteinAlphabetBuilder()

        mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        conv.SetInFormat('smiles')
        conv.ReadString(mol, 'O=C[C@H](CCCCNC(=O)[C@H](CCCC[NH3+])[NH3+])[NH3+]')

        atom = mol.GetAtom(26)
        self.assertTrue(builder.is_n_terminus(atom))

    def test_test_ProteinAlphabetBuilder_is_terminus(self):
        builder = protein.ProteinAlphabetBuilder()

        mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        conv.SetInFormat('smiles')
        conv.ReadString(mol, 'O=C[C@H](CCCCNC(=O)[C@H](CCCC[NH3+])[NH3+])[NH3+]')

        self.assertTrue(builder.is_terminus(mol.GetAtom(26), mol.GetAtom(2)))
        self.assertFalse(builder.is_terminus(mol.GetAtom(26), mol.GetAtom(10)))

    def test_ProteinAlphabetBuilder_get_monomer_details(self):
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

    def test_ProteinAlphabetBuilder_get_monomer_structure(self):
        path = os.path.join(self.dirname, 'alphabet.yml')

        structure, index_n, index_c = protein.ProteinAlphabetBuilder().get_monomer_structure(
            'AA0005', self.tmp_pdbfile, ph=7.4, major_tautomer=True)

        # check if correct index for N and C atoms
        self.assertEqual(index_n, 5)
        self.assertEqual(index_c, 9)
        # just in case check that original structure has not been modified
        self.assertEqual(OpenBabelUtils.get_inchi(structure),
                         'InChI=1S/C3H7NOS/c4-3(1-5)2-6/h1,3,6H,2,4H2/p+1/t3-/m1/s1')

    def test_ProteinAlphabetBuilder(self):
        pdb_dir = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.pdb'))
        if os.path.isdir(pdb_dir):
            shutil.rmtree(pdb_dir)

        path = os.path.join(self.dirname, 'alphabet.yml')

        alphabet = protein.ProteinAlphabetBuilder().run(path=path)
        self.assertEqual(alphabet.monomers.F.get_formula(), EmpiricalFormula('C9H12NO'))

        self.assertTrue(os.path.isfile(path))
        alphabet = Alphabet().from_yaml(path)
        self.assertEqual(alphabet.monomers.F.get_formula(), EmpiricalFormula('C9H12NO'))

    def test_validate_canonical_linkages(self):
        validate_bpform_linkages(protein.CanonicalProteinForm)
        validate_bpform_linkages(protein.ProteinForm)

        self.validate_alphabet(protein.CanonicalProteinForm().alphabet)
        self.validate_alphabet(protein.ProteinForm().alphabet)

    def validate_alphabet(self, alphabet):
        errors = []
        builder = protein.ProteinAlphabetBuilder()
        for monomer in alphabet.monomers.values():
            if monomer.left_bond_atoms:
                atom_n = monomer.structure.GetAtom(monomer.left_bond_atoms[0].position)            
                if not builder.is_n_terminus(atom_n):
                    errors.append('Monomer {} does not have a N-terminus'.format(monomer.id))

            atom_c = monomer.structure.GetAtom(monomer.right_bond_atoms[0].position)
            if not builder.is_c_terminus(atom_c):
                errors.append('Monomer {} does not have a C-terminus'.format(monomer.id))
            
            if monomer.left_bond_atoms:
                if not builder.is_terminus(atom_n, atom_c):
                    errors.append('Monomer {} does not have termini'.format(monomer.id))
        if errors:
            raise ValueError('Alphabet has invalid monomer(s):\n  {}'.format('\n  '.join(
                errors)))

    def test_validate_form(self):
        form = protein.ProteinForm()
        form.from_str('ARG')
        self.assertEqual(form.validate(), [])
