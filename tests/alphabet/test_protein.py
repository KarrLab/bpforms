""" Test of bpforms.alphabet.protein

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Alphabet, Identifier, IdentifierSet, Monomer
from bpforms.alphabet import protein
from bpforms.util import validate_bpform_bonds
from wc_utils.util.chem import EmpiricalFormula, OpenBabelUtils
import openbabel
import os.path
import pkg_resources
import shutil
import requests
import tempfile
import unittest

ALA_smiles = 'C[C@H]([NH3+])C(=O)O'
di_ALA_smiles = 'C[C@H]([NH3+])C(=O)N[C@@H](C)C(=O)O'
tri_ALA_smiles = 'C[C@H]([NH3+])C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O'


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
        self.assertEqual(protein.protein_alphabet.monomers.F.get_formula(),
                         EmpiricalFormula('C9H12NO') + EmpiricalFormula('O'))
        self.assertEqual(protein.canonical_protein_alphabet.monomers.F.get_formula(),
                         EmpiricalFormula('C9H12NO') + EmpiricalFormula('O'))

    def test_ProteinForm_init(self):
        protein.ProteinForm()
        protein.CanonicalProteinForm()

    def test_ProteinForm_properties(self):
        monomers = protein.canonical_protein_alphabet.monomers

        form = protein.CanonicalProteinForm().from_str('A')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C3H7NO2') + EmpiricalFormula('H'))
        self.assertEqual(form.get_charge(), 1)
        self.assertEqual(form.export('smiles'), ALA_smiles)

        form = protein.CanonicalProteinForm().from_str('AA')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C6H12N2O3') + EmpiricalFormula('H'))
        self.assertEqual(form.get_charge(), 1)
        self.assertEqual(form.export('smiles'), di_ALA_smiles)

        form = protein.CanonicalProteinForm().from_str('AAA')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C9H17N3O4') + EmpiricalFormula('H'))
        self.assertEqual(form.get_charge(), 1)
        self.assertEqual(form.export('smiles'), tri_ALA_smiles)

        form = protein.ProteinForm().from_str('AAA')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C9H17N3O4') + EmpiricalFormula('H'))
        self.assertEqual(form.get_charge(), 1)
        self.assertEqual(form.export('smiles'), tri_ALA_smiles)

    def test_ProteinAlphabetBuilder_is_termini(self):
        builder = protein.ProteinAlphabetBuilder()

        mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        conv.SetInFormat('smiles')
        conv.ReadString(mol, 'C[C@H]([NH3+])C=O')

        for i_atom in range(mol.NumAtoms()):
            if i_atom + 1 == 4:
                self.assertTrue(builder.is_n_terminus(mol, mol.GetAtom(i_atom + 1)))
            else:
                self.assertFalse(builder.is_n_terminus(mol, mol.GetAtom(i_atom + 1)))

            if i_atom + 1 == 8:
                self.assertTrue(builder.is_c_terminus(mol, mol.GetAtom(i_atom + 1)))
            else:
                self.assertFalse(builder.is_c_terminus(mol, mol.GetAtom(i_atom + 1)))

        self.assertFalse(builder.is_n_terminus(mol, None))

        conv.ReadString(mol, 'N')
        atom = mol.GetAtom(1)
        self.assertFalse(builder.is_n_terminus(mol, atom))

        conv.ReadString(mol, '[NH3+]O')
        atom = mol.GetAtom(1)
        self.assertFalse(builder.is_n_terminus(mol, atom))

        conv.ReadString(mol, '[NH4+]')
        atom = mol.GetAtom(1)
        self.assertFalse(builder.is_n_terminus(mol, atom))

        conv.ReadString(mol, 'C[NH+](C)C')
        atom = mol.GetAtom(2)
        self.assertFalse(builder.is_n_terminus(mol, atom))

        self.assertFalse(builder.is_c_terminus(mol, None))

        conv.ReadString(mol, 'C[C@H]([NH3+])C=O')
        atom = mol.GetAtom(8)
        atom.SetFormalCharge(1)
        self.assertFalse(builder.is_c_terminus(mol, atom))

        conv.ReadString(mol, 'CC(C)O')
        atom = mol.GetAtom(2)
        self.assertFalse(builder.is_c_terminus(mol, atom))

        conv.ReadString(mol, 'CCO')
        atom = mol.GetAtom(2)
        self.assertFalse(builder.is_c_terminus(mol, atom))

        conv.ReadString(mol, 'C=C=O')
        atom = mol.GetAtom(2)
        self.assertFalse(builder.is_c_terminus(mol, atom))

    def test_ProteinAlphabetBuilder_is_termini_2(self):
        builder = protein.ProteinAlphabetBuilder()

        mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        conv.SetInFormat('smiles')
        conv.ReadString(mol, 'O=C[C@H](CCCCNC(=O)[C@H](CCCC[NH3+])[NH3+])[NH3+]')

        atom = mol.GetAtom(26)
        self.assertTrue(builder.is_n_terminus(mol, atom))

    def test_test_ProteinAlphabetBuilder_is_terminus(self):
        builder = protein.ProteinAlphabetBuilder()

        mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        conv.SetInFormat('smiles')
        conv.ReadString(mol, 'O=C[C@H](CCCCNC(=O)[C@H](CCCC[NH3+])[NH3+])[NH3+]')

        self.assertTrue(builder.is_terminus(mol.GetAtom(26), mol.GetAtom(2)))
        self.assertFalse(builder.is_terminus(mol.GetAtom(26), mol.GetAtom(10)))

    def test_ProteinAlphabetBuilder_get_resid_monomer_details(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        session = requests.Session()
        self.assertEqual(protein.ProteinAlphabetBuilder().get_resid_monomer_details('AA0037', session)[0], 'S')
        self.assertEqual(protein.ProteinAlphabetBuilder().get_resid_monomer_details('AA0037', session)[3], {'AA0016'})

        identifiers = protein.ProteinAlphabetBuilder().get_resid_monomer_details('AA0037', session)[2]
        identifiers_test = set([
            Identifier('pdb.ligand', 'SEP'),
            Identifier('mod', 'MOD:00046'),
            Identifier('go', 'GO:0018105'),
            Identifier('cas', '407-41-0'),
            Identifier('chebi', 'CHEBI:45522'),
            Identifier('resid', 'AA0037'),
        ])
        self.assertEqual(identifiers, identifiers_test)

    def test_ProteinAlphabetBuilder_get_resid_monomer_structure(self):
        path = os.path.join(self.dirname, 'alphabet.yml')

        structure, index_n, index_c = protein.ProteinAlphabetBuilder().get_resid_monomer_structure(
            'AA0005', self.tmp_pdbfile, ph=7.4, major_tautomer=True)

        # just in case check that original structure has not been modified
        self.assertEqual(OpenBabelUtils.export(structure, 'smiles'),
                         'OC(=O)[C@@H]([NH3+])CS')

        # check if correct index for N and C atoms
        self.assertEqual(index_n, 6)
        self.assertEqual(index_c, 2)

    def test_ProteinAlphabetBuilder(self):
        pdb_dir = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.pdb'))
        if os.path.isdir(pdb_dir):
            shutil.rmtree(pdb_dir)

        path = os.path.join(self.dirname, 'alphabet.yml')

        alphabet = protein.ProteinAlphabetBuilder(_max_monomers=10).run(ph=7.4, path=path)
        self.assertEqual(alphabet.monomers.F.get_formula(),
                         EmpiricalFormula('C9H12NO') + EmpiricalFormula('O'))

        self.assertTrue(os.path.isfile(path))
        alphabet = Alphabet().from_yaml(path)
        self.assertEqual(alphabet.monomers.F.get_formula(),
                         EmpiricalFormula('C9H12NO') + EmpiricalFormula('O'))

    def test_validate_bonds(self):
        validate_bpform_bonds(protein.CanonicalProteinForm)
        validate_bpform_bonds(protein.ProteinForm)

        self.validate_alphabet(protein.CanonicalProteinForm().alphabet)
        self.validate_alphabet(protein.ProteinForm().alphabet)

    def validate_alphabet(self, alphabet):
        errors = []
        builder = protein.ProteinAlphabetBuilder()
        for monomer in alphabet.monomers.values():
            if monomer.l_bond_atoms:
                atom_n = monomer.structure.GetAtom(monomer.l_bond_atoms[0].position)
                if not builder.is_n_terminus(monomer.structure, atom_n):
                    errors.append('Monomer {} does not have a N-terminus'.format(monomer.id))

            if monomer.r_bond_atoms:
                atom_c = monomer.structure.GetAtom(monomer.r_bond_atoms[0].position)
                if not builder.is_c_terminus(monomer.structure, atom_c, residue=False):
                    errors.append('Monomer {} does not have a C-terminus'.format(monomer.id))

            # if monomer.l_bond_atoms and monomer.r_bond_atoms:
            #    if not builder.is_terminus(atom_n, atom_c):
            #        errors.append('Monomer {} does not have termini'.format(monomer.id))
        if errors:
            raise ValueError('Alphabet has invalid monomer(s):\n  {}'.format('\n  '.join(
                errors)))

    def test_validate_form(self):
        form = protein.ProteinForm()
        form.from_str('ARG')
        self.assertEqual(form.validate(), [])

    def test_n_terminus_only(self):
        # AA0318
        form = protein.ProteinForm().from_str('{AA0318}')
        self.assertEqual(form.validate(), [])
        self.assertEqual(form.export('smiles'), form.seq[0].export('smiles'))
        self.assertEqual(form.get_formula(), form.seq[0].get_formula())
        self.assertEqual(form.get_charge(), form.seq[0].get_charge())

        form = protein.ProteinForm().from_str('A{AA0318}')
        self.assertEqual(form.validate(), [])

        form = protein.ProteinForm().from_str('{AA0318}A')
        self.assertNotEqual(form.validate(), [])

    def test_c_terminus_only(self):
        # AA0062
        form = protein.ProteinForm().from_str('{AA0062}')
        self.assertEqual(form.validate(), [])
        self.assertEqual(form.export('smiles'), form.seq[0].export('smiles'))
        self.assertEqual(form.get_formula(), form.seq[0].get_formula())
        self.assertEqual(form.get_charge(), form.seq[0].get_charge())

        form = protein.ProteinForm().from_str('A{AA0062}')
        self.assertNotEqual(form.validate(), [])

        form = protein.ProteinForm().from_str('{AA0062}A')
        self.assertEqual(form.validate(), [])

    def test_no_termini(self):
        # AA0344
        monomer = Monomer(structure='OC[C@H](N1SC[C@@H](C1=O)[NH3+])C=O')
        form = protein.ProteinForm()
        form.seq.append(monomer)
        self.assertEqual(form.validate(), [])
        self.assertEqual(form.export('smiles'), form.seq[0].export('smiles'))
        self.assertEqual(form.get_formula(), form.seq[0].get_formula())
        self.assertEqual(form.get_charge(), form.seq[0].get_charge())

        form = protein.ProteinForm()
        form.seq.append(monomer)
        form.seq.append(protein.protein_alphabet.monomers.A)
        self.assertNotEqual(form.validate(), [])

        form = protein.ProteinForm()
        form.seq.append(protein.protein_alphabet.monomers.A)
        form.seq.append(monomer)
        self.assertNotEqual(form.validate(), [])

    def test_incomplete_termini(self):
        monomer1 = protein.ProteinForm().from_str('{AA0062}')
        monomer2 = protein.ProteinForm().from_str('{AA0318}')

        dimer = protein.ProteinForm().from_str('{AA0062}{AA0318}')
        self.assertEqual(dimer.validate(), [])
        self.assertEqual(dimer.get_formula(),
                         monomer1.get_formula() + monomer2.get_formula() - EmpiricalFormula('OH3'))
        self.assertEqual(dimer.get_charge(),
                         monomer1.get_charge() + monomer2.get_charge() - 1)

        dimer = protein.ProteinForm().from_str('{AA0318}{AA0062}')
        self.assertNotEqual(dimer.validate(), [])
