""" Test of bpforms.alphabet.rna

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Identifier, IdentifierSet, Monomer
from bpforms.alphabet import rna
from bpforms.util import validate_bpform_linkages
from wc_utils.util.chem import EmpiricalFormula
import mock
import openbabel
import os.path
import requests
import shutil
import tempfile
import unittest

amp_smiles = 'C(C3(C(C(C(N2(C1(=C(C(=NC=N1)N)N=C2)))O3)O)O))OP([O-])([O-])=O'
cmp_smiles = 'C(C2(C(C(C(N1(C(N=C(C=C1)N)=O))O2)O)O))OP([O-])([O-])=O'
gmp_smiles = 'C(OP(=O)([O-])[O-])C1(OC(C(O)C(O)1)N3(C=NC2(C(=O)NC(N)=NC=23)))'
ump_smiles = 'C(OP(=O)([O-])[O-])C1(OC(C(O)C(O)1)N2(C=CC(=O)NC(=O)2))'
amp_inchi = ('InChI=1S/C10H14N5O7P/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(22-10)1-21-23(18,19)20'
             '/h2-4,6-7,10,16-17H,1H2,(H2,11,12,13)(H2,18,19,20)/p-2')
cmp_inchi = ('InChI=1S/C9H14N3O8P/c10-5-1-2-12(9(15)11-5)8-7(14)6(13)4(20-8)3-19-21(16,17)18'
             '/h1-2,4,6-8,13-14H,3H2,(H2,10,11,15)(H2,16,17,18)/p-2')
gmp_inchi = ('InChI=1S/C10H14N5O8P/c11-10-13-7-4(8(18)14-10)12-2-15(7)9-6(17)5(16)3(23-9)1-22-24(19,20)21'
             '/h2-3,5-6,9,16-17H,1H2,(H2,19,20,21)(H3,11,13,14,18)/p-2')
ump_inchi = ('InChI=1S/C9H13N2O9P/c12-5-1-2-11(9(15)10-5)8-7(14)6(13)4(20-8)3-19-21(16,17)18'
             '/h1-2,4,6-8,13-14H,3H2,(H,10,12,15)(H2,16,17,18)/p-2')
di_amp_inchi = ('InChI=1S'
                '/C20H26N10O13P2'
                '/c21-15-9-17(25-3-23-15)29(5-27-9)19-12(32)11(31)7(41-19)1-40-45(37,38)43-14-8(2-39-44(34,35)36)'
                '42-20(13(14)33)30-6-28-10-16(22)24-4-26-18(10)30'
                '/h3-8,11-14,19-20,31-33H,1-2H2,(H,37,38)(H2,21,23,25)(H2,22,24,26)(H2,34,35,36)'
                '/p-3')
AC_inchi = ('InChI=1S'
            '/C19H26N8O14P2'
            '/c20-9-1-2-26(19(31)25-9)17-12(29)11(28)7(39-17)3-38-43(35,36)41-14-8(4-37-42(32,33)34)40-18(13(14)30)'
            '27-6-24-10-15(21)22-5-23-16(10)27'
            '/h1-2,5-8,11-14,17-18,28-30H,3-4H2,(H,35,36)(H2,20,25,31)(H2,21,22,23)(H2,32,33,34)'
            '/p-3')
ACG_inchi = ('InChI=1S'
             '/C29H38N13O21P3'
             '/c30-12-1-2-40(29(48)37-12)26-17(45)20(63-65(52,53)57-3-9-15(43)16(44)25(59-9)42-8-36-14-23(42)38-28(32)39-24'
             '(14)47)11(60-26)5-58-66(54,55)62-19-10(4-56-64(49,50)51)61-27(18(19)46)41-7-35-13-21(31)33-6-34-22(13)41'
             '/h1-2,6-11,15-20,25-27,43-46H,3-5H2,(H,52,53)(H,54,55)(H2,30,37,48)(H2,31,33,34)(H2,49,50,51)(H3,32,38,39,47)'
             '/p-4')


class RnaTestCase(unittest.TestCase):
    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_rna_alphabet(self):
        self.assertEqual(rna.rna_alphabet.monomers.A.get_formula(), EmpiricalFormula('C10O4N5H13'))
        self.assertEqual(rna.rna_alphabet.monomers.C.get_formula(), EmpiricalFormula('C9O5N3H13'))
        self.assertEqual(rna.rna_alphabet.monomers.G.get_formula(), EmpiricalFormula('C10O5N5H13'))
        self.assertEqual(rna.rna_alphabet.monomers.U.get_formula(), EmpiricalFormula('C9O6N2H12'))

        self.assertEqual(rna.canonical_rna_alphabet.monomers.A.get_formula(), EmpiricalFormula('C10O4N5H13'))
        self.assertEqual(rna.canonical_rna_alphabet.monomers.C.get_formula(), EmpiricalFormula('C9O5N3H13'))
        self.assertEqual(rna.canonical_rna_alphabet.monomers.G.get_formula(), EmpiricalFormula('C10O5N5H13'))
        self.assertEqual(rna.canonical_rna_alphabet.monomers.U.get_formula(), EmpiricalFormula('C9O6N2H12'))

    def test_RnaForm_init(self):
        rna.RnaForm()
        rna.CanonicalRnaForm()

    def test_RnaForm_properties(self):
        monomers = rna.canonical_rna_alphabet.monomers

        form = rna.CanonicalRnaForm().from_str('A')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C10H12N5O7P'))
        self.assertEqual(form.get_charge(), -2)
        self.assertEqual(form.export('inchi'), amp_inchi)

        form = rna.CanonicalRnaForm().from_str('C')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C9H12N3O8P'))
        self.assertEqual(form.get_charge(), -2)
        self.assertEqual(form.export('inchi'), cmp_inchi)

        form = rna.CanonicalRnaForm().from_str('G')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C10H12N5O8P'))
        self.assertEqual(form.get_charge(), -2)
        self.assertEqual(form.export('inchi'), gmp_inchi)

        form = rna.CanonicalRnaForm().from_str('U')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C9H11N2O9P'))
        self.assertEqual(form.get_charge(), -2)
        self.assertEqual(form.export('inchi'), ump_inchi)

        form = rna.CanonicalRnaForm().from_str('AA')
        self.assertEqual(form.get_formula(), EmpiricalFormula('H-1O3P') * 2 + EmpiricalFormula('OH') * -1
                         + monomers.A.get_formula() * 2)
        self.assertEqual(form.get_formula(), EmpiricalFormula('C20H23N10O13P2'))
        self.assertEqual(form.get_charge(), -3)
        print(form.export('inchi'))
        self.assertEqual(form.export('inchi'), di_amp_inchi)

        form = rna.CanonicalRnaForm().from_str('AC')
        self.assertEqual(form.get_formula(), EmpiricalFormula('H-1O3P') * 2 + EmpiricalFormula('OH') * -1
                         + monomers.A.get_formula()
                         + monomers.C.get_formula())
        self.assertEqual(form.get_formula(), EmpiricalFormula('C19H23O14N8P2'))
        self.assertEqual(form.get_charge(), -3)
        self.assertEqual(form.export('inchi'), AC_inchi)

        form = rna.CanonicalRnaForm().from_str('ACG')
        self.assertEqual(form.get_formula(), EmpiricalFormula('H-1O3P') * 3 + EmpiricalFormula('OH') * -2
                         + monomers.A.get_formula()
                         + monomers.C.get_formula()
                         + monomers.G.get_formula())
        self.assertEqual(form.get_formula(), EmpiricalFormula('C29H34O21N13P3'))
        self.assertEqual(form.get_charge(), -4)
        self.assertEqual(form.export('inchi'), ACG_inchi)

        form = rna.RnaForm().from_str('ACG')
        self.assertEqual(form.get_formula(), EmpiricalFormula('H-1O3P') * 3 + EmpiricalFormula('OH') * -2
                         + monomers.A.get_formula()
                         + monomers.C.get_formula()
                         + monomers.G.get_formula())
        self.assertEqual(form.get_formula(), EmpiricalFormula('C29H34O21N13P3'))
        self.assertEqual(form.get_charge(), -4)
        self.assertEqual(form.export('inchi'), ACG_inchi)

        form = rna.CanonicalRnaForm().from_str('ACGU')
        self.assertEqual(form.get_formula(), EmpiricalFormula('H-1O3P') * 4 + EmpiricalFormula('OH') * -3
                         + monomers.A.get_formula()
                         + monomers.C.get_formula()
                         + monomers.G.get_formula()
                         + monomers.U.get_formula())
        self.assertEqual(form.get_formula(), EmpiricalFormula('C38H44O29N15P4'))
        self.assertEqual(form.get_charge(), -5)

        form = rna.RnaForm().from_str('ACGU')
        self.assertEqual(form.get_formula(), EmpiricalFormula('H-1O3P') * 4 + EmpiricalFormula('OH') * -3
                         + monomers.A.get_formula()
                         + monomers.C.get_formula()
                         + monomers.G.get_formula()
                         + monomers.U.get_formula())
        self.assertEqual(form.get_formula(), EmpiricalFormula('C38H44O29N15P4'))
        self.assertEqual(form.get_charge(), -5)

    def test_RnaAlphabetBuilder_list_monomers(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        with mock.patch.object(rna.RnaAlphabetBuilder, 'get_nucleoside_details', return_value=(None, IdentifierSet())):
            alphabet = rna.RnaAlphabetBuilder().run(path=path)
        self.assertEqual(alphabet.monomers.A.get_formula(), EmpiricalFormula('C10O4N5H13'))
        self.assertTrue(os.path.isfile(path))

    def test_RnaAlphabetBuilder_get_nucleoside_details(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        session = requests.Session()
        self.assertEqual(rna.RnaAlphabetBuilder().get_nucleoside_details('A', session)[1],
                         IdentifierSet([Identifier('pubchem.compound', '60961')]))
        self.assertEqual(rna.RnaAlphabetBuilder().get_nucleoside_details('(pN)2′3′>p', session), (None, IdentifierSet()))
        self.assertEqual(rna.RnaAlphabetBuilder().get_nucleoside_details('Xm', session), (None, IdentifierSet()))

    def test_RnaAlphabetBuilder_is_valid_nucleoside(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        session = requests.Session()
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside(
            '01A', 'm1Am', '1,2′-O-dimethyladenosine', Monomer(structure=(
                'C3(OC(CO)C(O)C(O(C))3)N2C=NC1=C2N=CN(C)C1=N'))),
            True)
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside('100G (base)', None, None, None), False)
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside('79553N', None, '7-methylguanosine cap (cap 0)', None), False)
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside('2551N', None, 'alpha-dimethylmonophosphate cap', None), False)
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside('455N', 'CoA(pN)', '5′ (3′ -dephospho-CoA)', None), False)
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside('00A', 'Ar(p)', '2′-O-ribosyladenosine (phosphate)', None), False)
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside('00A', '5′-OH-N', '5′ hydroxyl end', Monomer()), False)
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside(
            '00A', '5′-OH-N', '5′ hydroxyl end', Monomer(structure='[O-]P([O-])=O')), False)

    def test_RnaAlphabetBuilder(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        alphabet = rna.RnaAlphabetBuilder(_max_monomers=3).run(path=path)
        alphabet = rna.RnaAlphabetBuilder().run(path=path)
        self.assertEqual(alphabet.monomers.A.get_formula(), EmpiricalFormula('C10O4N5H13'))
        self.assertTrue(os.path.isfile(path))

    def test_validate_linkages(self):
        validate_bpform_linkages(rna.CanonicalRnaForm)
        validate_bpform_linkages(rna.RnaForm)

        for monomer in rna.CanonicalRnaForm().alphabet.monomers.values():
            b_atom = monomer.structure.GetAtom(monomer.monomer_bond_atoms[0].position)
            r_atom = monomer.structure.GetAtom(monomer.left_bond_atoms[0].position)
            assert self.are_atoms_valid(b_atom, r_atom)

        for monomer in rna.RnaForm().alphabet.monomers.values():
            b_atom = monomer.structure.GetAtom(monomer.monomer_bond_atoms[0].position)
            r_atom = monomer.structure.GetAtom(monomer.left_bond_atoms[0].position)
            assert self.are_atoms_valid(b_atom, r_atom)

    def are_atoms_valid(self, b_atom, r_atom):
        assert b_atom.GetAtomicNum() == 8
        assert b_atom.GetFormalCharge() == 0

        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(b_atom)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(b_atom)])
        other_atoms += [1] * (2 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        assert other_atoms == [1, 6]

        # get first C
        for other_atom in openbabel.OBAtomAtomIter(b_atom):
            if other_atom.GetAtomicNum() == 6:
                c_1 = other_atom
        assert c_1.GetFormalCharge() == 0
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(c_1)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(c_1)])
        other_atoms += [1] * (4 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        assert other_atoms == [1, 1, 6, 8]

        # get second C
        for other_atom in openbabel.OBAtomAtomIter(c_1):
            if other_atom.GetAtomicNum() == 6:
                c_2 = other_atom
        assert c_2.GetFormalCharge() == 0
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(c_2)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(c_2)])
        other_atoms += [1] * (4 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        assert other_atoms == [1, 6, 6, 8]

        # get third C
        for other_atom in openbabel.OBAtomAtomIter(c_2):
            if other_atom.GetAtomicNum() == 6 and other_atom.GetIdx() != c_1.GetIdx():
                c_3 = other_atom
        assert c_3.GetFormalCharge() == 0
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(c_3)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(c_3)])
        other_atoms += [1] * (4 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        assert other_atoms == [1, 6, 6, 8]

        # get second O
        for other_atom in openbabel.OBAtomAtomIter(c_3):
            if other_atom.GetAtomicNum() == 8:
                r_atom_2 = other_atom
        assert r_atom_2.GetFormalCharge() == 0
        assert r_atom_2.GetIdx() == r_atom.GetIdx()
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(r_atom_2)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(r_atom_2)])
        other_atoms += [1] * (2 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        assert other_atoms == [1, 6]

        return True
