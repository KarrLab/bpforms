""" Test of bpforms.alphabet.dna

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.alphabet import dna
from bpforms.util import validate_bpform_bonds
from wc_utils.util.chem import EmpiricalFormula
import mock
import openbabel
import os.path
import shutil
import tempfile
import unittest

damp_smiles = 'NC1=NC=NC2=C1N=CN2C1CC(O)C(COP([O-])([O-])=O)O1'
dcmp_smiles = 'NC1=NC(=O)N(C=C1)C1CC(O)C(COP([O-])([O-])=O)O1'
dgmp_smiles = 'NC1=NC2=C(N=CN2C2CC(O)C(COP([O-])([O-])=O)O2)C(=O)N1'
dtmp_smiles = 'CC1=CN(C2CC(O)C(COP([O-])([O-])=O)O2)C(=O)NC1=O'
damp_inchi = ('InChI=1S'
              '/C10H14N5O6P'
              '/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(16)6(21-7)2-20-22(17,18)19'
              '/h3-7,16H,1-2H2,(H2,11,12,13)(H2,17,18,19)'
              '/p-2')
dcmp_inchi = ('InChI=1S'
              '/C9H14N3O7P'
              '/c10-7-1-2-12(9(14)11-7)8-3-5(13)6(19-8)4-18-20(15,16)17'
              '/h1-2,5-6,8,13H,3-4H2,(H2,10,11,14)(H2,15,16,17)'
              '/p-2')
dgmp_inchi = ('InChI=1S'
              '/C10H14N5O7P'
              '/c11-10-13-8-7(9(17)14-10)12-3-15(8)6-1-4(16)5(22-6)2-21-23(18,19)20'
              '/h3-6,16H,1-2H2,(H2,18,19,20)(H3,11,13,14,17)'
              '/p-2')
dtmp_inchi = ('InChI=1S'
              '/C10H15N2O8P'
              '/c1-5-3-12(10(15)11-9(5)14)8-2-6(13)7(20-8)4-19-21(16,17)18'
              '/h3,6-8,13H,2,4H2,1H3,(H,11,14,15)(H2,16,17,18)'
              '/p-2')
di_damp_inchi = ('InChI=1S'
                 '/C20H26N10O11P2'
                 '/c21-17-15-19(25-5-23-17)29(7-27-15)13-1-9(31)11(39-13)3-38-43(35,36)'
                 '41-10-2-14(40-12(10)4-37-42(32,33)34)30-8-28-16-18(22)24-6-26-20(16)30'
                 '/h5-14,31H,1-4H2,(H,35,36)(H2,21,23,25)(H2,22,24,26)(H2,32,33,34)'
                 '/p-3')
dAdC_inchi = ('InChI=1S'
              '/C19H26N8O12P2'
              '/c20-13-1-2-26(19(29)25-13)14-3-9(28)11(37-14)5-36-41(33,34)39-10-4-15'
              '(38-12(10)6-35-40(30,31)32)27-8-24-16-17(21)22-7-23-18(16)27'
              '/h1-2,7-12,14-15,28H,3-6H2,(H,33,34)(H2,20,25,29)(H2,21,22,23)(H2,30,31,32)'
              '/p-3')
dAdCdG_inchi = ('InChI=1S'
                '/C29H38N13O18P3'
                '/c30-18-1-2-40(29(45)37-18)20-4-13(59-62(49,50)54-6-15-12(43)3-19(56-15)42-11-36-23-26'
                '(42)38-28(32)39-27(23)44)17(57-20)8-55-63(51,52)60-14-5-21(58-16(14)7-53-61(46,47)48)'
                '41-10-35-22-24(31)33-9-34-25(22)41'
                '/h1-2,9-17,19-21,43H,3-8H2,(H,49,50)(H,51,52)(H2,30,37,45)(H2,31,33,34)(H2,46,47,48)(H3,32,38,39,44)'
                '/p-4')


class DnaTestCase(unittest.TestCase):
    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_get_dnamod(self):
        filename = os.path.join(self.dirname, 'dnamod.sqlite')
        dna.get_dnamod(filename)
        self.assertTrue(os.path.isfile(filename))

        dna.get_dnamod(filename)

    def test_dna_alphabet(self):
        self.assertEqual(dna.dna_alphabet.monomers.A.get_formula(), EmpiricalFormula('C5H5N5') + EmpiricalFormula('C5H7O6P'))
        self.assertEqual(dna.dna_alphabet.monomers.C.get_formula(), EmpiricalFormula('C4H5N3O') + EmpiricalFormula('C5H7O6P'))
        self.assertEqual(dna.dna_alphabet.monomers.G.get_formula(), EmpiricalFormula('C5H5N5O') + EmpiricalFormula('C5H7O6P'))
        self.assertEqual(dna.dna_alphabet.monomers.T.get_formula(), EmpiricalFormula('C5H6N2O2') + EmpiricalFormula('C5H7O6P'))

        self.assertEqual(dna.canonical_dna_alphabet.monomers.A.get_formula(), EmpiricalFormula('C5H5N5') + EmpiricalFormula('C5H7O6P'))
        self.assertEqual(dna.canonical_dna_alphabet.monomers.C.get_formula(), EmpiricalFormula('C4H5N3O') + EmpiricalFormula('C5H7O6P'))
        self.assertEqual(dna.canonical_dna_alphabet.monomers.G.get_formula(), EmpiricalFormula('C5H5N5O') + EmpiricalFormula('C5H7O6P'))
        self.assertEqual(dna.canonical_dna_alphabet.monomers.T.get_formula(), EmpiricalFormula('C5H6N2O2') + EmpiricalFormula('C5H7O6P'))

    def test_DnaForm_init(self):
        dna.DnaForm()
        dna.CanonicalDnaForm()

    def test_DnaForm_properties(self):
        monomers = dna.canonical_dna_alphabet.monomers

        form = dna.CanonicalDnaForm().from_str('A')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C10H12N5O6P'))
        self.assertEqual(form.get_charge(), -2)
        self.assertEqual(form.export('inchi'), damp_inchi)

        form = dna.CanonicalDnaForm().from_str('C')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C9H12N3O7P'))
        self.assertEqual(form.get_charge(), -2)
        self.assertEqual(form.export('inchi'), dcmp_inchi)

        form = dna.CanonicalDnaForm().from_str('G')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C10H12N5O7P'))
        self.assertEqual(form.get_charge(), -2)
        self.assertEqual(form.export('inchi'), dgmp_inchi)

        form = dna.CanonicalDnaForm().from_str('T')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C10H13N2O8P'))
        self.assertEqual(form.get_charge(), -2)
        self.assertEqual(form.export('inchi'), dtmp_inchi)

        form = dna.CanonicalDnaForm().from_str('AA')
        self.assertEqual(form.get_formula(),
                         EmpiricalFormula('OH') * -1 + monomers.A.get_formula() * 2)
        self.assertEqual(form.get_formula(), EmpiricalFormula('C20H23N10O11P2'))
        self.assertEqual(form.get_charge(), -3)
        self.assertEqual(form.export('inchi'), di_damp_inchi)

        form = dna.CanonicalDnaForm().from_str('AC')
        self.assertEqual(form.get_formula(), EmpiricalFormula('OH') * -1
                         + monomers.A.get_formula()
                         + monomers.C.get_formula())
        self.assertEqual(form.get_formula(), EmpiricalFormula('C19H23O12N8P2'))
        self.assertEqual(form.get_charge(), -3)
        self.assertEqual(form.export('inchi'), dAdC_inchi)

        form = dna.CanonicalDnaForm().from_str('ACG')
        self.assertEqual(form.get_formula(), EmpiricalFormula('OH') * -2
                         + monomers.A.get_formula()
                         + monomers.C.get_formula()
                         + monomers.G.get_formula())
        self.assertEqual(form.get_formula(), EmpiricalFormula('C29H34O18N13P3'))
        self.assertEqual(form.get_charge(), -4)
        self.assertEqual(form.export('inchi'), dAdCdG_inchi)

        form = dna.DnaForm().from_str('ACG')
        self.assertEqual(form.get_formula(), EmpiricalFormula('OH') * -2
                         + monomers.A.get_formula()
                         + monomers.C.get_formula()
                         + monomers.G.get_formula())
        self.assertEqual(form.get_formula(), EmpiricalFormula('C29H34O18N13P3'))
        self.assertEqual(form.get_charge(), -4)
        self.assertEqual(form.export('inchi'), dAdCdG_inchi)

        form = dna.CanonicalDnaForm().from_str('ACGT')
        self.assertEqual(form.get_formula(), EmpiricalFormula('OH') * -3
                         + monomers.A.get_formula()
                         + monomers.C.get_formula()
                         + monomers.G.get_formula()
                         + monomers.T.get_formula())
        self.assertEqual(form.get_formula(), EmpiricalFormula('C39H46O25N15P4'))
        self.assertEqual(form.get_charge(), -5)

        form = dna.DnaForm().from_str('ACGT')
        self.assertEqual(form.get_formula(), EmpiricalFormula('OH') * -3
                         + monomers.A.get_formula()
                         + monomers.C.get_formula()
                         + monomers.G.get_formula()
                         + monomers.T.get_formula())
        self.assertEqual(form.get_formula(), EmpiricalFormula('C39H46O25N15P4'))
        self.assertEqual(form.get_charge(), -5)

    def test_DnaAlphabetBuilder(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        alphabet = dna.DnaAlphabetBuilder(_max_monomers=3).run(ph=7.4, path=path)
        alphabet = dna.DnaAlphabetBuilder().run(ph=7.4, path=path)
        self.assertEqual(alphabet.monomers.A.get_formula(),
                         EmpiricalFormula('C5H5N5') + EmpiricalFormula('C5H7O6P'))
        self.assertTrue(os.path.isfile(path))

    def test_validate_bonds(self):
        validate_bpform_bonds(dna.CanonicalDnaForm)
        validate_bpform_bonds(dna.DnaForm)

    def test_validate_form(self):
        form = dna.DnaForm()
        form.from_str('ACGT')
        self.assertEqual(form.validate(), [])
