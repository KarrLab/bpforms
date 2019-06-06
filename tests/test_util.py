""" Test of bpforms.util

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.alphabet import dna
from bpforms.alphabet import rna
from bpforms.alphabet import protein
from bpforms import core
from bpforms import util
from wc_utils.util.chem import EmpiricalFormula
import mock
import os
import shutil
import tempfile
import unittest


class UtilTestCase(unittest.TestCase):
    def setUp(self):
        os.rename(dna.filename, dna.filename + '.save')
        os.rename(rna.filename, rna.filename + '.save')
        os.rename(protein.filename, protein.filename + '.save')
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        os.rename(dna.filename + '.save', dna.filename)
        os.rename(rna.filename + '.save', rna.filename)
        os.rename(protein.filename + '.save', protein.filename)
        shutil.rmtree(self.tempdir)

    def test_get_alphabets(self):
        self.assertEqual(util.get_alphabets()['dna'], dna.dna_alphabet)

    def test_get_alphabet(self):
        self.assertEqual(util.get_alphabet('dna'), dna.dna_alphabet)
        self.assertEqual(util.get_alphabet('rna'), rna.rna_alphabet)
        self.assertEqual(util.get_alphabet('protein'), protein.protein_alphabet)
        with self.assertRaises(ValueError):
            util.get_alphabet('lipid')

    def test_get_form(self):
        self.assertEqual(util.get_form('dna'), dna.DnaForm)
        self.assertEqual(util.get_form('rna'), rna.RnaForm)
        self.assertEqual(util.get_form('protein'), protein.ProteinForm)
        self.assertEqual(util.get_form('canonical_dna'), dna.CanonicalDnaForm)
        self.assertEqual(util.get_form('canonical_rna'), rna.CanonicalRnaForm)
        self.assertEqual(util.get_form('canonical_protein'), protein.CanonicalProteinForm)
        with self.assertRaises(ValueError):
            util.get_form('lipid')

    def test_build_alphabets(self):
        self.assertFalse(os.path.isfile(dna.filename))
        self.assertFalse(os.path.isfile(rna.filename))
        self.assertFalse(os.path.isfile(protein.filename))

        util.build_alphabets(_max_monomers=3)

        self.assertTrue(os.path.isfile(dna.filename))
        self.assertTrue(os.path.isfile(rna.filename))
        self.assertTrue(os.path.isfile(protein.filename))

    def test_gen_html_viz_alphabet(self):
        filename = os.path.join(self.tempdir, 'alphabet.html')
        util.gen_html_viz_alphabet(dna.CanonicalDnaForm, filename)
        self.assertTrue(os.path.isfile(filename))

        alphabet = core.Alphabet(monomers={'A': core.Monomer()})

        class TestForm(core.BpForm):
            def __init__(self, seq=None, circular=False):
                super(TestForm, self).__init__(
                    seq=seq, alphabet=alphabet,
                    backbone=core.Backbone(
                        structure='OC1CCOC1COP([O-])([O-])=O',
                        monomer_bond_atoms=[core.Atom(core.Backbone, element='C', position=4)],
                        monomer_displaced_atoms=[core.Atom(core.Backbone, element='H', position=4)]),
                    bond=core.Bond(
                        left_bond_atoms=[core.Atom(core.Backbone, element='O', position=1)],
                        right_bond_atoms=[core.Atom(core.Backbone, element='P', position=9)],
                        left_displaced_atoms=[core.Atom(core.Backbone, element='H', position=1)],
                        right_displaced_atoms=[core.Atom(core.Backbone, element='O', position=11, charge=-1)]),
                    circular=circular)
        util.gen_html_viz_alphabet(TestForm, filename)

    def test_validate_bpform_linkages(self):
        class TestDnaForm(core.BpForm):
            def __init__(self, seq=None):
                super(TestDnaForm, self).__init__(
                    seq=seq, alphabet=dna.canonical_dna_alphabet,
                    backbone=core.Backbone(
                        structure='OC1CCOC1COP([O-])([O-])=O',
                        monomer_bond_atoms=[core.Atom(core.Backbone, element='C', position=4)],
                        monomer_displaced_atoms=[core.Atom(core.Backbone, element='H', position=4)]),
                    bond=core.Bond(
                        left_bond_atoms=[core.Atom(core.Backbone, element='O', position=1)],
                        right_bond_atoms=[core.Atom(core.Backbone, element='P', position=9)],
                        left_displaced_atoms=[core.Atom(core.Backbone, element='H', position=1)],
                        right_displaced_atoms=[core.Atom(core.Backbone, element='O', position=11, charge=-1)]))
        util.validate_bpform_linkages(TestDnaForm)

        class TestDnaForm(core.BpForm):
            def __init__(self, seq=None):
                super(TestDnaForm, self).__init__(
                    seq=seq, alphabet=dna.canonical_dna_alphabet,
                    backbone=core.Backbone(
                        structure='OC1CCOC1COP([O-])([O-])=O',
                        monomer_bond_atoms=[core.Atom(core.Backbone, element='C', position=1000)],
                        monomer_displaced_atoms=[core.Atom(core.Backbone, element='H', position=4)]),
                    bond=core.Bond(
                        left_bond_atoms=[core.Atom(core.Backbone, element='O', position=1)],
                        right_bond_atoms=[core.Atom(core.Backbone, element='P', position=9)],
                        left_displaced_atoms=[core.Atom(core.Backbone, element='H', position=1)],
                        right_displaced_atoms=[core.Atom(core.Backbone, element='O', position=11, charge=-1)]))
        with self.assertRaises(ValueError):
            util.validate_bpform_linkages(TestDnaForm)

        class TestDnaForm(core.BpForm):
            def __init__(self, seq=None):
                super(TestDnaForm, self).__init__(
                    seq=seq, alphabet=dna.canonical_dna_alphabet,
                    backbone=core.Backbone(
                        structure='OC1CCOC1COP([O-])([O-])=O',
                        monomer_bond_atoms=[core.Atom(core.Backbone, element='E', position=4)],
                        monomer_displaced_atoms=[core.Atom(core.Backbone, element='H', position=4)]),
                    bond=core.Bond(
                        left_bond_atoms=[core.Atom(core.Backbone, element='O', position=1)],
                        right_bond_atoms=[core.Atom(core.Backbone, element='P', position=9)],
                        left_displaced_atoms=[core.Atom(core.Backbone, element='H', position=1)],
                        right_displaced_atoms=[core.Atom(core.Backbone, element='O', position=11, charge=-1)]))
        with self.assertRaises(ValueError):
            util.validate_bpform_linkages(TestDnaForm)

        alphabet = core.Alphabet()
        alphabet.monomers.A = core.Monomer(id='adenine', structure='NC1=C2N=CNC2=NC=N1',
                                           backbone_bond_atoms=[core.Atom(molecule=core.Monomer, element='N', position=6)],
                                           backbone_displaced_atoms=[core.Atom(molecule=core.Monomer, element='H', position=6)],
                                           )

        class TestDnaForm(core.BpForm):
            def __init__(self, seq=None):
                super(TestDnaForm, self).__init__(
                    seq=seq, alphabet=alphabet,
                    backbone=core.Backbone(
                        structure='OC1CCOC1COP([O-])([O-])=O',
                        monomer_bond_atoms=[core.Atom(core.Backbone, element='C', position=4)],
                        monomer_displaced_atoms=[core.Atom(core.Backbone, element='H', position=4)]),
                    bond=core.Bond(
                        left_bond_atoms=[core.Atom(core.Backbone, element='O', position=1)],
                        right_bond_atoms=[core.Atom(core.Backbone, element='P', position=9)],
                        left_displaced_atoms=[core.Atom(core.Backbone, element='H', position=1)],
                        right_displaced_atoms=[core.Atom(core.Backbone, element='O', position=11, charge=-1)]))
        util.validate_bpform_linkages(TestDnaForm)

        alphabet = core.Alphabet()
        alphabet.monomers.A = core.Monomer(id='adenine', structure='NC1=C2N=CNC2=NC=N1',
                                           backbone_bond_atoms=[core.Atom(molecule=core.Monomer, element='N', position=1000)],
                                           backbone_displaced_atoms=[core.Atom(molecule=core.Monomer, element='H', position=6)],
                                           )

        class TestDnaForm(core.BpForm):
            def __init__(self, seq=None):
                super(TestDnaForm, self).__init__(
                    seq=seq, alphabet=alphabet,
                    backbone=core.Backbone(
                        structure='OC1CCOC1COP([O-])([O-])=O',
                        monomer_bond_atoms=[core.Atom(core.Backbone, element='C', position=4)],
                        monomer_displaced_atoms=[core.Atom(core.Backbone, element='H', position=4)]),
                    bond=core.Bond(
                        left_bond_atoms=[core.Atom(core.Backbone, element='O', position=1)],
                        right_bond_atoms=[core.Atom(core.Backbone, element='P', position=9)],
                        left_displaced_atoms=[core.Atom(core.Backbone, element='H', position=1)],
                        right_displaced_atoms=[core.Atom(core.Backbone, element='O', position=11, charge=-1)]))
        with self.assertRaises(ValueError):
            util.validate_bpform_linkages(TestDnaForm)

        alphabet = core.Alphabet()
        alphabet.monomers.A = core.Monomer(id='adenine', structure='NC1=C2N=CNC2=NC=N1',
                                           backbone_bond_atoms=[core.Atom(molecule=core.Monomer, element='E', position=6)],
                                           backbone_displaced_atoms=[core.Atom(molecule=core.Monomer, element='H', position=6)],
                                           )

        class TestDnaForm(core.BpForm):
            def __init__(self, seq=None):
                super(TestDnaForm, self).__init__(
                    seq=seq, alphabet=alphabet,
                    backbone=core.Backbone(
                        structure='OC1CCOC1COP([O-])([O-])=O',
                        monomer_bond_atoms=[core.Atom(core.Backbone, element='C', position=4)],
                        monomer_displaced_atoms=[core.Atom(core.Backbone, element='H', position=4)]),
                    bond=core.Bond(
                        left_bond_atoms=[core.Atom(core.Backbone, element='O', position=1)],
                        right_bond_atoms=[core.Atom(core.Backbone, element='P', position=9)],
                        left_displaced_atoms=[core.Atom(core.Backbone, element='H', position=1)],
                        right_displaced_atoms=[core.Atom(core.Backbone, element='O', position=11, charge=-1)]))
        with self.assertRaises(ValueError):
            util.validate_bpform_linkages(TestDnaForm)

        alphabet = core.Alphabet()
        alphabet.monomers.A = core.Monomer(id='adenine', structure='NC1=C2N=CNC2=NC=N1',
                                           backbone_bond_atoms=[core.Atom(molecule=core.Monomer, element='N', position=6)],
                                           backbone_displaced_atoms=[core.Atom(molecule=core.Monomer, element='H', position=6)],
                                           )

        class TestDnaForm(core.BpForm):
            def __init__(self, seq=None):
                super(TestDnaForm, self).__init__(
                    seq=seq, alphabet=alphabet,
                    backbone=core.Backbone(
                        structure='OC1CCOC1COP([O-])([O-])=O',
                        monomer_bond_atoms=[core.Atom(core.Backbone, element='C', position=4)],
                        monomer_displaced_atoms=[core.Atom(core.Backbone, element='H', position=4)]),
                    bond=core.Bond(
                        left_bond_atoms=[core.Atom(core.Backbone, element='O', position=1)],
                        right_bond_atoms=[core.Atom(core.Backbone, element='P', position=9)],
                        left_displaced_atoms=[core.Atom(core.Backbone, element='H', position=1)],
                        right_displaced_atoms=[core.Atom(core.Backbone, element='O', position=11, charge=-1)]))
        with self.assertRaises(ValueError):
            with mock.patch('wc_utils.util.chem.OpenBabelUtils.export', side_effect=['smiles1', 'inchi1', 'inchi2']):
                util.validate_bpform_linkages(TestDnaForm)

        with self.assertRaises(ValueError):
            with mock.patch.object(core.BpForm, 'get_formula', return_value=EmpiricalFormula()):
                util.validate_bpform_linkages(TestDnaForm)

        with self.assertRaises(ValueError):
            with mock.patch.object(core.BpForm, 'get_formula', side_effect=[EmpiricalFormula('C10H12N5O6P'), EmpiricalFormula()]):
                util.validate_bpform_linkages(TestDnaForm)

        with self.assertRaises(ValueError):
            with mock.patch.object(core.BpForm, 'get_charge', return_value=-100):
                util.validate_bpform_linkages(TestDnaForm)

        with self.assertRaises(ValueError):
            with mock.patch.object(core.BpForm, 'get_charge', side_effect=[-2, -100]):
                util.validate_bpform_linkages(TestDnaForm)
