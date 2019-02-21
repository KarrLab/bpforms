""" Test of bpforms.alphabet.rna

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Identifier, IdentifierSet, Monomer
from bpforms.alphabet import rna
from wc_utils.util.chem import EmpiricalFormula
import mock
import os.path
import requests
import shutil
import tempfile
import unittest


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

        form = rna.CanonicalRnaForm().from_str('ACG')
        self.assertEqual(form.get_formula(), EmpiricalFormula('C29H34O21N13P3'))
        self.assertEqual(form.get_charge(), -4)

        form = rna.CanonicalRnaForm().from_str('ACGU')
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
                'InChI=1S/C12H17N5O4'
                '/c1-16-4-15-11-7(10(16)13)14-5-17(11)12-9(20-2)8(19)6(3-18)21-12'
                '/h4-6,8-9,12-13,18-19H,3H2,1-2H3'))),
            True)
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside('100G (base)', None, None, None), False)
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside('79553N', None, '7-methylguanosine cap (cap 0)', None), False)
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside('2551N', None, 'alpha-dimethylmonophosphate cap', None), False)
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside('455N', 'CoA(pN)', '5′ (3′ -dephospho-CoA)', None), False)
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside('00A', 'Ar(p)', '2′-O-ribosyladenosine (phosphate)', None), False)
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside('00A', '5′-OH-N', '5′ hydroxyl end', Monomer()), False)
        self.assertEqual(rna.RnaAlphabetBuilder().is_valid_nucleoside(
            '00A', '5′-OH-N', '5′ hydroxyl end', Monomer(structure='InChI=1S/H3O3P/c1-4(2)3/h4H,(H2,1,2,3)/p-2')), False)

    def test_RnaAlphabetBuilder(self):
        path = os.path.join(self.dirname, 'alphabet.yml')
        alphabet = rna.RnaAlphabetBuilder(_max_monomers=3).run(path=path)
        alphabet = rna.RnaAlphabetBuilder().run(path=path)
        self.assertEqual(alphabet.monomers.A.get_formula(), EmpiricalFormula('C10O4N5H13'))
        self.assertTrue(os.path.isfile(path))
