""" Test of bpforms.core

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms import core
from bpforms import dna
from bpforms import rna
from wc_utils.util.chem import EmpiricalFormula
import copy
import lark.exceptions
import mock
import openbabel
import os
import shutil
import tempfile
import unittest

dAMP_inchi = dna.dna_alphabet.A.get_inchi()
dCMP_inchi = dna.dna_alphabet.C.get_inchi()
dGMP_inchi = dna.dna_alphabet.G.get_inchi()
dIMP_inchi = 'InChI=1S/C10H12N4O4/c15-2-6-5(16)1-7(18-6)14-4-13-8-9(14)11-3-12-10(8)17/h3-7,15-16H,1-2H2,(H,11,12,17)/t5-,6+,7+/m0/s1'


class IdentifierTestCase(unittest.TestCase):
    def test_constructor(self):
        id = core.Identifier('ec-code', '1.1.1.2')
        self.assertEqual(id.ns, 'ec-code')
        self.assertEqual(id.id, '1.1.1.2')

        id = core.Identifier('kegg.compound', 'C00001')
        self.assertEqual(id.ns, 'kegg.compound')
        self.assertEqual(id.id, 'C00001')

        id = core.Identifier('chebi', 'CHEBI:57566')
        self.assertEqual(id.ns, 'chebi')
        self.assertEqual(id.id, 'CHEBI:57566')

        id = core.Identifier('metacyc.compound', 'DAMP')
        self.assertEqual(id.ns, 'metacyc.compound')
        self.assertEqual(id.id, 'DAMP')

        id = core.Identifier('pubchem.compound', '22848660')
        self.assertEqual(id.ns, 'pubchem.compound')
        self.assertEqual(id.id, '22848660')

        with self.assertRaises(ValueError):
            id.ns = ''
        with self.assertRaises(ValueError):
            id.ns = '"pubchem"'
        with self.assertRaises(ValueError):
            id.ns = 0
        with self.assertRaises(ValueError):
            id.ns = None
        with self.assertRaises(ValueError):
            id.id = ''
        with self.assertRaises(ValueError):
            id.id = 0
        with self.assertRaises(ValueError):
            id.id = None
        with self.assertRaises(ValueError):
            id.id = '"22848660"'
        self.assertEqual(id.ns, 'pubchem.compound')
        self.assertEqual(id.id, '22848660')

    def test_eq(self):
        id_1 = core.Identifier('ec-code', '1.1.1.2')
        id_2 = core.Identifier('ec-code', '1.1.1.2')
        id_3 = core.Identifier('ec-code', '1.1.1.1')
        id_4 = core.Identifier('kegg.compound', '1.1.1.2')

        self.assertEqual(id_1, id_2)
        self.assertEqual(id_2, id_1)
        self.assertNotEqual(id_1, id_3)
        self.assertNotEqual(id_1, id_4)
        self.assertNotEqual(id_1, '')
        self.assertNotEqual(id_1, 0)

    def test_hash(self):
        id_1 = core.Identifier('ec-code', '1.1.1.2')
        id_2 = core.Identifier('ec-code', '1.1.1.2')
        id_3 = core.Identifier('ec-code', '1.1.1.1')
        id_4 = core.Identifier('kegg.compound', '1.1.1.2')

        self.assertEqual(id_1.__hash__(), id_2.__hash__())
        self.assertNotEqual(id_1.__hash__(), id_3.__hash__())
        self.assertNotEqual(id_1.__hash__(), id_4.__hash__())

    def test_set(self):
        id_1 = core.Identifier('ec-code', '1.1.1.2')
        id_2 = core.Identifier('ec-code', '1.1.1.2')
        id_3 = core.Identifier('ec-code', '1.1.1.1')
        id_4 = core.Identifier('kegg.compound', '1.1.1.2')

        ids = set([id_1, id_2, id_3, id_4])
        self.assertEqual(len(ids), 3)


class IdentifiersTestCase(unittest.TestCase):

    def test_constructor(self):
        id_1 = core.Identifier('ec-code', '1.1.1.2')
        id_2 = core.Identifier('ec-code', '1.1.1.1')
        id_3 = core.Identifier('kegg.compound', '1.1.1.2')

        ids = core.IdentifierSet()
        self.assertEqual(len(ids), 0)

        ids = core.IdentifierSet([id_1, id_2, id_3])
        self.assertEqual(len(ids), 3)

        with self.assertRaises(ValueError):
            core.IdentifierSet(['ec-code'])

    def test_add(self):
        id_1 = core.Identifier('ec-code', '1.1.1.2')
        id_2 = core.Identifier('ec-code', '1.1.1.1')
        id_3 = core.Identifier('kegg.compound', '1.1.1.2')

        ids = core.IdentifierSet()
        ids.add(id_1)
        ids.add(id_2)
        ids.add(id_3)
        self.assertEqual(len(ids), 3)

        with self.assertRaises(ValueError):
            ids.add(('ec-code', '1.1.1.2'))

    def test_update(self):
        id_1 = core.Identifier('ec-code', '1.1.1.2')
        id_2 = core.Identifier('ec-code', '1.1.1.1')
        id_3 = core.Identifier('kegg.compound', '1.1.1.2')

        ids = core.IdentifierSet()
        ids.update([id_1, id_2])
        self.assertEqual(len(ids), 2)

        with self.assertRaises(ValueError):
            ids.update([('ec-code', '1.1.1.2')])

    def test_symmetric_difference_update(self):
        id_1 = core.Identifier('ec-code', '1.1.1.2')
        id_2 = core.Identifier('ec-code', '1.1.1.1')
        id_3 = core.Identifier('kegg.compound', '1.1.1.2')

        ids_1 = core.IdentifierSet([id_1, id_2])
        ids_2 = core.IdentifierSet([id_1, id_3])
        ids_1.symmetric_difference_update(ids_2)
        self.assertEqual(ids_1, core.IdentifierSet([id_2, id_3]))

        ids_1 = core.IdentifierSet([id_1, id_2])
        ids_2 = set([id_1, id_3])
        ids_1.symmetric_difference_update(ids_2)
        self.assertEqual(ids_1, core.IdentifierSet([id_2, id_3]))

        with self.assertRaises(TypeError):
            ids_1.symmetric_difference_update(id_1)


class SynonymsTestCase(unittest.TestCase):

    def test_constructor(self):
        syn_1 = 'a'
        syn_2 = 'b'
        syn_3 = 'c'

        syns = core.SynonymSet()
        self.assertEqual(len(syns), 0)

        syns = core.SynonymSet([syn_1, syn_2, syn_3])
        self.assertEqual(len(syns), 3)

        with self.assertRaises(ValueError):
            core.SynonymSet('')
        with self.assertRaises(ValueError):
            core.SynonymSet('a')
        with self.assertRaises(ValueError):
            core.SynonymSet([''])
        with self.assertRaises(ValueError):
            core.SynonymSet([0])
        with self.assertRaises(ValueError):
            core.SynonymSet([None])

    def test_add(self):
        syn_1 = 'a'
        syn_2 = 'b'
        syn_3 = 'c'

        syns = core.SynonymSet()
        syns.add(syn_1)
        syns.add(syn_2)
        syns.add(syn_3)
        self.assertEqual(len(syns), 3)

        with self.assertRaises(ValueError):
            syns.add(0)
        with self.assertRaises(ValueError):
            syns.add(None)
        with self.assertRaises(ValueError):
            syns.add('')

    def test_update(self):
        syn_1 = 'a'
        syn_2 = 'b'
        syn_3 = 'c'

        syns = core.SynonymSet()
        syns.update([syn_1, syn_2])
        self.assertEqual(len(syns), 2)

        with self.assertRaises(ValueError):
            syns.update([0])

    def test_symmetric_difference_update(self):
        syn_1 = 'a'
        syn_2 = 'b'
        syn_3 = 'c'

        syns_1 = core.SynonymSet([syn_1, syn_2])
        syns_2 = core.SynonymSet([syn_1, syn_3])
        syns_1.symmetric_difference_update(syns_2)
        self.assertEqual(syns_1, core.SynonymSet([syn_2, syn_3]))

        syns_1 = core.SynonymSet([syn_1, syn_2])
        syns_2 = set([syn_1, syn_3])
        syns_1.symmetric_difference_update(syns_2)
        self.assertEqual(syns_1, core.SynonymSet([syn_2, syn_3]))

        with self.assertRaises(TypeError):
            syns_1.symmetric_difference_update(0)


class BaseTestCase(unittest.TestCase):
    def test_init(self):
        identifiers = set([
            core.Identifier('chebi', 'CHEBI:58245'),
            core.Identifier('pubchem.compound', '22848660'),
            core.Identifier('metacyc.compound', 'DAMP'),
        ])
        synonyms = set(['A', 'dAMP', 'deoxyadenosine monophosphate'])
        base = core.Base(id='dAMP', name='deoxyadenosine monophosphate', synonyms=synonyms, identifiers=identifiers,
                         structure=dAMP_inchi, delta_mass=1., delta_charge=-1, start_position=2, end_position=10, comments='Long string')
        self.assertEqual(base.id, 'dAMP')
        self.assertEqual(base.name, 'deoxyadenosine monophosphate')
        self.assertEqual(base.synonyms, synonyms)
        self.assertEqual(base.identifiers, identifiers)
        self.assertEqual(base.get_inchi(), dAMP_inchi)
        self.assertEqual(base.delta_mass, 1.)
        self.assertEqual(base.delta_charge, -1)
        self.assertEqual(base.start_position, 2)
        self.assertEqual(base.end_position, 10)
        self.assertEqual(base.comments, 'Long string')

    def test_id_setter(self):
        base = core.Base()
        base.id = None
        base.id = ''
        base.id = 'A'
        with self.assertRaises(ValueError):
            base.id = 1

    def test_name_setter(self):
        base = core.Base()
        base.name = None
        base.name = ''
        base.name = 'A'
        with self.assertRaises(ValueError):
            base.name = 1

    def test_synonyms_setter(self):
        base = core.Base()
        base.synonyms = core.SynonymSet()
        base.synonyms = set(['A'])
        base.synonyms = ['A']
        with self.assertRaises(ValueError):
            base.synonyms = None
        with self.assertRaises(ValueError):
            base.synonyms = 'A'

    def test_identifiers_setter(self):
        base = core.Base()
        base.identifiers = core.IdentifierSet()
        base.identifiers = set([core.Identifier('ns', 'id')])
        base.identifiers = [core.Identifier('ns', 'id')]
        with self.assertRaises(ValueError):
            base.identifiers = None
        with self.assertRaises(ValueError):
            base.identifiers = 'A'
        with self.assertRaises(TypeError):
            base.identifiers = core.Identifier('ns', 'id')

    def test_structure_setter(self):
        base = core.Base()
        base.structure = dAMP_inchi

        ob_mol = openbabel.OBMol()
        conversion = openbabel.OBConversion()
        conversion.SetInFormat('inchi')
        conversion.ReadString(ob_mol, dAMP_inchi)
        base.structure = ob_mol

        base.structure = ''
        base.structure = None

        with self.assertRaises(ValueError):
            base.structure = 'InChI'

    def test_delta_mass_setter(self):
        base = core.Base()
        base.delta_mass = None
        base.delta_mass = 1
        base.delta_mass = 1.
        with self.assertRaises(ValueError):
            base.delta_mass = 'a'

    def test_delta_charge_setter(self):
        base = core.Base()
        base.delta_charge = None
        base.delta_charge = 1
        base.delta_charge = 1.
        with self.assertRaises(ValueError):
            base.delta_charge = 1.5
        with self.assertRaises(ValueError):
            base.delta_charge = 'a'

    def test_start_position_setter(self):
        base = core.Base()
        base.start_position = None
        base.start_position = 1
        base.start_position = 1.
        with self.assertRaises(ValueError):
            base.start_position = 1.5
        with self.assertRaises(ValueError):
            base.start_position = -1
        with self.assertRaises(ValueError):
            base.start_position = 'a'

    def test_end_position_setter(self):
        base = core.Base()
        base.end_position = None
        base.end_position = 1
        base.end_position = 1.
        with self.assertRaises(ValueError):
            base.end_position = 1.5
        with self.assertRaises(ValueError):
            base.end_position = -1
        with self.assertRaises(ValueError):
            base.end_position = 'a'

    def test_comments_setter(self):
        base = core.Base()
        base.comments = None
        base.comments = '1'
        with self.assertRaises(ValueError):
            base.comments = 1

    def test_protonate(self):
        base = core.Base()
        base.protonate(7.)

        base = core.Base(structure=dAMP_inchi)
        base.protonate(7.)
        base.protonate(10.)

    def test_get_inchi(self):
        base = core.Base()
        self.assertEqual(base.get_inchi(), None)

        base = core.Base(structure=dAMP_inchi)
        self.assertEqual(base.get_inchi(), dAMP_inchi)

    def test_get_formula(self):
        base = core.Base()
        self.assertEqual(base.get_formula(), None)

        base = core.Base(structure=dAMP_inchi)
        self.assertEqual(base.get_formula(), EmpiricalFormula('C10H12N5O6P'))

    def test_get_mol_wt(self):
        base = core.Base()
        self.assertEqual(base.get_mol_wt(), None)

        base = core.Base(structure=dAMP_inchi)
        self.assertEqual(base.get_mol_wt(), 329.208761998)

        base.delta_mass = 1.
        self.assertEqual(base.get_mol_wt(), 330.208761998)

    def test_get_charge(self):
        base = core.Base()
        self.assertEqual(base.get_charge(), None)

        base = core.Base(structure=dAMP_inchi)
        self.assertEqual(base.get_charge(), -2)

        base = core.Base(structure=dAMP_inchi)
        base.delta_charge = 1
        self.assertEqual(base.get_charge(), -1)

    def test_str(self):
        base = core.Base()

        base.id = 'dAMP'
        self.assertEqual(str(base), '[id: "dAMP"]')

        base.name = 'deoxyadenosine monophosphate'
        self.assertEqual(str(base), '[id: "dAMP" | name: "deoxyadenosine monophosphate"]')

        base.synonyms = set(['A', 'dAMP'])
        self.assertIn(' | synonym: "A"', str(base))
        self.assertIn(' | synonym: "dAMP"', str(base))

        base.identifiers = set([core.Identifier('chebi', 'CHEBI:58245'), core.Identifier('biocyc.compound', 'DAMP')])
        self.assertIn(' | identifier: chebi/CHEBI:58245', str(base))
        self.assertIn(' | identifier: biocyc.compound/DAMP', str(base))

        base.structure = dAMP_inchi
        self.assertIn(' | structure: {}]'.format(dAMP_inchi), str(base))

        base.delta_mass = 1.
        base.delta_charge = -1
        self.assertIn(' | delta-mass: 1', str(base))
        self.assertIn(' | delta-charge: -1', str(base))

        base.start_position = 3
        self.assertIn(' | position: 3-]', str(base))
        base.end_position = 5
        self.assertIn(' | position: 3-5]', str(base))
        base.start_position = None
        self.assertIn(' | position: -5]', str(base))

        base.comments = 'help "me"'
        self.assertIn(' | comments: "help \\"me\\""', str(base))

    def test_is_equal(self):
        base_1 = core.Base(id='A', structure=dAMP_inchi)
        base_2 = core.Base(id='A', structure=dAMP_inchi)
        base_3 = core.Base(id='B', structure=dAMP_inchi)
        base_4 = core.Base(id='A', structure=dCMP_inchi)

        self.assertTrue(base_1.is_equal(base_1))
        self.assertTrue(base_1.is_equal(base_2))
        self.assertTrue(base_2.is_equal(base_1))
        self.assertFalse(base_1.is_equal(mock.Mock(id='A', structure=dAMP_inchi)))
        self.assertFalse(base_1.is_equal(base_3))
        self.assertFalse(base_1.is_equal(base_4))


class BaseSequenceTestCase(unittest.TestCase):
    def test_init(self):
        seq = core.BaseSequence(None)
        seq = core.BaseSequence()
        self.assertEqual(len(seq), 0)

        seq = core.BaseSequence([core.Base(), core.Base()])
        self.assertEqual(len(seq), 2)

        with self.assertRaises(ValueError):
            core.BaseSequence('A')
        with self.assertRaises(ValueError):
            core.BaseSequence(['A'])

    def test_append(self):
        seq = core.BaseSequence()
        seq.append(core.Base())
        seq.append(core.Base())
        self.assertEqual(len(seq), 2)
        with self.assertRaises(ValueError):
            seq.append('A')

    def test_extend(self):
        seq = core.BaseSequence()
        seq.extend([core.Base(), core.Base()])
        self.assertEqual(len(seq), 2)
        with self.assertRaises(ValueError):
            seq.extend(['A'])

    def test_insert(self):
        seq = core.BaseSequence()
        seq.insert(0, core.Base())
        with self.assertRaises(ValueError):
            seq.insert(0, 'A')

    def test_setitem(self):
        seq = core.BaseSequence([core.Base(), core.Base()])

        seq[0] = core.Base()
        seq[0:1] = [core.Base()]
        seq[0:1] = core.BaseSequence([core.Base()])

        with self.assertRaises(ValueError):
            seq[0] = 'A'
        with self.assertRaises(ValueError):
            seq[0] = ['A']
        with self.assertRaises(ValueError):
            seq[0:1] = 'A'
        with self.assertRaises(ValueError):
            seq[0:1] = ['A']

    def test_get_base_counts(self):
        base_1 = core.Base(id='A')
        base_2 = core.Base(id='A')
        base_3 = core.Base(id='A')
        seq = core.BaseSequence([base_1, base_2, base_3, base_3, base_3, base_2, base_2, base_3])
        self.assertEqual(seq.get_base_counts(), {
            base_1: 1,
            base_2: 3,
            base_3: 4,
        })

    def test_is_equal(self):
        seq_1 = core.BaseSequence([core.Base(id='A'), core.Base(id='B')])
        seq_2 = core.BaseSequence([core.Base(id='A'), core.Base(id='B')])
        seq_3 = core.BaseSequence([core.Base(id='A'), core.Base(id='C')])
        self.assertTrue(seq_1.is_equal(seq_1))
        self.assertTrue(seq_1.is_equal(seq_2))
        self.assertFalse(seq_1.is_equal([]))
        self.assertFalse(seq_1.is_equal(seq_3))


class BpFormTestCase(unittest.TestCase):
    def test_init(self):
        bp_form = core.BpForm()
        self.assertEqual(bp_form.base_seq, core.BaseSequence())
        self.assertEqual(bp_form.alphabet, {})
        self.assertEqual(bp_form.bond_formula, EmpiricalFormula())
        self.assertEqual(bp_form.bond_charge, 0)

    def test_set_base_seq(self):
        bp_form = core.BpForm()

        bp_form.base_seq = core.BaseSequence()
        self.assertEqual(len(bp_form.base_seq), 0)

        bp_form.base_seq = [core.Base(), core.Base()]
        self.assertIsInstance(bp_form.base_seq, core.BaseSequence)
        self.assertEqual(len(bp_form.base_seq), 2)

        with self.assertRaises(ValueError):
            bp_form.base_seq = None
        with self.assertRaises(ValueError):
            bp_form.base_seq = 'A'

    def test_set_alphabet(self):
        bp_form = core.BpForm()

        bp_form.alphabet = dna.dna_alphabet
        self.assertEqual(len(bp_form.alphabet), 4)

        with self.assertRaises(ValueError):
            bp_form.alphabet = None
        with self.assertRaises(ValueError):
            bp_form.alphabet = 'A'

    def test_set_bond_formula(self):
        bp_form = core.BpForm()

        bp_form.bond_formula = EmpiricalFormula('CHO')
        self.assertEqual(bp_form.bond_formula, EmpiricalFormula('CHO'))

        bp_form.bond_formula = 'CHO'
        self.assertEqual(bp_form.bond_formula, EmpiricalFormula('CHO'))

        with self.assertRaises(ValueError):
            bp_form.bond_formula = '123'
        with self.assertRaises(TypeError):
            bp_form.bond_formula = 123
        with self.assertRaises(TypeError):
            bp_form.bond_formula = None

    def test_set_bond_charge(self):
        bp_form = core.BpForm()

        bp_form.bond_charge = 1
        self.assertEqual(bp_form.bond_charge, 1)

        bp_form.bond_charge = 1.
        self.assertEqual(bp_form.bond_charge, 1)

        bp_form.bond_charge = -1
        self.assertEqual(bp_form.bond_charge, -1)

        with self.assertRaises(ValueError):
            bp_form.bond_charge = 1.5

        with self.assertRaises(ValueError):
            bp_form.bond_charge = None

    def test_is_equal(self):
        bp_form_1 = core.BpForm(base_seq=core.BaseSequence([core.Base(id='A'), core.Base(id='B')]))
        bp_form_2 = core.BpForm(base_seq=core.BaseSequence([core.Base(id='A'), core.Base(id='B')]))
        bp_form_3 = None
        bp_form_4 = core.BpForm(base_seq=core.BaseSequence([core.Base(id='A'), core.Base(id='B')]), alphabet=dna.dna_alphabet)
        bp_form_5 = core.BpForm(base_seq=core.BaseSequence([core.Base(id='A'), core.Base(id='B')]), bond_charge=-1)
        bp_form_6 = core.BpForm(base_seq=core.BaseSequence(
            [core.Base(id='A'), core.Base(id='B')]), bond_formula=EmpiricalFormula('H'))
        self.assertTrue(bp_form_1.is_equal(bp_form_1))
        self.assertTrue(bp_form_1.is_equal(bp_form_2))
        self.assertFalse(bp_form_1.is_equal(bp_form_3))
        self.assertFalse(bp_form_1.is_equal(bp_form_4))
        self.assertFalse(bp_form_1.is_equal(bp_form_5))
        self.assertFalse(bp_form_1.is_equal(bp_form_6))

    def test_getitem(self):
        base_1 = core.Base(id='A')
        base_2 = core.Base(id='B')
        base_3 = core.Base(id='C')
        bp_form = core.BpForm([base_1, base_2, base_3])
        self.assertEqual(bp_form[0], base_1)
        self.assertEqual(bp_form[1], base_2)
        self.assertEqual(bp_form[0:1], [base_1])

    def test_setitem(self):
        base_1 = core.Base(id='A')
        base_2 = core.Base(id='B')
        base_3 = core.Base(id='C')
        bp_form = core.BpForm([base_1, base_2, base_3])

        self.assertEqual(bp_form[0], base_1)

        bp_form[0] = base_2
        self.assertEqual(bp_form[0], base_2)

        bp_form[0:1] = [base_3]
        self.assertEqual(bp_form[0], base_3)

    def test_delitem(self):
        base_1 = core.Base(id='A')
        base_2 = core.Base(id='B')
        base_3 = core.Base(id='C')
        bp_form = core.BpForm([base_1, base_2, base_3])
        del(bp_form[1])
        self.assertTrue(bp_form.is_equal(core.BpForm([base_1, base_3])))

    def test_iter(self):
        base_1 = core.Base(id='A')
        base_2 = core.Base(id='B')
        base_3 = core.Base(id='C')
        bp_form = core.BpForm([base_1, base_2, base_3])
        for i_base, base in enumerate(bp_form):
            if i_base == 0:
                self.assertEqual(base, base_1)
            if i_base == 1:
                self.assertEqual(base, base_2)
            if i_base == 2:
                self.assertEqual(base, base_3)

    def test_reversed(self):
        base_1 = core.Base(id='A')
        base_2 = core.Base(id='B')
        base_3 = core.Base(id='C')
        bp_form = core.BpForm([base_1, base_2, base_3])
        for i_base, base in enumerate(reversed(bp_form)):
            if i_base == 2:
                self.assertEqual(base, base_1)
            if i_base == 1:
                self.assertEqual(base, base_2)
            if i_base == 0:
                self.assertEqual(base, base_3)

    def test_contains(self):
        base_1 = core.Base(id='A')
        base_2 = core.Base(id='B')
        base_3 = core.Base(id='C')
        bp_form = core.BpForm([base_1, base_2])
        self.assertIn(base_1, bp_form)
        self.assertIn(base_2, bp_form)
        self.assertNotIn(base_3, bp_form)

    def test_len(self):
        bp_form = core.BpForm()
        self.assertEqual(len(bp_form), 0)

        bp_form = core.BpForm(base_seq=[core.Base(), core.Base()])
        self.assertEqual(len(bp_form), 2)

    def test_get_base_counts(self):
        base_1 = core.Base(id='A')
        base_2 = core.Base(id='B')
        base_3 = core.Base(id='C')
        bp_form = core.BpForm([base_1, base_2, base_1, base_1, base_1, base_2, base_2, base_3])
        self.assertEqual(bp_form.get_base_counts(), {
            base_1: 4,
            base_2: 3,
            base_3: 1,
        })

    def test_get_formula_mol_wt_charge(self):
        base_A = core.Base(id='A', structure=dAMP_inchi)
        base_C = core.Base(id='C', structure=dCMP_inchi)

        bp_form = core.BpForm([base_A])
        self.assertEqual(bp_form.get_formula(), base_A.get_formula())
        self.assertEqual(bp_form.get_mol_wt(), base_A.get_mol_wt())
        self.assertEqual(bp_form.get_charge(), base_A.get_charge())

        bp_form = core.BpForm([base_C])
        self.assertEqual(bp_form.get_formula(), base_C.get_formula())
        self.assertEqual(bp_form.get_mol_wt(), base_C.get_mol_wt())
        self.assertEqual(bp_form.get_charge(), base_C.get_charge())

        bp_form = core.BpForm([base_A, base_C])
        self.assertEqual(bp_form.get_formula(), base_A.get_formula() + base_C.get_formula())
        self.assertEqual(bp_form.get_mol_wt(), base_A.get_mol_wt() + base_C.get_mol_wt())
        self.assertEqual(bp_form.get_charge(), base_A.get_charge() + base_C.get_charge())

        bp_form = core.BpForm([base_A, base_C], bond_formula=EmpiricalFormula('H') * -1, bond_charge=1)
        self.assertEqual(bp_form.get_formula(), base_A.get_formula() + base_C.get_formula() - EmpiricalFormula('H'))
        self.assertEqual(bp_form.get_mol_wt(), base_A.get_mol_wt() + base_C.get_mol_wt() -
                         EmpiricalFormula('H').get_molecular_weight())
        self.assertEqual(bp_form.get_charge(), base_A.get_charge() + base_C.get_charge() + 1)

        bp_form = core.BpForm([base_A, base_A, base_C, base_C, base_C], bond_formula=EmpiricalFormula('H') * -1, bond_charge=1)
        self.assertEqual(bp_form.get_formula(), base_A.get_formula() * 2 + base_C.get_formula() * 3 - EmpiricalFormula('H') * 4)
        self.assertEqual(bp_form.get_mol_wt(), base_A.get_mol_wt() * 2 + base_C.get_mol_wt()
                         * 3 - EmpiricalFormula('H').get_molecular_weight() * 4)
        self.assertEqual(bp_form.get_charge(), base_A.get_charge() * 2 + base_C.get_charge() * 3 + 1 * 4)

    def test_protonate(self):
        base_1 = core.Base(id='A')
        base_2 = core.Base(id='B')
        bp_form = core.BpForm([base_1, base_2])
        bp_form.protonate(7.)

    def test_str(self):
        base_A = core.Base(id='A', structure=dAMP_inchi)
        base_C = core.Base(id='C', structure=dCMP_inchi)
        base_G = core.Base(id='G', structure=dGMP_inchi)

        bp_form = core.BpForm([base_A, base_C, base_G, base_A])
        self.assertEqual(str(bp_form), '{}{}{}{}'.format(str(base_A), str(base_C), str(base_G), str(base_A)))

        bp_form = core.BpForm([base_A, base_C, base_G, base_A], alphabet=core.Alphabet({
            'A': base_A,
            'C': base_C,
        }))
        self.assertEqual(str(bp_form), '{}{}{}{}'.format('A', 'C', str(base_G), 'A'))
        self.assertEqual(str(bp_form), '{}{}{}{}'.format('A', 'C', '[id: "{}" | structure: {}]'.format('G', dGMP_inchi), 'A'))

    def test_from_str(self):
        self.assertTrue(dna.DnaForm.from_str('AAA').is_equal(dna.DnaForm([
            dna.dna_alphabet.A, dna.dna_alphabet.A, dna.dna_alphabet.A,
        ])))

        self.assertTrue(dna.DnaForm.from_str('ACTG').is_equal(dna.DnaForm([
            dna.dna_alphabet.A, dna.dna_alphabet.C, dna.dna_alphabet.T, dna.dna_alphabet.G,
        ])))

        with self.assertRaisesRegex(lark.exceptions.VisitError, 'not in alphabet'):
            self.assertTrue(dna.DnaForm.from_str('UAA').is_equal(dna.DnaForm([
                dna.dna_alphabet.A, dna.dna_alphabet.A, dna.dna_alphabet.A,
            ])))

        self.assertTrue(dna.DnaForm.from_str(
            'AA[id: "dI"'
            + ' | name: "2\'-deoxyinosine"'
            + ' | synonym: "2\'-deoxyinosine, 9-[(2R,4S,5R)-4-hydroxy-5-(hydroxymethyl)tetrahydrofuran-2-yl]-9H-purin-6-ol"'
            + ' | identifier: chebi/CHEBI:28997'
            + ' | structure: ' + dIMP_inchi
            + ' | delta-mass: -2.5'
            + ' | delta-charge: 3'
            + ' | position: 3-5'
            + ' | comments: "A purine 2\'-deoxyribonucleoside that is inosine ..."]A').is_equal(dna.DnaForm([
                dna.dna_alphabet.A,
                dna.dna_alphabet.A,
                core.Base(
                    id='dI',
                    name="2'-deoxyinosine",
                    synonyms=core.SynonymSet(
                        ["2'-deoxyinosine, 9-[(2R,4S,5R)-4-hydroxy-5-(hydroxymethyl)tetrahydrofuran-2-yl]-9H-purin-6-ol"]),
                    identifiers=core.IdentifierSet([core.Identifier('chebi', 'CHEBI:28997')]),
                    structure=dIMP_inchi,
                    delta_mass=-2.5,
                    delta_charge=3,
                    start_position=3,
                    end_position=5,
                    comments="A purine 2'-deoxyribonucleoside that is inosine ...",
                ),
                dna.dna_alphabet.A,
            ])))

        self.assertTrue(dna.DnaForm.from_str(
            'AA[id: "dI"'
            ' | position: 3-]A').is_equal(dna.DnaForm([
                dna.dna_alphabet.A,
                dna.dna_alphabet.A,
                core.Base(
                    id='dI',
                    start_position=3,
                    end_position=None,
                ),
                dna.dna_alphabet.A,
            ])))

        self.assertTrue(dna.DnaForm.from_str(
            'AA[id: "dI"'
            ' | position: -5]A').is_equal(dna.DnaForm([
                dna.dna_alphabet.A,
                dna.dna_alphabet.A,
                core.Base(
                    id='dI',
                    start_position=None,
                    end_position=5,
                ),
                dna.dna_alphabet.A,
            ])))

        with self.assertRaisesRegex(lark.exceptions.VisitError, 'cannot be repeated'):
            dna.DnaForm.from_str(
                'AA[id: "dI"'
                ' | name: "2\'-deoxyinosine"'
                ' | name: "2\'-deoxyinosine"]A')


class AlphabetTestCase(unittest.TestCase):
    def setUp(self):
        self.dir_path = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir_path)

    def test_getitem(self):
        self.assertEqual(dna.dna_alphabet.A.get_inchi(), dAMP_inchi)

    def test_protonate(self):
        alphabet = core.Alphabet({
            'A': core.Base(structure=dAMP_inchi),
            'A': core.Base(structure=dCMP_inchi),
        })
        alphabet.protonate(7.)

    def test_is_equal(self):
        self.assertTrue(dna.dna_alphabet.is_equal(dna.dna_alphabet))
        self.assertFalse(dna.dna_alphabet.is_equal(rna.rna_alphabet))
        self.assertFalse(dna.dna_alphabet.is_equal(None))
        self.assertFalse(dna.dna_alphabet.is_equal(core.Alphabet()))

        dna_alphabet = core.Alphabet({
            'A': dna.dna_alphabet.A,
            'C': dna.dna_alphabet.C,
            'G': dna.dna_alphabet.G,
            'T': dna.dna_alphabet.T,
        })
        self.assertTrue(dna_alphabet.is_equal(dna.dna_alphabet))

    def test_to_from_yaml(self):
        dna_alphabet = core.Alphabet({
            'A': dna.dna_alphabet.A,
            'C': dna.dna_alphabet.C,
            'G': dna.dna_alphabet.G,
            'T': dna.dna_alphabet.T,
        })
        path = os.path.join(self.dir_path, 'alphabet.yml')
        dna_alphabet.to_yaml(path)
        dna_alphabet_2 = core.Alphabet.from_yaml(path)
        self.assertTrue(dna_alphabet_2.is_equal(dna_alphabet))
