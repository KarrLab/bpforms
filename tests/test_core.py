""" Test of bpforms.core

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms import core
from bpforms.alphabet import dna
from bpforms.alphabet import protein
from bpforms.alphabet import rna
from wc_utils.util.chem import EmpiricalFormula
import copy
import lark.exceptions
import mock
import openbabel
import os
import requests
import shutil
import tempfile
import unittest

dAMP_inchi = dna.canonical_dna_alphabet.monomers.A.export('inchi')
dCMP_inchi = dna.canonical_dna_alphabet.monomers.C.export('inchi')
dGMP_inchi = dna.canonical_dna_alphabet.monomers.G.export('inchi')
dAMP_smiles = dna.canonical_dna_alphabet.monomers.A.export('smiles')
dCMP_smiles = dna.canonical_dna_alphabet.monomers.C.export('smiles')
dGMP_smiles = dna.canonical_dna_alphabet.monomers.G.export('smiles')
dIMP_smiles = 'OCC1OC(CC1O)N1C=NC2=C1N=CN=C2O'


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
            id.ns = 0
        with self.assertRaises(ValueError):
            id.ns = None
        with self.assertRaises(ValueError):
            id.id = ''
        with self.assertRaises(ValueError):
            id.id = 0
        with self.assertRaises(ValueError):
            id.id = None
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


class MonomerTestCase(unittest.TestCase):
    def test_init(self):
        identifiers = set([
            core.Identifier('chebi', 'CHEBI:58245'),
            core.Identifier('pubchem.compound', '22848660'),
            core.Identifier('metacyc.compound', 'DAMP'),
        ])
        synonyms = set(['A', 'dAMP', 'deoxyadenosine monophosphate'])
        monomer_0 = core.Monomer()
        monomer = core.Monomer(id='dAMP', name='deoxyadenosine monophosphate', synonyms=synonyms, identifiers=identifiers,
                               structure=dAMP_smiles, delta_mass=1., delta_charge=-1, start_position=2, end_position=10,
                               base_monomers=[monomer_0],
                               comments='Long string')
        self.assertEqual(monomer.id, 'dAMP')
        self.assertEqual(monomer.name, 'deoxyadenosine monophosphate')
        self.assertEqual(monomer.synonyms, synonyms)
        self.assertEqual(monomer.identifiers, identifiers)
        self.assertEqual(monomer.export('inchi'), dAMP_inchi)
        self.assertEqual(monomer.delta_mass, 1.)
        self.assertEqual(monomer.delta_charge, -1)
        self.assertEqual(monomer.start_position, 2)
        self.assertEqual(monomer.end_position, 10)
        self.assertEqual(monomer.base_monomers, set([monomer_0]))
        self.assertEqual(monomer.comments, 'Long string')

    def test_id_setter(self):
        monomer = core.Monomer()
        monomer.id = None
        monomer.id = ''
        monomer.id = 'A'
        with self.assertRaises(ValueError):
            monomer.id = 1

    def test_name_setter(self):
        monomer = core.Monomer()
        monomer.name = None
        monomer.name = ''
        monomer.name = 'A'
        with self.assertRaises(ValueError):
            monomer.name = 1

    def test_synonyms_setter(self):
        monomer = core.Monomer()
        monomer.synonyms = core.SynonymSet()
        monomer.synonyms = set(['A'])
        monomer.synonyms = ['A']
        with self.assertRaises(ValueError):
            monomer.synonyms = None
        with self.assertRaises(ValueError):
            monomer.synonyms = 'A'

    def test_identifiers_setter(self):
        monomer = core.Monomer()
        monomer.identifiers = core.IdentifierSet()
        monomer.identifiers = set([core.Identifier('ns', 'id')])
        monomer.identifiers = [core.Identifier('ns', 'id')]
        with self.assertRaises(ValueError):
            monomer.identifiers = None
        with self.assertRaises(ValueError):
            monomer.identifiers = 'A'
        with self.assertRaises(TypeError):
            monomer.identifiers = core.Identifier('ns', 'id')

    def test_structure_setter(self):
        monomer = core.Monomer()
        monomer.structure = dAMP_smiles

        ob_mol = openbabel.OBMol()
        conversion = openbabel.OBConversion()
        conversion.SetInFormat('smi')
        conversion.ReadString(ob_mol, dAMP_smiles)
        monomer.structure = ob_mol

        monomer.structure = ''
        monomer.structure = None

        with self.assertRaises(ValueError):
            monomer.structure = 'InChI'

    def test_delta_mass_setter(self):
        monomer = core.Monomer()
        monomer.delta_mass = None
        monomer.delta_mass = 1
        monomer.delta_mass = 1.
        with self.assertRaises(ValueError):
            monomer.delta_mass = 'a'

    def test_delta_charge_setter(self):
        monomer = core.Monomer()
        monomer.delta_charge = None
        monomer.delta_charge = 1
        monomer.delta_charge = 1.
        with self.assertRaises(ValueError):
            monomer.delta_charge = 1.5
        with self.assertRaises(ValueError):
            monomer.delta_charge = 'a'

    def test_start_position_setter(self):
        monomer = core.Monomer()
        monomer.start_position = None
        monomer.start_position = 1
        monomer.start_position = 1.
        with self.assertRaises(ValueError):
            monomer.start_position = 1.5
        with self.assertRaises(ValueError):
            monomer.start_position = -1
        with self.assertRaises(ValueError):
            monomer.start_position = 'a'

    def test_end_position_setter(self):
        monomer = core.Monomer()
        monomer.end_position = None
        monomer.end_position = 1
        monomer.end_position = 1.
        with self.assertRaises(ValueError):
            monomer.end_position = 1.5
        with self.assertRaises(ValueError):
            monomer.end_position = -1
        with self.assertRaises(ValueError):
            monomer.end_position = 'a'

    def test_base_monomers_setter(self):
        monomer = core.Monomer()
        monomer.base_monomers = []
        monomer.base_monomers = set([core.Monomer()])
        with self.assertRaises(ValueError):
            monomer.base_monomers = 'A'

    def test_set_monomer_bond_atoms(self):
        monomer = core.Monomer()
        monomer.monomer_bond_atoms = []
        monomer.monomer_bond_atoms = core.AtomList()
        with self.assertRaises(ValueError):
            monomer.monomer_bond_atoms = None

    def test_set_monomer_displaced_atoms(self):
        monomer = core.Monomer()
        monomer.monomer_displaced_atoms = []
        monomer.monomer_displaced_atoms = core.AtomList()
        with self.assertRaises(ValueError):
            monomer.monomer_displaced_atoms = None

    def test_set_left_bond_atoms(self):
        monomer = core.Monomer()
        monomer.left_bond_atoms = []
        monomer.left_bond_atoms = core.AtomList()
        with self.assertRaises(ValueError):
            monomer.left_bond_atoms = None

    def test_set_bond_bond_atoms(self):
        monomer = core.Monomer()
        monomer.right_bond_atoms = []
        monomer.right_bond_atoms = core.AtomList()
        with self.assertRaises(ValueError):
            monomer.right_bond_atoms = None

    def test_set_left_displaced_atoms(self):
        monomer = core.Monomer()
        monomer.left_displaced_atoms = []
        monomer.left_displaced_atoms = core.AtomList()
        with self.assertRaises(ValueError):
            monomer.left_displaced_atoms = None

    def test_set_bond_displaced_atoms(self):
        monomer = core.Monomer()
        monomer.right_displaced_atoms = []
        monomer.right_displaced_atoms = core.AtomList()
        with self.assertRaises(ValueError):
            monomer.right_displaced_atoms = None

    def test_comments_setter(self):
        monomer = core.Monomer()
        monomer.comments = None
        monomer.comments = '1'
        with self.assertRaises(ValueError):
            monomer.comments = 1

    def test_get_major_micro_species(self):
        monomer = core.Monomer()
        monomer.get_major_micro_species(7.)

        monomer = core.Monomer(structure=dAMP_smiles)
        monomer.get_major_micro_species(7.)
        monomer.get_major_micro_species(10.)

    def test_export(self):
        monomer = core.Monomer()
        self.assertEqual(monomer.export('inchi'), None)

        monomer = core.Monomer(structure=dAMP_smiles)
        self.assertEqual(monomer.export('inchi'), dAMP_inchi)

    def test_get_formula(self):
        monomer = core.Monomer(structure=dAMP_smiles)
        self.assertEqual(monomer.get_formula(), EmpiricalFormula('C5H5N5'))

        with self.assertRaises(ValueError):
            monomer = core.Monomer()
            monomer.get_formula()

    def test_get_mol_wt(self):
        monomer = core.Monomer()
        self.assertEqual(monomer.get_mol_wt(), None)

        monomer = core.Monomer(structure=dAMP_smiles)
        self.assertEqual(monomer.get_mol_wt(), 135.13)

        monomer.delta_mass = 1.
        self.assertEqual(monomer.get_mol_wt(), 136.13)

    def test_get_charge(self):
        monomer = core.Monomer(structure=dAMP_smiles)
        self.assertEqual(monomer.get_charge(), 0)

        monomer = core.Monomer(structure=dAMP_smiles)
        monomer.delta_charge = 1
        self.assertEqual(monomer.get_charge(), 1)

        with self.assertRaises(ValueError):
            monomer = core.Monomer()
            monomer.get_charge()

    def test_to_from_dict(self):
        alphabet = dna.dna_alphabet

        monomer = core.Monomer()
        monomer.base_monomers = [alphabet.monomers.A]
        monomer.monomer_bond_atoms = [core.Atom(core.Monomer, 'C', charge=3), core.Atom(core.Monomer, 'H', position=3, charge=2)]
        monomer_as_dict = monomer.to_dict(alphabet=alphabet)
        self.assertEqual(monomer_as_dict, {
            'base_monomers': ['A'],
            'monomer_bond_atoms': [
                {'molecule': 'Monomer', 'element': 'C', 'charge': 3},
                {'molecule': 'Monomer', 'element': 'H', 'position': 3, 'charge': 2},
            ],
        })

        monomer_2 = core.Monomer()
        monomer_2.from_dict(monomer_as_dict, alphabet=alphabet)
        self.assertEqual(monomer_2.base_monomers, set([alphabet.monomers.A]))
        self.assertTrue(monomer.is_equal(monomer_2))

    def test_str(self):
        monomer = core.Monomer()
        self.assertEqual(str(monomer), '[]')

        monomer.id = 'dAMP'
        self.assertEqual(str(monomer), '[id: "dAMP"]')

        monomer.name = 'deoxyadenosine monophosphate'
        self.assertEqual(str(monomer), '[id: "dAMP" | name: "deoxyadenosine monophosphate"]')

        monomer.synonyms = set(['A', 'dAMP'])
        self.assertIn(' | synonym: "A"', str(monomer))
        self.assertIn(' | synonym: "dAMP"', str(monomer))

        monomer.identifiers = set([core.Identifier('chebi', 'CHEBI:58245'), core.Identifier('biocyc.compound', 'DAMP')])
        self.assertIn(' | identifier: "chebi" / "CHEBI:58245"', str(monomer))
        self.assertIn(' | identifier: "biocyc.compound" / "DAMP"', str(monomer))

        monomer.structure = dAMP_smiles
        self.assertIn(' | structure: "{}"]'.format('Nc1ncnc2c1nc[nH]2'), str(monomer))

        monomer.monomer_bond_atoms.append(core.Atom(core.Monomer, 'C', 2, -3))
        self.assertIn(' | monomer-bond-atom: Monomer / C / 2 / -3]', str(monomer))

        monomer.delta_mass = 1.
        monomer.delta_charge = -1
        self.assertIn(' | delta-mass: 1', str(monomer))
        self.assertIn(' | delta-charge: -1', str(monomer))

        monomer.start_position = 3
        self.assertIn(' | position: 3-]', str(monomer))
        monomer.end_position = 5
        self.assertIn(' | position: 3-5]', str(monomer))
        monomer.start_position = None
        self.assertIn(' | position: -5]', str(monomer))

        monomer.comments = 'help "me"'
        self.assertIn(' | comments: "help \\"me\\""', str(monomer))

    def test_is_equal(self):
        monomer_1 = core.Monomer(id='A', structure=dAMP_smiles)
        monomer_2 = core.Monomer(id='A', structure=dAMP_smiles)
        monomer_3 = core.Monomer(id='B', structure=dAMP_smiles)
        monomer_4 = core.Monomer(id='A', structure=dCMP_smiles)
        monomer_5 = core.Monomer(id='A', structure=dAMP_smiles, base_monomers=[core.Monomer(id='A')])
        monomer_6 = core.Monomer(id='A', structure=dAMP_smiles, base_monomers=[core.Monomer(id='A')])
        monomer_7 = core.Monomer(id='A', structure=dAMP_smiles, base_monomers=[core.Monomer(id='B')])
        monomer_8 = core.Monomer(id='A', structure=dAMP_smiles, monomer_bond_atoms=[core.Atom(None, 'S')])

        self.assertTrue(monomer_1.is_equal(monomer_1))
        self.assertTrue(monomer_1.is_equal(monomer_2))
        self.assertTrue(monomer_2.is_equal(monomer_1))
        self.assertFalse(monomer_1.is_equal(mock.Mock(id='A', structure=dAMP_smiles)))
        self.assertFalse(monomer_1.is_equal(monomer_3))
        self.assertFalse(monomer_1.is_equal(monomer_4))
        self.assertFalse(monomer_1.is_equal(monomer_5))
        self.assertTrue(monomer_5.is_equal(monomer_6))
        self.assertFalse(monomer_5.is_equal(monomer_7))
        self.assertFalse(monomer_1.is_equal(monomer_8))

    def test_get_image_url(self):
        url = dna.dna_alphabet.monomers.A.get_image_url()
        self.assertNotEqual(url, None)
        response = requests.get(url)
        response.raise_for_status()

        self.assertEqual(core.Monomer().get_image_url(), None)

    def test_get_image(self):
        svg = dna.dna_alphabet.monomers.A.get_image()
        self.assertTrue(svg.startswith('<?xml'))

        svg = dna.dna_alphabet.monomers.A.get_image(include_xml_header=False)
        self.assertTrue(svg.startswith('<svg'))

        self.assertEqual(core.Monomer().get_image(), None)


class MonomerSequenceTestCase(unittest.TestCase):
    def test_init(self):
        seq = core.MonomerSequence(None)
        seq = core.MonomerSequence()
        self.assertEqual(len(seq), 0)

        seq = core.MonomerSequence([core.Monomer(), core.Monomer()])
        self.assertEqual(len(seq), 2)

        with self.assertRaises(ValueError):
            core.MonomerSequence('A')
        with self.assertRaises(ValueError):
            core.MonomerSequence(['A'])

    def test_append(self):
        seq = core.MonomerSequence()
        seq.append(core.Monomer())
        seq.append(core.Monomer())
        self.assertEqual(len(seq), 2)
        with self.assertRaises(ValueError):
            seq.append('A')

    def test_extend(self):
        seq = core.MonomerSequence()
        seq.extend([core.Monomer(), core.Monomer()])
        self.assertEqual(len(seq), 2)
        with self.assertRaises(ValueError):
            seq.extend(['A'])

    def test_insert(self):
        seq = core.MonomerSequence()
        seq.insert(0, core.Monomer())
        with self.assertRaises(ValueError):
            seq.insert(0, 'A')

    def test_setitem(self):
        seq = core.MonomerSequence([core.Monomer(), core.Monomer()])

        seq[0] = core.Monomer()
        seq[0:1] = [core.Monomer()]
        seq[0:1] = core.MonomerSequence([core.Monomer()])

        with self.assertRaises(ValueError):
            seq[0] = 'A'
        with self.assertRaises(ValueError):
            seq[0] = ['A']
        with self.assertRaises(ValueError):
            seq[0:1] = 'A'
        with self.assertRaises(ValueError):
            seq[0:1] = ['A']

    def test_get_monomer_counts(self):
        monomer_1 = core.Monomer(id='A')
        monomer_2 = core.Monomer(id='A')
        monomer_3 = core.Monomer(id='A')
        seq = core.MonomerSequence([monomer_1, monomer_2, monomer_3, monomer_3, monomer_3, monomer_2, monomer_2, monomer_3])
        self.assertEqual(seq.get_monomer_counts(), {
            monomer_1: 1,
            monomer_2: 3,
            monomer_3: 4,
        })

    def test_is_equal(self):
        seq_1 = core.MonomerSequence([core.Monomer(id='A'), core.Monomer(id='B')])
        seq_2 = core.MonomerSequence([core.Monomer(id='A'), core.Monomer(id='B')])
        seq_3 = core.MonomerSequence([core.Monomer(id='A'), core.Monomer(id='C')])
        self.assertTrue(seq_1.is_equal(seq_1))
        self.assertTrue(seq_1.is_equal(seq_2))
        self.assertFalse(seq_1.is_equal([]))
        self.assertFalse(seq_1.is_equal(seq_3))


class AtomTestCase(unittest.TestCase):
    def test_molecule_setter(self):
        atom = core.Atom(core.Monomer, 'C')

        atom.molecule = None
        self.assertEqual(atom.molecule, None)

        atom.molecule = core.Monomer
        self.assertEqual(atom.molecule, core.Monomer)

        atom.molecule = core.Backbone
        self.assertEqual(atom.molecule, core.Backbone)

        with self.assertRaises(ValueError):
            atom.molecule = core.Atom

    def test_element_setter(self):
        atom = core.Atom(core.Monomer, 'C')
        atom.element = 'C'
        self.assertEqual(atom.element, 'C')
        with self.assertRaises(ValueError):
            atom.element = 1

    def test_position_setter(self):
        atom = core.Atom(core.Monomer, 'C')

        atom.position = 2
        self.assertEqual(atom.position, 2)

        atom.position = 2.
        self.assertEqual(atom.position, 2)

        atom.position = None
        self.assertEqual(atom.position, None)

        with self.assertRaises(ValueError):
            atom.position = 'a'
        with self.assertRaises(ValueError):
            atom.position = 2.5
        with self.assertRaises(ValueError):
            atom.position = -1

    def test_charge_setter(self):
        atom = core.Atom(core.Monomer, 'C')

        atom.charge = 2
        self.assertEqual(atom.charge, 2)

        atom.charge = -3
        self.assertEqual(atom.charge, -3)

        with self.assertRaises(ValueError):
            atom.charge = None

        with self.assertRaises(ValueError):
            atom.charge = 'a'

        with self.assertRaises(ValueError):
            atom.charge = 2.5

    def test_is_equal(self):
        atom_1 = core.Atom(core.Monomer, 'C', position=2, charge=-3)
        self.assertTrue(atom_1.is_equal(atom_1))
        self.assertTrue(atom_1.is_equal(core.Atom(core.Monomer, 'C', position=2, charge=-3)))
        self.assertFalse(atom_1.is_equal({}))
        self.assertFalse(atom_1.is_equal(core.Atom(core.Monomer, 'H', position=2, charge=-3)))
        self.assertFalse(atom_1.is_equal(core.Atom(core.Monomer, 'C', position=3, charge=-3)))
        self.assertFalse(atom_1.is_equal(core.Atom(core.Monomer, 'C', position=2, charge=-2)))

    def test_to_from_dict(self):
        atom_1 = core.Atom(core.Monomer, 'C', position=None, charge=-3)
        atom_1_dict = atom_1.to_dict()
        self.assertEqual(atom_1_dict, {'molecule': 'Monomer', 'element': 'C', 'charge': -3})
        atom_2 = core.Atom(None, '').from_dict(atom_1_dict)
        self.assertTrue(atom_1.is_equal(atom_2))

        atom_1 = core.Atom(core.Monomer, 'C', position=2, charge=-3)
        atom_1_dict = atom_1.to_dict()
        self.assertEqual(atom_1_dict, {'molecule': 'Monomer', 'element': 'C', 'charge': -3, 'position': 2})
        atom_2 = core.Atom(None, '').from_dict(atom_1.to_dict())
        self.assertTrue(atom_1.is_equal(atom_2))

        atom_1 = core.Atom(core.Backbone, 'C', position=2, charge=-3)
        atom_1_dict = atom_1.to_dict()
        self.assertEqual(atom_1_dict, {'molecule': 'Backbone', 'element': 'C', 'charge': -3, 'position': 2})
        atom_2 = core.Atom(None, '').from_dict(atom_1.to_dict())
        self.assertTrue(atom_1.is_equal(atom_2))

        atom_1 = core.Atom(None, 'C', position=2, charge=-3)
        atom_1_dict = atom_1.to_dict()
        self.assertEqual(atom_1_dict, {'element': 'C', 'charge': -3, 'position': 2})
        atom_2 = core.Atom(None, '').from_dict(atom_1.to_dict())
        self.assertTrue(atom_1.is_equal(atom_2))


class AtomListTestCase(unittest.TestCase):
    def test_init(self):
        atom_1 = core.Atom(core.Monomer, 'C')
        atom_2 = core.Atom(core.Monomer, 'H')
        atom_list = core.AtomList([atom_1, atom_2])

    def test_append(self):
        atom_1 = core.Atom(core.Monomer, 'C')
        atom_list = core.AtomList()
        atom_list.append(atom_1)
        self.assertEqual(atom_list, core.AtomList([atom_1]))
        with self.assertRaises(ValueError):
            atom_list.append('C')

    def test_extend(self):
        atom_1 = core.Atom(core.Monomer, 'C')
        atom_2 = core.Atom(core.Monomer, 'H')
        atom_list = core.AtomList()
        atom_list.extend([atom_1, atom_2])
        self.assertEqual(atom_list, core.AtomList([atom_1, atom_2]))

    def test_insert(self):
        atom_1 = core.Atom(core.Monomer, 'C')
        atom_2 = core.Atom(core.Monomer, 'H')
        atom_3 = core.Atom(core.Monomer, 'O')
        atom_list = core.AtomList([atom_1, atom_2])

        atom_list.insert(1, atom_3)
        self.assertEqual(atom_list, core.AtomList([atom_1, atom_3, atom_2]))

        with self.assertRaises(ValueError):
            atom_list.insert(1, 'C')

    def test_set_item(self):
        atom_1 = core.Atom(core.Monomer, 'C')
        atom_2 = core.Atom(core.Monomer, 'H')
        atom_3 = core.Atom(core.Monomer, 'O')

        atom_list = core.AtomList([atom_1, atom_2, atom_3])
        atom_list[0] = atom_3
        self.assertEqual(atom_list, core.AtomList([atom_3, atom_2, atom_3]))

        atom_list = core.AtomList([atom_1, atom_2, atom_3])
        atom_list[0:1] = [atom_3]
        self.assertEqual(atom_list, core.AtomList([atom_3, atom_2, atom_3]))

        atom_list = core.AtomList([atom_1, atom_2, atom_3])
        with self.assertRaises(ValueError):
            atom_list[0] = 'C'

        atom_list = core.AtomList([atom_1, atom_2, atom_3])
        with self.assertRaises(ValueError):
            atom_list[0:1] = ['C']

    def test_is_equal(self):
        atom_1 = core.Atom(core.Monomer, 'C')
        atom_2 = core.Atom(core.Monomer, 'H')
        atom_3 = core.Atom(core.Monomer, 'O')

        atom_list_1 = core.AtomList([core.Atom(core.Monomer, 'C'), core.Atom(core.Monomer, 'H'), core.Atom(core.Monomer, 'O')])
        atom_list_2 = core.AtomList([core.Atom(core.Monomer, 'C'), core.Atom(core.Monomer, 'H'), core.Atom(core.Monomer, 'O')])
        atom_list_3 = core.AtomList([core.Atom(core.Monomer, 'C'), core.Atom(core.Monomer, 'N'), core.Atom(core.Monomer, 'O')])
        atom_list_4 = core.AtomList([core.Atom(core.Backbone, 'C'), core.Atom(core.Monomer, 'H'), core.Atom(core.Monomer, 'O')])
        self.assertTrue(atom_list_1.is_equal(atom_list_1))
        self.assertTrue(atom_list_1.is_equal(atom_list_2))
        self.assertFalse(atom_list_1.is_equal({}))
        self.assertFalse(atom_list_1.is_equal(atom_list_3))
        self.assertFalse(atom_list_1.is_equal(atom_list_4))

    def test_to_from_list(self):
        atom_list_1 = core.AtomList([core.Atom(core.Monomer, 'C'), core.Atom(core.Monomer, 'H'), core.Atom(core.Monomer, 'O')])
        atom_list_1_as_list = atom_list_1.to_list()
        self.assertEqual(atom_list_1_as_list, [
            {'molecule': 'Monomer', 'element': 'C'},
            {'molecule': 'Monomer', 'element': 'H'},
            {'molecule': 'Monomer', 'element': 'O'},
        ])
        atom_list_2 = core.AtomList().from_list(atom_list_1_as_list)
        self.assertTrue(atom_list_1.is_equal(atom_list_2))


class BackboneTestCase(unittest.TestCase):
    def test_set_structure(self):
        backbone = core.Backbone()
        backbone.structure = None
        backbone.structure = dAMP_smiles
        with self.assertRaises(ValueError):
            backbone.structure = 'dAMP'

    def test_set_monomer_bond_atoms(self):
        backbone = core.Backbone()
        backbone.monomer_bond_atoms = []
        backbone.monomer_bond_atoms = core.AtomList()
        with self.assertRaises(ValueError):
            backbone.monomer_bond_atoms = None

    def test_set_backbone_bond_atoms(self):
        backbone = core.Backbone()
        backbone.backbone_bond_atoms = []
        backbone.backbone_bond_atoms = core.AtomList()
        with self.assertRaises(ValueError):
            backbone.backbone_bond_atoms = None

    def test_set_monomer_displaced_atoms(self):
        backbone = core.Backbone()
        backbone.monomer_displaced_atoms = []
        backbone.monomer_displaced_atoms = core.AtomList()
        with self.assertRaises(ValueError):
            backbone.monomer_displaced_atoms = None

    def test_set_backbone_displaced_atoms(self):
        backbone = core.Backbone()
        backbone.backbone_displaced_atoms = []
        backbone.backbone_displaced_atoms = core.AtomList()
        with self.assertRaises(ValueError):
            backbone.backbone_displaced_atoms = None

    def test_export(self):
        backbone = core.Backbone()

        backbone.structure = dAMP_smiles
        self.assertEqual(backbone.export('inchi'), dAMP_inchi)
        self.assertEqual(backbone.export('smiles'), dAMP_smiles)

        backbone.structure = None
        self.assertEqual(backbone.export('inchi'), None)
        self.assertEqual(backbone.export('smiles'), None)

    def test_get_formula(self):
        backbone = core.Backbone()

        backbone.structure = dAMP_smiles
        self.assertEqual(backbone.get_formula(), EmpiricalFormula('C5H5N5'))

        backbone.monomer_displaced_atoms = core.AtomList([core.Atom(core.Monomer, 'C')])
        backbone.backbone_displaced_atoms = core.AtomList([core.Atom(core.Monomer, 'H'), core.Atom(core.Monomer, 'H')])
        self.assertEqual(backbone.get_formula(), EmpiricalFormula('C4H3N5'))

        backbone.structure = None
        self.assertEqual(backbone.get_formula(), EmpiricalFormula('CH2') * -1)

    def test_get_mol_wt(self):
        backbone = core.Backbone()

        backbone.structure = dAMP_smiles
        self.assertEqual(backbone.get_mol_wt(), 135.13)

    def test_get_charge(self):
        backbone = core.Backbone()

        backbone.structure = dAMP_smiles
        self.assertEqual(backbone.get_charge(), 0)

        backbone.structure = None
        self.assertEqual(backbone.get_charge(), 0)

        backbone.monomer_displaced_atoms = core.AtomList([core.Atom(core.Monomer, 'C', charge=2)])
        backbone.backbone_displaced_atoms = core.AtomList([core.Atom(core.Monomer, 'H', charge=3)])
        self.assertEqual(backbone.get_charge(), -5)

    def test_is_equal(self):
        backbone_1 = core.Backbone()
        backbone_2 = core.Backbone()
        backbone_3 = core.Backbone(structure=dAMP_smiles)
        backbone_4 = core.Backbone(monomer_bond_atoms=[core.Atom(core.Monomer, 'H')])
        backbone_5 = core.Backbone(backbone_bond_atoms=[core.Atom(core.Monomer, 'H')])
        backbone_6 = core.Backbone(monomer_displaced_atoms=[core.Atom(core.Monomer, 'H')])
        backbone_7 = core.Backbone(backbone_displaced_atoms=[core.Atom(core.Monomer, 'H')])
        self.assertTrue(backbone_1.is_equal(backbone_1))
        self.assertTrue(backbone_1.is_equal(backbone_2))
        self.assertFalse(backbone_1.is_equal({}))
        self.assertFalse(backbone_1.is_equal(backbone_3))
        self.assertFalse(backbone_1.is_equal(backbone_4))
        self.assertFalse(backbone_1.is_equal(backbone_5))
        self.assertFalse(backbone_1.is_equal(backbone_6))
        self.assertFalse(backbone_1.is_equal(backbone_7))


class BondTestCase(unittest.TestCase):
    def test_set_left_bond_atoms(self):
        bond = core.Bond()
        bond.left_bond_atoms = []
        bond.left_bond_atoms = core.AtomList()
        with self.assertRaises(ValueError):
            bond.left_bond_atoms = None

    def test_set_bond_bond_atoms(self):
        bond = core.Bond()
        bond.right_bond_atoms = []
        bond.right_bond_atoms = core.AtomList()
        with self.assertRaises(ValueError):
            bond.right_bond_atoms = None

    def test_set_left_displaced_atoms(self):
        bond = core.Bond()
        bond.left_displaced_atoms = []
        bond.left_displaced_atoms = core.AtomList()
        with self.assertRaises(ValueError):
            bond.left_displaced_atoms = None

    def test_set_bond_displaced_atoms(self):
        bond = core.Bond()
        bond.right_displaced_atoms = []
        bond.right_displaced_atoms = core.AtomList()
        with self.assertRaises(ValueError):
            bond.right_displaced_atoms = None

    def test_get_formula(self):
        bond = core.Bond()

        bond.left_displaced_atoms = core.AtomList([core.Atom(core.Monomer, 'C')])
        bond.right_displaced_atoms = core.AtomList([core.Atom(core.Monomer, 'H'), core.Atom(core.Monomer, 'H')])
        self.assertEqual(bond.get_formula(), EmpiricalFormula('CH2') * -1)

    def test_get_mol_wt(self):
        bond = core.Bond()

        self.assertEqual(bond.get_mol_wt(), 0.)

    def test_get_charge(self):
        bond = core.Bond()

        bond.left_displaced_atoms = core.AtomList([core.Atom(core.Monomer, 'C', charge=2)])
        bond.right_displaced_atoms = core.AtomList([core.Atom(core.Monomer, 'H', charge=3)])
        self.assertEqual(bond.get_charge(), -5)

    def test_is_equal(self):
        bond_1 = core.Bond()
        bond_2 = core.Bond()
        bond_3 = core.Bond(left_bond_atoms=[core.Atom(core.Monomer, 'H')])
        bond_4 = core.Bond(right_bond_atoms=[core.Atom(core.Monomer, 'H')])
        bond_5 = core.Bond(left_displaced_atoms=[core.Atom(core.Monomer, 'H')])
        bond_6 = core.Bond(right_displaced_atoms=[core.Atom(core.Monomer, 'H')])
        self.assertTrue(bond_1.is_equal(bond_1))
        self.assertTrue(bond_1.is_equal(bond_2))
        self.assertFalse(bond_1.is_equal({}))
        self.assertFalse(bond_1.is_equal(bond_3))
        self.assertFalse(bond_1.is_equal(bond_4))
        self.assertFalse(bond_1.is_equal(bond_5))
        self.assertFalse(bond_1.is_equal(bond_6))


class BpFormTestCase(unittest.TestCase):
    def test_init(self):
        bp_form = core.BpForm()
        self.assertEqual(bp_form.monomer_seq, core.MonomerSequence())
        self.assertEqual(bp_form.alphabet.monomers, {})
        self.assertEqual(bp_form.backbone.get_formula(), EmpiricalFormula())
        self.assertEqual(bp_form.backbone.get_charge(), 0)
        self.assertEqual(bp_form.bond.get_formula(), EmpiricalFormula())
        self.assertEqual(bp_form.bond.get_charge(), 0)

    def test_set_monomer_seq(self):
        bp_form = core.BpForm()

        bp_form.monomer_seq = core.MonomerSequence()
        self.assertEqual(len(bp_form.monomer_seq), 0)

        bp_form.monomer_seq = [core.Monomer(), core.Monomer()]
        self.assertIsInstance(bp_form.monomer_seq, core.MonomerSequence)
        self.assertEqual(len(bp_form.monomer_seq), 2)

        with self.assertRaises(ValueError):
            bp_form.monomer_seq = None
        with self.assertRaises(ValueError):
            bp_form.monomer_seq = 'A'

    def test_set_alphabet(self):
        bp_form = core.BpForm()

        bp_form.alphabet = dna.canonical_dna_alphabet
        self.assertEqual(len(bp_form.alphabet.monomers), 6)

        with self.assertRaises(ValueError):
            bp_form.alphabet = None
        with self.assertRaises(ValueError):
            bp_form.alphabet = 'A'

    def test_set_backbone(self):
        bp_form = core.BpForm()

        bp_form.backbone = core.Backbone()

        with self.assertRaises(ValueError):
            bp_form.backbone = None
        with self.assertRaises(ValueError):
            bp_form.backbone = '123'

    def test_set_bond(self):
        bp_form = core.BpForm()

        bp_form.bond = core.Bond()

        with self.assertRaises(ValueError):
            bp_form.bond = None
        with self.assertRaises(ValueError):
            bp_form.bond = '123'

    def test_set_circular(self):
        bp_form = core.BpForm()

        bp_form.circular = True
        self.assertEqual(bp_form.circular, True)

        bp_form.circular = False
        self.assertEqual(bp_form.circular, False)

        with self.assertRaises(ValueError):
            bp_form.circular = None

    def test_is_equal(self):
        bp_form_1 = core.BpForm(monomer_seq=core.MonomerSequence(
            [core.Monomer(id='A'), core.Monomer(id='B')]))
        bp_form_2 = core.BpForm(monomer_seq=core.MonomerSequence(
            [core.Monomer(id='A'), core.Monomer(id='B')]))
        bp_form_3 = None
        bp_form_4 = core.BpForm(monomer_seq=core.MonomerSequence(
            [core.Monomer(id='A'), core.Monomer(id='B')]), alphabet=dna.canonical_dna_alphabet)
        bp_form_5 = core.BpForm(monomer_seq=core.MonomerSequence(
            [core.Monomer(id='A'), core.Monomer(id='B')]), backbone=core.Backbone(structure='O'))
        bp_form_6 = core.BpForm(monomer_seq=core.MonomerSequence(
            [core.Monomer(id='A'), core.Monomer(id='B')]), bond=core.Bond(left_bond_atoms=[core.Atom(core.Monomer, 'C')]))
        bp_form_7 = core.BpForm(monomer_seq=core.MonomerSequence(
            [core.Monomer(id='A'), core.Monomer(id='B')]), circular=True)
        self.assertTrue(bp_form_1.is_equal(bp_form_1))
        self.assertTrue(bp_form_1.is_equal(bp_form_2))
        self.assertFalse(bp_form_1.is_equal(bp_form_3))
        self.assertFalse(bp_form_1.is_equal(bp_form_4))
        self.assertFalse(bp_form_1.is_equal(bp_form_5))
        self.assertFalse(bp_form_1.is_equal(bp_form_6))
        self.assertFalse(bp_form_1.is_equal(bp_form_7))

    def test_getitem(self):
        monomer_1 = core.Monomer(id='A')
        monomer_2 = core.Monomer(id='B')
        monomer_3 = core.Monomer(id='C')
        bp_form = core.BpForm([monomer_1, monomer_2, monomer_3])
        self.assertEqual(bp_form[0], monomer_1)
        self.assertEqual(bp_form[1], monomer_2)
        self.assertEqual(bp_form[0:1], [monomer_1])

    def test_setitem(self):
        monomer_1 = core.Monomer(id='A')
        monomer_2 = core.Monomer(id='B')
        monomer_3 = core.Monomer(id='C')
        bp_form = core.BpForm([monomer_1, monomer_2, monomer_3])

        self.assertEqual(bp_form[0], monomer_1)

        bp_form[0] = monomer_2
        self.assertEqual(bp_form[0], monomer_2)

        bp_form[0:1] = [monomer_3]
        self.assertEqual(bp_form[0], monomer_3)

    def test_delitem(self):
        monomer_1 = core.Monomer(id='A')
        monomer_2 = core.Monomer(id='B')
        monomer_3 = core.Monomer(id='C')
        bp_form = core.BpForm([monomer_1, monomer_2, monomer_3])
        del(bp_form[1])
        self.assertTrue(bp_form.is_equal(core.BpForm([monomer_1, monomer_3])))

    def test_iter(self):
        monomer_1 = core.Monomer(id='A')
        monomer_2 = core.Monomer(id='B')
        monomer_3 = core.Monomer(id='C')
        bp_form = core.BpForm([monomer_1, monomer_2, monomer_3])
        for i_monomer, monomer in enumerate(bp_form):
            if i_monomer == 0:
                self.assertEqual(monomer, monomer_1)
            if i_monomer == 1:
                self.assertEqual(monomer, monomer_2)
            if i_monomer == 2:
                self.assertEqual(monomer, monomer_3)

    def test_reversed(self):
        monomer_1 = core.Monomer(id='A')
        monomer_2 = core.Monomer(id='B')
        monomer_3 = core.Monomer(id='C')
        bp_form = core.BpForm([monomer_1, monomer_2, monomer_3])
        for i_monomer, monomer in enumerate(reversed(bp_form)):
            if i_monomer == 2:
                self.assertEqual(monomer, monomer_1)
            if i_monomer == 1:
                self.assertEqual(monomer, monomer_2)
            if i_monomer == 0:
                self.assertEqual(monomer, monomer_3)

    def test_contains(self):
        monomer_1 = core.Monomer(id='A')
        monomer_2 = core.Monomer(id='B')
        monomer_3 = core.Monomer(id='C')
        bp_form = core.BpForm([monomer_1, monomer_2])
        self.assertIn(monomer_1, bp_form)
        self.assertIn(monomer_2, bp_form)
        self.assertNotIn(monomer_3, bp_form)

    def test_len(self):
        bp_form = core.BpForm()
        self.assertEqual(len(bp_form), 0)

        bp_form = core.BpForm(monomer_seq=[core.Monomer(), core.Monomer()])
        self.assertEqual(len(bp_form), 2)

    def test_get_monomer_counts(self):
        monomer_1 = core.Monomer(id='A')
        monomer_2 = core.Monomer(id='B')
        monomer_3 = core.Monomer(id='C')
        bp_form = core.BpForm([monomer_1, monomer_2, monomer_1, monomer_1, monomer_1, monomer_2, monomer_2, monomer_3])
        self.assertEqual(bp_form.get_monomer_counts(), {
            monomer_1: 4,
            monomer_2: 3,
            monomer_3: 1,
        })

    def test_get_formula_mol_wt_charge(self):
        monomer_A = core.Monomer(id='A', structure=dAMP_smiles)
        monomer_C = core.Monomer(id='C', structure=dCMP_smiles)

        bp_form = core.BpForm([monomer_A])
        self.assertEqual(bp_form.get_formula(), monomer_A.get_formula())
        self.assertEqual(bp_form.get_mol_wt(), monomer_A.get_mol_wt())
        self.assertEqual(bp_form.get_charge(), monomer_A.get_charge())

        bp_form = core.BpForm([monomer_C])
        self.assertEqual(bp_form.get_formula(), monomer_C.get_formula())
        self.assertEqual(bp_form.get_mol_wt(), monomer_C.get_mol_wt())
        self.assertEqual(bp_form.get_charge(), monomer_C.get_charge())

        bp_form = core.BpForm([monomer_A, monomer_C])
        self.assertEqual(bp_form.get_formula(), monomer_A.get_formula() + monomer_C.get_formula())
        self.assertEqual(bp_form.get_mol_wt(), monomer_A.get_mol_wt() + monomer_C.get_mol_wt())
        self.assertEqual(bp_form.get_charge(), monomer_A.get_charge() + monomer_C.get_charge())

        bp_form = core.BpForm([monomer_A, monomer_C],
                              bond=core.Bond(right_displaced_atoms=[core.Atom(core.Monomer, 'H', charge=-1, position=None)]))
        self.assertEqual(bp_form.get_formula(), monomer_A.get_formula() + monomer_C.get_formula() - EmpiricalFormula('H'))
        self.assertEqual(bp_form.get_mol_wt(), monomer_A.get_mol_wt() + monomer_C.get_mol_wt() -
                         EmpiricalFormula('H').get_molecular_weight())
        self.assertEqual(bp_form.get_charge(), monomer_A.get_charge() + monomer_C.get_charge() + 1)

        bp_form = core.BpForm([monomer_A, monomer_A, monomer_C, monomer_C, monomer_C],
                              bond=core.Bond(right_displaced_atoms=[core.Atom(core.Monomer, 'H', charge=-1, position=None)]))
        self.assertEqual(bp_form.get_formula(), monomer_A.get_formula() * 2 + monomer_C.get_formula() * 3 - EmpiricalFormula('H') * 4)
        self.assertEqual(bp_form.get_mol_wt(), monomer_A.get_mol_wt() * 2 + monomer_C.get_mol_wt()
                         * 3 - EmpiricalFormula('H').get_molecular_weight() * 4)
        self.assertEqual(bp_form.get_charge(), monomer_A.get_charge() * 2 + monomer_C.get_charge() * 3 + 1 * 4)

    def test_get_major_micro_species(self):
        monomer_1 = core.Monomer(id='A', structure=dAMP_smiles)
        monomer_2 = core.Monomer(id='C', structure=dCMP_smiles)
        bp_form = core.BpForm([monomer_1, monomer_2])
        bp_form.get_major_micro_species(7.)

    def test_str(self):
        monomer_A = core.Monomer(id='A', structure=dAMP_smiles)
        monomer_C = core.Monomer(id='C', structure=dCMP_smiles)
        monomer_G = core.Monomer(id='G', structure=dGMP_smiles)

        bp_form = core.BpForm([monomer_A, monomer_C, monomer_G, monomer_A])
        self.assertEqual(str(bp_form), '{}{}{}{}'.format(str(monomer_A), str(monomer_C), str(monomer_G), str(monomer_A)))

        bp_form = core.BpForm([monomer_A, monomer_C, monomer_G, monomer_A], alphabet=core.Alphabet(monomers={
            'A': monomer_A,
            'C': monomer_C,
        }))
        self.assertEqual(str(bp_form), '{}{}{}{}'.format('A', 'C', str(monomer_G), 'A'))
        dGMP_smiles_2 = 'Nc1[nH]c(=O)c2c(n1)[nH]cn2'
        self.assertEqual(str(bp_form), '{}{}{}{}'.format('A', 'C', '[id: "{}" | structure: "{}"]'.format('G', dGMP_smiles_2), 'A'))

    def test_from_str(self):
        self.assertTrue(dna.DnaForm().from_str('AAA').is_equal(dna.DnaForm([
            dna.canonical_dna_alphabet.monomers.A, dna.canonical_dna_alphabet.monomers.A, dna.canonical_dna_alphabet.monomers.A,
        ])))

        self.assertTrue(dna.DnaForm().from_str('ACTG').is_equal(dna.DnaForm([
            dna.canonical_dna_alphabet.monomers.A, dna.canonical_dna_alphabet.monomers.C,
            dna.canonical_dna_alphabet.monomers.T, dna.canonical_dna_alphabet.monomers.G,
        ])))

        with self.assertRaisesRegex(lark.exceptions.VisitError, 'not in alphabet'):
            self.assertTrue(dna.DnaForm().from_str('EAA').is_equal(dna.DnaForm([
                dna.canonical_dna_alphabet.monomers.A, dna.canonical_dna_alphabet.monomers.A,
                dna.canonical_dna_alphabet.monomers.A,
            ])))

        dna_form_1 = ('AA[id: "dI"'
                      + ' | name: "2\'-deoxyinosine"'
                      + ' | synonym: "2\'-deoxyinosine, 9-[(2R,4S,5R)-4-hydroxy-5-(hydroxymethyl)tetrahydrofuran-2-yl]-9H-purin-6-ol"'
                      + ' | identifier: "chebi" / "CHEBI:28997"'
                      + ' | structure: "' + dIMP_smiles + '"'
                      + ' | monomer-bond-atom: Monomer / C / 1 / -1'
                      + ' | monomer-bond-atom: Monomer / D / 2 / -2'
                      + ' | monomer-displaced-atom: Backbone / D / 2 / -2'
                      + ' | left-bond-atom: Monomer / E / 3 / -3'
                      + ' | left-displaced-atom: Backbone / F / 4 / -4'
                      + ' | right-bond-atom: Monomer / G / 5 / -5'
                      + ' | right-displaced-atom: Backbone / H / 6 / -6'
                      + ' | delta-mass: -2.5'
                      + ' | delta-charge: 3'
                      + ' | position: 3-5'
                      + ' | base-monomer: "A"'
                      + ' | comments: "A purine 2\'-deoxyribonucleoside that is inosine ..."]A')
        dna_form_2 = dna.DnaForm().from_str(dna_form_1)
        dna_form_3 = dna.DnaForm([
            dna.canonical_dna_alphabet.monomers.A,
            dna.canonical_dna_alphabet.monomers.A,
            core.Monomer(
                id='dI',
                name="2'-deoxyinosine",
                synonyms=core.SynonymSet(
                    ["2'-deoxyinosine, 9-[(2R,4S,5R)-4-hydroxy-5-(hydroxymethyl)tetrahydrofuran-2-yl]-9H-purin-6-ol"]),
                identifiers=core.IdentifierSet([core.Identifier('chebi', 'CHEBI:28997')]),
                structure=dIMP_smiles,
                monomer_bond_atoms=[
                    core.Atom(core.Monomer, 'C', position=1, charge=-1),
                    core.Atom(core.Monomer, 'D', position=2, charge=-2),
                ],
                monomer_displaced_atoms=[core.Atom(core.Backbone, 'D', position=2, charge=-2)],
                left_bond_atoms=[core.Atom(core.Monomer, 'E', position=3, charge=-3)],
                left_displaced_atoms=[core.Atom(core.Backbone, 'F', position=4, charge=-4)],
                right_bond_atoms=[core.Atom(core.Monomer, 'G', position=5, charge=-5)],
                right_displaced_atoms=[core.Atom(core.Backbone, 'H', position=6, charge=-6)],
                delta_mass=-2.5,
                delta_charge=3,
                start_position=3,
                end_position=5,
                base_monomers=[dna.canonical_dna_alphabet.monomers.A],
                comments="A purine 2'-deoxyribonucleoside that is inosine ...",
            ),
            dna.canonical_dna_alphabet.monomers.A,
        ])
        self.assertEqual(str(dna_form_2), dna_form_1.replace(dIMP_smiles, 'OCC1OC(CC1O)n1cnc2c1ncnc2O'))
        self.assertTrue(dna_form_2.is_equal(dna_form_3))

        self.assertTrue(dna.DnaForm().from_str(
            'AA[id: "dI"'
            ' | position: 3-]A').is_equal(dna.DnaForm([
                dna.canonical_dna_alphabet.monomers.A,
                dna.canonical_dna_alphabet.monomers.A,
                core.Monomer(
                    id='dI',
                    start_position=3,
                    end_position=None,
                ),
                dna.canonical_dna_alphabet.monomers.A,
            ])))

        self.assertTrue(dna.DnaForm().from_str(
            'AA[id: "dI"'
            ' | position: -5]A').is_equal(dna.DnaForm([
                dna.canonical_dna_alphabet.monomers.A,
                dna.canonical_dna_alphabet.monomers.A,
                core.Monomer(
                    id='dI',
                    start_position=None,
                    end_position=5,
                ),
                dna.canonical_dna_alphabet.monomers.A,
            ])))

        alphabet = core.Alphabet()
        alphabet.monomers['aA'] = core.Monomer()
        alphabet.monomers['Aa'] = core.Monomer()
        alphabet.monomers['A'] = core.Monomer()
        alphabet.monomers['AA'] = core.Monomer()
        alphabet.monomers['*'] = core.Monomer()
        alphabet.monomers[' '] = core.Monomer()
        alphabet.monomers[' A'] = core.Monomer()
        self.assertTrue(core.BpForm(alphabet=alphabet).from_str(
            'AAA{AA}AA{aA}{Aa}AA* { A}').is_equal(core.BpForm([
                alphabet.monomers['A'], alphabet.monomers['A'], alphabet.monomers['A'],
                alphabet.monomers['AA'],
                alphabet.monomers['A'], alphabet.monomers['A'],
                alphabet.monomers['aA'],
                alphabet.monomers['Aa'],
                alphabet.monomers['A'], alphabet.monomers['A'],
                alphabet.monomers['*'],
                alphabet.monomers[' '],
                alphabet.monomers[' A'],
            ], alphabet=alphabet)))

        as_str = 'AAA{AA}AA{aA}{Aa}AA'
        form = core.BpForm(alphabet=alphabet).from_str(as_str)
        self.assertEqual(str(form), as_str)

        alphabet = core.Alphabet()
        alphabet.monomers['aA'] = core.Monomer()
        alphabet.monomers['Aa'] = core.Monomer()
        alphabet.monomers['A'] = core.Monomer()
        alphabet.monomers['AA'] = core.Monomer()
        with self.assertRaises(lark.exceptions.VisitError):
            core.BpForm(alphabet=alphabet).from_str('AAA{(AA}AA{aA}{Aa}AA')
        with self.assertRaises(lark.exceptions.UnexpectedCharacters):
            core.BpForm(alphabet=alphabet).from_str('AAA{AA}AA{aA}{[Aa}AA')
        with self.assertRaises(lark.exceptions.VisitError):
            core.BpForm(alphabet=alphabet).from_str('AAA[base-monomer: "C"]')

        with self.assertRaisesRegex(lark.exceptions.VisitError, 'cannot be repeated'):
            dna.DnaForm().from_str(
                'AA[id: "dI"'
                ' | name: "2\'-deoxyinosine"'
                ' | name: "2\'-deoxyinosine"]A')

        dna_form_1 = dna.DnaForm().from_str('[structure: "' + dIMP_smiles + '"]')
        dna_form_2 = dna.DnaForm().from_str('[structure: "' + dIMP_smiles + '"]')
        self.assertTrue(dna_form_1.is_equal(dna_form_2))

    def test_export(self):
        form = dna.CanonicalDnaForm().from_str('A')
        self.assertEqual(form.export('inchi'), ('InChI=1S'
                                                '/C10H14N5O6P'
                                                '/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(16)6(21-7)2-20-22(17,18)19'
                                                '/h3-7,16H,1-2H2,(H2,11,12,13)(H2,17,18,19)'
                                                '/p-2'))

        form = dna.CanonicalDnaForm().from_str('AA')
        self.assertEqual(form.export('inchi'), ('InChI=1S'
                                                '/C20H26N10O11P2'
                                                '/c21-17-15-19(25-5-23-17)29(7-27-15)13-1-9(31)11(39-13)3-38-43(35,36)'
                                                '41-10-2-14(40-12(10)4-37-42(32,33)34)30-8-28-16-18(22)24-6-26-20(16)30'
                                                '/h5-14,31H,1-4H2,(H,35,36)(H2,21,23,25)(H2,22,24,26)(H2,32,33,34)'
                                                '/p-3'))

        form = rna.CanonicalRnaForm().from_str('AA')
        self.assertEqual(form.export('inchi'), ('InChI=1S'
                                                '/C20H26N10O13P2'
                                                '/c21-15-9-17(25-3-23-15)29(5-27-9)19-12(32)11(31)7(41-19)1-40-45'
                                                '(37,38)43-14-8(2-39-44(34,35)36)'
                                                '42-20(13(14)33)30-6-28-10-16(22)24-4-26-18(10)30'
                                                '/h3-8,11-14,19-20,31-33H,1-2H2,(H,37,38)(H2,21,23,25)(H2,22,24,26)(H2,34,35,36)'
                                                '/p-3'))

        form = protein.CanonicalProteinForm().from_str('AA')
        self.assertEqual(form.export('inchi'), 'InChI=1S/C6H12N2O3/c1-3(7)5(9)8-4(2)6(10)11/h3-4H,7H2,1-2H3,(H,8,9)(H,10,11)/p+1')

        form = dna.CanonicalDnaForm()
        self.assertEqual(form.export('inchi'), None)

        form = dna.CanonicalDnaForm().from_str('A[structure: "NC1=C2N=CNC2=NC=N1"]')
        with self.assertRaises(ValueError):
            form.export('inchi')

    def test_circular_export(self):
        form = dna.CanonicalDnaForm(circular=True).from_str('A')
        self.assertEqual(form.export('inchi'), ('InChI=1S/C10H12N5O5P'
                                                '/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5-6(19-7)2-18-21(16,17)20-5'
                                                '/h3-7H,1-2H2,(H,16,17)(H2,11,12,13)'
                                                '/p-1'))

        form = dna.CanonicalDnaForm(circular=True).from_str('AA')
        self.assertEqual(form.export('inchi'), ('InChI=1S/C20H24N10O10P2'
                                                '/c21-17-15-19(25-5-23-17)29(7-27-15)13-1-9-11(37-13)3-35-42(33,34)'
                                                '40-10-2-14(38-12(10)4-36-41(31,32)39-9)30-8-28-16-18(22)24-6-26-20(16)30'
                                                '/h5-14H,1-4H2,(H,31,32)(H,33,34)(H2,21,23,25)(H2,22,24,26)/p-2'))

        form = rna.CanonicalRnaForm(circular=True).from_str('AA')
        self.assertEqual(form.export('inchi'), ('InChI=1S/C20H24N10O12P2'
                                                '/c21-15-9-17(25-3-23-15)29(5-27-9)19-11(31)13-7(39-19)1-37-43(33,34)42-14-8(2-38-44(35,36)41-13)40-20'
                                                '(12(14)32)30-6-28-10-16(22)24-4-26-18(10)30'
                                                '/h3-8,11-14,19-20,31-32H,1-2H2,(H,33,34)(H,35,36)(H2,21,23,25)(H2,22,24,26)/p-2'))

    def test_get_fasta(self):
        alphabet = core.Alphabet()
        alphabet.monomers.A = core.Monomer()
        alphabet.monomers.C = core.Monomer()
        alphabet.monomers.G = core.Monomer()
        alphabet.monomers.T = core.Monomer()
        alphabet.monomers.m2A = core.Monomer(base_monomers=[alphabet.monomers.A])
        alphabet.monomers.m22A = core.Monomer(base_monomers=[alphabet.monomers.m2A])
        alphabet.monomers.m222A = core.Monomer(base_monomers=[alphabet.monomers.m22A])
        alphabet.monomers.m2222A = core.Monomer(base_monomers=[alphabet.monomers.A, alphabet.monomers.m222A])
        alphabet.monomers.m2222C = core.Monomer(base_monomers=[alphabet.monomers.C, alphabet.monomers.m222A])

        self.assertEqual(alphabet.monomers.A.get_root_monomers(), set([alphabet.monomers.A]))
        self.assertEqual(alphabet.monomers.C.get_root_monomers(), set([alphabet.monomers.C]))
        self.assertEqual(alphabet.monomers.G.get_root_monomers(), set([alphabet.monomers.G]))
        self.assertEqual(alphabet.monomers.T.get_root_monomers(), set([alphabet.monomers.T]))
        self.assertEqual(alphabet.monomers.m2A.get_root_monomers(), set([alphabet.monomers.A]))
        self.assertEqual(alphabet.monomers.m22A.get_root_monomers(), set([alphabet.monomers.A]))
        self.assertEqual(alphabet.monomers.m2222A.get_root_monomers(), set([alphabet.monomers.A]))
        self.assertEqual(alphabet.monomers.m2222C.get_root_monomers(), set([alphabet.monomers.A, alphabet.monomers.C]))

        bpform = core.BpForm(alphabet=alphabet, monomer_seq=[
            alphabet.monomers.A, alphabet.monomers.C, alphabet.monomers.G, alphabet.monomers.T,
            alphabet.monomers.m2A, alphabet.monomers.m22A, alphabet.monomers.m222A,
            alphabet.monomers.m2222A, alphabet.monomers.m2222C,
        ])
        self.assertEqual(bpform.get_fasta(), 'ACGTAAAA?')

        class CustomBpForm(core.BpForm):
            DEFAULT_FASTA_CODE = 'X'
        bpform = CustomBpForm(alphabet=alphabet, monomer_seq=[
            alphabet.monomers.A, alphabet.monomers.C, alphabet.monomers.G, alphabet.monomers.T,
            alphabet.monomers.m2A, alphabet.monomers.m22A, alphabet.monomers.m222A,
            alphabet.monomers.m2222A, alphabet.monomers.m2222C,
        ])
        self.assertEqual(bpform.get_fasta(), 'ACGTAAAAX')


class AlphabetTestCase(unittest.TestCase):
    def setUp(self):
        self.dir_path = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir_path)

    def test_set_monomers(self):
        alphabet = core.Alphabet()
        alphabet.monomers = core.MonomerDict()
        alphabet.monomers = {}
        with self.assertRaises(ValueError):
            alphabet.monomers = None

    def test_getitem(self):
        self.assertEqual(dna.canonical_dna_alphabet.monomers.A.export('inchi'), dAMP_inchi)

    def test_setitem(self):
        alphabet = core.Alphabet()
        alphabet.monomers['A'] = core.Monomer(structure=dAMP_smiles)
        alphabet.monomers['aA'] = core.Monomer(structure=dAMP_smiles)
        alphabet.monomers['Aa'] = core.Monomer(structure=dAMP_smiles)
        alphabet.monomers['*'] = core.Monomer(structure=dAMP_smiles)
        with self.assertRaises(ValueError):
            alphabet.monomers['{aa'] = core.Monomer(structure=dAMP_smiles)
        with self.assertRaises(ValueError):
            alphabet.monomers['A]'] = core.Monomer(structure=dAMP_smiles)

    def test_get_monomer_code(self):
        alphabet = dna.canonical_dna_alphabet
        self.assertEqual(alphabet.get_monomer_code(alphabet.monomers.A), 'A')
        with self.assertRaises(ValueError):
            alphabet.get_monomer_code(core.Monomer())

    def test_get_major_micro_species(self):
        alphabet = core.Alphabet(monomers={
            'A': core.Monomer(structure=dAMP_smiles),
            'C': core.Monomer(structure=dCMP_smiles),
        })
        alphabet.get_major_micro_species(7.)

    def test_is_equal(self):
        self.assertTrue(dna.canonical_dna_alphabet.is_equal(dna.canonical_dna_alphabet))
        self.assertFalse(dna.canonical_dna_alphabet.is_equal(rna.rna_alphabet))
        self.assertFalse(dna.canonical_dna_alphabet.is_equal(None))
        self.assertFalse(dna.canonical_dna_alphabet.is_equal(core.Alphabet()))

        dna_alphabet_1 = core.Alphabet(id=dna.canonical_dna_alphabet.id,
                                       name=dna.canonical_dna_alphabet.name,
                                       description=dna.canonical_dna_alphabet.description,
                                       monomers={
                                           'A': dna.canonical_dna_alphabet.monomers.A,
                                           'C': dna.canonical_dna_alphabet.monomers.C,
                                           'G': dna.canonical_dna_alphabet.monomers.G,
                                           'T': dna.canonical_dna_alphabet.monomers.T,
                                       })
        dna_alphabet_2 = core.Alphabet(id=dna.canonical_dna_alphabet.id,
                                       name=dna.canonical_dna_alphabet.name,
                                       description=dna.canonical_dna_alphabet.description,
                                       monomers={
                                           'A': dna.canonical_dna_alphabet.monomers.A,
                                           'C': dna.canonical_dna_alphabet.monomers.C,
                                           'G': dna.canonical_dna_alphabet.monomers.G,
                                           'T': dna.canonical_dna_alphabet.monomers.T,
                                       })
        self.assertTrue(dna_alphabet_1.is_equal(dna_alphabet_2))

        dna_alphabet_1 = core.Alphabet(id=dna.canonical_dna_alphabet.id,
                                       name=dna.canonical_dna_alphabet.name,
                                       description=dna.canonical_dna_alphabet.description,
                                       monomers={
                                           'A': dna.canonical_dna_alphabet.monomers.A,
                                           'C': dna.canonical_dna_alphabet.monomers.C,
                                           'G': dna.canonical_dna_alphabet.monomers.G,
                                           'T': dna.canonical_dna_alphabet.monomers.T,
                                       })
        dna_alphabet_2 = core.Alphabet(id=dna.canonical_dna_alphabet.id,
                                       name=dna.canonical_dna_alphabet.name,
                                       description=dna.canonical_dna_alphabet.description,
                                       monomers={
                                           'A': dna.canonical_dna_alphabet.monomers.A,
                                           'C': dna.canonical_dna_alphabet.monomers.C,
                                           'G': dna.canonical_dna_alphabet.monomers.G,
                                       })
        self.assertFalse(dna_alphabet_1.is_equal(dna_alphabet_2))

        dna_alphabet_1 = core.Alphabet(id=dna.canonical_dna_alphabet.id,
                                       name=dna.canonical_dna_alphabet.name,
                                       description=dna.canonical_dna_alphabet.description,
                                       monomers={
                                           'A': dna.canonical_dna_alphabet.monomers.A,
                                           'C': dna.canonical_dna_alphabet.monomers.C,
                                           'G': dna.canonical_dna_alphabet.monomers.G,
                                           'T': dna.canonical_dna_alphabet.monomers.T,
                                       })
        dna_alphabet_2 = core.Alphabet(id=dna.canonical_dna_alphabet.id,
                                       name=dna.canonical_dna_alphabet.name,
                                       description=dna.canonical_dna_alphabet.description,
                                       monomers={
                                           'C': dna.canonical_dna_alphabet.monomers.A,
                                           'A': dna.canonical_dna_alphabet.monomers.C,
                                           'G': dna.canonical_dna_alphabet.monomers.G,
                                           'T': dna.canonical_dna_alphabet.monomers.T,
                                       })
        self.assertFalse(dna_alphabet_1.is_equal(dna_alphabet_2))

    def test_to_from_yaml(self):
        dna_alphabet = core.Alphabet(id='test',
                                     name='Test',
                                     description='Test description',
                                     monomers={
                                         'A': dna.canonical_dna_alphabet.monomers.A,
                                         'C': dna.canonical_dna_alphabet.monomers.C,
                                         'G': dna.canonical_dna_alphabet.monomers.G,
                                         'T': dna.canonical_dna_alphabet.monomers.T,
                                     })
        path = os.path.join(self.dir_path, 'alphabet.yml')
        dna_alphabet.to_yaml(path)
        dna_alphabet_2 = core.Alphabet().from_yaml(path)
        self.assertTrue(dna_alphabet_2.is_equal(dna_alphabet))

    def test_to_from_yaml_without_name(self):
        dna_alphabet = core.Alphabet(id='test',
                                     monomers={
                                         'A': dna.canonical_dna_alphabet.monomers.A,
                                         'C': dna.canonical_dna_alphabet.monomers.C,
                                         'G': dna.canonical_dna_alphabet.monomers.G,
                                         'T': dna.canonical_dna_alphabet.monomers.T,
                                     })
        path = os.path.join(self.dir_path, 'alphabet.yml')
        dna_alphabet.to_yaml(path)
        dna_alphabet_2 = core.Alphabet().from_yaml(path)
        self.assertTrue(dna_alphabet_2.is_equal(dna_alphabet))


class AlphabetBuilderTestCase(unittest.TestCase):
    def setUp(self):
        self.dir_path = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir_path)

    def test(self):
        class AlphabetBuilder(core.AlphabetBuilder):
            def build(self):
                return core.Alphabet()
        alphabet = AlphabetBuilder().run()
        self.assertIsInstance(alphabet, core.Alphabet)

        path = os.path.join(self.dir_path, 'alphabet.yml')
        alphabet = AlphabetBuilder().run(ph=7.4, major_tautomer=True, path=path)
        self.assertTrue(os.path.isfile(path))
