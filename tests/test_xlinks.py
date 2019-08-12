""" Test of bpforms.xlink

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-08-12
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.alphabet import protein
from bpforms.core import get_hydrogen_atom
from bpforms.xlink.core import crosslinks_onto, parse_atom
import openbabel
import unittest


class XlinkTestCase(unittest.TestCase):
    def test_parse_atom(self):
        self.assertEqual(parse_atom('C1'), ('C', 1, 0))
        self.assertEqual(parse_atom('C1+1'), ('C', 1, 1))
        self.assertEqual(parse_atom('C1-1'), ('C', 1, -1))
        self.assertEqual(parse_atom('Ca3+2'), ('Ca', 3, 2))

    def test_bonds(self):
        el_table = openbabel.OBElementTable()

        for id, bond in crosslinks_onto.items():
            for direction in ['l', 'r']:
                monomer = getattr(bond, direction + '_monomer')
                h_atoms = []
                for atom_type in ['bond_atoms', 'displaced_atoms']:
                    for i_atom, atom_md in enumerate(getattr(bond, direction + '_' + atom_type)):
                        atom = monomer.structure.GetAtom(atom_md.position)
                        if atom_md.element == 'H':
                            atom = get_hydrogen_atom(atom, h_atoms, 0)

                        if atom is not None:
                            el = el_table.GetSymbol(atom.GetAtomicNum())
                            self.assertEqual(el, atom_md.element, "{}_{} {} of {} doesn't match structure of {} {}".format(
                                direction, atom_type, i_atom, id, monomer.id, monomer.export('smiles')))

                        elif atom_md.element != 'H':
                            raise Exception("{}_{} {} of {} doesn't match structure of {} {}".format(
                                direction, atom_type, i_atom, id, monomer.id, monomer.export('smiles')))
