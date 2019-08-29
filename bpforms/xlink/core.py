""" Classes to represent an ontology of crosslinks

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-08-12
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Atom, Bond, BondOrder, BondStereo, Monomer, parse_yaml
from bpforms.util import get_alphabet
from ruamel import yaml
import importlib
import os
import pkg_resources
import re


class CrosslinksOnto(dict):
    pass


onto_filename = pkg_resources.resource_filename('bpforms', os.path.join('xlink', 'xlink.yml'))


def load_onto(filename=onto_filename):
    """ Load an ontology of crosslinks from a YAML file

    Args:
        filename (:obj:`str`): path to ontology

    Returns:
        :obj:`CrosslinksOnto`: ontology of crosslinks
    """
    onto_dict = parse_yaml(filename)
    onto = CrosslinksOnto()
    for id, bond_dict in onto_dict.items():
        bond = Bond()
        onto[id] = bond

        bond.id = id
        bond.name = bond_dict.get('name', None)
        bond.synonyms = bond_dict.get('synonyms', [])        

        l_alph = get_alphabet(bond_dict['l_monomer_alphabet'])
        r_alph = get_alphabet(bond_dict['r_monomer_alphabet'])
        bond.l_monomer = l_alph.monomers.get(bond_dict['l_monomer'])
        bond.r_monomer = r_alph.monomers.get(bond_dict['r_monomer'])

        for atom_type in ['l_bond_atoms', 'r_bond_atoms', 'l_displaced_atoms', 'r_displaced_atoms']:
            for atom in bond_dict[atom_type]:
                element, position, charge = parse_atom(atom)
                getattr(bond, atom_type).append(Atom(Monomer, element, position=position, charge=charge))

        bond.order = BondOrder[bond_dict.get('order', 'single')]
        stereo = bond_dict.get('stereo', None)
        if stereo is None:
            bond.stereo = None
        else:
            bond.stereo = BondStereo[stereo]

        bond.comments = bond_dict.get('comments', None)

    crosslink_to_id = {xlink: id for id, xlink in onto.items()}

    return onto, crosslink_to_id


def parse_atom(atom):
    """ Parse description of atom in ontology

    Args:
        atom (:obj:`str`)

    Returns:
        :obj:`tuple`: element, position, charge
    """
    match = re.match(r'([A-Z][a-z]?)(\d+)(([\+\-]\d+)?)', atom)
    return (match.group(1), int(float(match.group(2))), int(float(match.group(3) or 0)))


crosslinks_onto, onto_crosslink_to_id = load_onto()
