""" Utilities for BpForms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from . import core
from .alphabet import dna
from .alphabet import rna
from .alphabet import protein
import openbabel


def get_alphabets():
    """ Get a list of available alphabets 

    Returns:
        :obj:`dict`: dictionary which maps the ids of alphabets to alphabets
    """
    alphabets = [
        dna.dna_alphabet,
        rna.rna_alphabet,
        protein.protein_alphabet,
        dna.canonical_dna_alphabet,
        rna.canonical_rna_alphabet,
        protein.canonical_protein_alphabet,
    ]
    return {alphabet.id: alphabet for alphabet in alphabets}


def get_alphabet(alphabet):
    """ Get an alphabet

    Args:
        alphabet (:obj:`str`): alphabet

    Returns:
        :obj:`core.Alphabet`: alphabet
    """
    alphabet_obj = get_alphabets().get(alphabet, None)
    if alphabet_obj is None:
        raise ValueError('Unknown alphabet "{}"'.format(alphabet))
    return alphabet_obj


def get_form(alphabet):
    """ Get a subclass of BpFrom

    Args:
        alphabet (:obj:`str`): alphabet

    Returns:
        :obj:`type`: subclass of BpForm
    """
    if alphabet == 'dna':
        return dna.DnaForm
    if alphabet == 'canonical_dna':
        return dna.CanonicalDnaForm

    if alphabet == 'rna':
        return rna.RnaForm
    if alphabet == 'canonical_rna':
        return rna.CanonicalRnaForm

    if alphabet == 'protein':
        return protein.ProteinForm
    if alphabet == 'canonical_protein':
        return protein.CanonicalProteinForm

    raise ValueError('Alphabet "{}" must be "dna", "rna", or "protein"'.format(alphabet))


def build_alphabets(ph=None, major_tautomer=False, _max_monomers=float('inf')):
    """ Build DNA, RNA, and protein alphabets

    Args:
        ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomer
        major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
        _max_monomers (:obj:`float`, optional): maximum number of monomers to build; used for testing
    """
    dna.DnaAlphabetBuilder(_max_monomers=_max_monomers).run(ph=ph, major_tautomer=major_tautomer)
    rna.RnaAlphabetBuilder(_max_monomers=_max_monomers).run(ph=ph, major_tautomer=major_tautomer)
    protein.ProteinAlphabetBuilder(_max_monomers=_max_monomers).run(ph=ph, major_tautomer=major_tautomer)


def gen_html_viz_alphabet(alphabet, filename):
    """ Create and save an HTML document with images of the monomers in an alphabet

    Args:
        alphabet (:obj:`Alphabet`): alphabet
        filename (:obj:`str`): path to save HTML document with images of monomers
    """
    doc = ''
    doc += '<html>'
    doc += '  <body>'
    doc += '    <table>'
    doc += '      <thead>'
    doc += '        <tr>'
    doc += '          <th>Code</th>'
    doc += '          <th>Structure</th>'
    doc += '        </tr>'
    doc += '      </thead>'
    for code, monomer in alphabet.monomers.items():
        doc += '        <tr>'
        doc += '          <td>{}</td>'.format(code)
        doc += '          <td>{}</td>'.format(
            monomer.get_image(include_xml_header=False))
        doc += '        </tr>'
    doc += '    </table>'
    doc += '  </body>'
    doc += '</html>'

    with open(filename, 'w') as file:
        file.write(doc)


def validate_bpform_linkages(form_type):
    """ Validate linkages in alphabet

    Args:
        form_type (:obj:`type`): type of BpForm

    Raises:
        :obj:`ValueError`: if any of the linkages are invalid
    """

    form = form_type()

    element_table = openbabel.OBElementTable()

    errors = []

    atom_types = [
        ['backbone', 'monomer_bond_atoms'],
        ['backbone', 'backbone_bond_atoms'],
        ['backbone', 'monomer_displaced_atoms'],
        ['backbone', 'backbone_displaced_atoms'],
        ['bond', 'left_bond_atoms'],
        ['bond', 'right_bond_atoms'],
        ['bond', 'left_displaced_atoms'],
        ['bond', 'right_displaced_atoms'],
    ]
    for molecule_md, atom_type in atom_types:
        molecule = getattr(form, molecule_md)
        selected_hydrogens = []
        for atom_md in getattr(molecule, atom_type):
            if atom_md.molecule == core.Backbone:
                if atom_md.position < 1 or atom_md.position > form.backbone.structure.NumAtoms():
                    errors.append('Invalid position {} for {}.{}'.format(atom_md.position, molecule_md, atom_type))
                    continue

                atom = form.backbone.structure.GetAtom(atom_md.position)
                if atom_md.element == 'H' and atom.GetAtomicNum() != 1:
                    atom = core.get_hydrogen_atom(atom, selected_hydrogens)
                    if atom is None:
                        continue

                if element_table.GetSymbol(atom.GetAtomicNum()) != atom_md.element:
                    errors.append('Invalid element at position {} for {}.{}'.format(atom_md.position, molecule_md, atom_type))

    atom_types = [
        'monomer_bond_atoms',
        'monomer_displaced_atoms',
        'left_bond_atoms',
        'right_bond_atoms',
        'left_displaced_atoms',
        'right_displaced_atoms',
    ]
    for monomer in form.alphabet.monomers.values():
        for atom_type in atom_types:
            selected_hydrogens = []
            for atom_md in getattr(monomer, atom_type):
                if atom_md.molecule == core.Monomer:
                    if atom_md.position < 1 or atom_md.position > monomer.structure.NumAtoms():
                        errors.append('Invalid position {} for Monomer:{} {}'.format(atom_md.position, monomer.id, atom_type))
                        continue

                    atom = monomer.structure.GetAtom(atom_md.position)
                    if atom_md.element == 'H' and atom.GetAtomicNum() != 1:
                        atom = core.get_hydrogen_atom(atom, selected_hydrogens)
                        if atom is None:
                            continue

                    if element_table.GetSymbol(atom.GetAtomicNum()) != atom_md.element:
                        errors.append('Invalid element at position {} for Monomer:{} {}'.format(
                            atom_md.position, monomer.id, atom_type))

    if errors:
        raise ValueError('BpForm {} is invalid:\n  {}'.format(form_type.__name__, '\n  '.join(errors)))
