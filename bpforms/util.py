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
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from wc_utils.util.chem import draw_molecule, OpenBabelUtils
import importlib
import openbabel
import shutil


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


def build_alphabets(ph=None, major_tautomer=False, dearomatize=False, _max_monomers=float('inf'), alphabets=None):
    """ Build DNA, RNA, and protein alphabets

    Args:
        ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomeric form
        major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
        dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule
        _max_monomers (:obj:`float`, optional): maximum number of monomeric forms to build; used for testing
        alphabets (:obj:`list` of :obj:`str` or :obj:`None`, optional): ids of alphabets to build. If :obj:`None`,
            build all alphabets
    """
    core.cache.clear()

    if alphabets is None or 'dna' in alphabets:
        dna.DnaAlphabetBuilder(_max_monomers=_max_monomers).run(ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)
    if alphabets is None or 'rna' in alphabets:
        rna.RnaAlphabetBuilder(_max_monomers=_max_monomers).run(ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)
    if alphabets is None or 'protein' in alphabets:
        protein.ProteinAlphabetBuilder(_max_monomers=_max_monomers).run(ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)


def gen_html_viz_alphabet(bpform_type, filename):
    """ Create and save an HTML document with images of the monomeric forms in an alphabet

    Args:
        bpform_type (:obj:`type`): subclass of :obj:`core.BpForm`
        filename (:obj:`str`): path to save HTML document with images of monomeric forms
    """
    width = 400
    height = 400

    bpform = bpform_type()
    alphabet = bpform.alphabet

    doc = ''
    doc += '<html>\n'
    doc += '  <style type="text/css">\n'
    doc += '  table {width: 100%;}\n'
    doc += '  thead th {background: #ccc}\n'
    doc += '  tbody tr:nth-child(even) {background: #dedede}\n'
    doc += '  tr td, tr th {padding:5px; text-align:center;}\n'
    doc += '  tr td:first-child, tr tr:first-child {padding-left:10px;}\n'
    doc += '  tr td:last-child, tr tr:last-child {padding-right:10px;}\n'
    doc += '  </style>\n'
    doc += '  <body>\n'
    doc += '    <table cellpadding="0" cellspacing="0">\n'
    doc += '      <thead>\n'
    doc += '        <tr>\n'
    doc += '          <th>Code</th>\n'
    doc += '          <th>Monomeric form</th>\n'
    doc += '          <th>Dimer</th>\n'
    doc += '          <th>Canonical SMILES</th>\n'
    doc += '          <th>Monomeric form bond atoms</th>\n'
    doc += '          <th>Monomeric form displaced atoms</th>\n'
    doc += '          <th>Left bond atoms</th>\n'
    doc += '          <th>Left displaced atoms</th>\n'
    doc += '          <th>Right bond atoms</th>\n'
    doc += '          <th>Right displaced atoms</th>\n'
    doc += '        </tr>\n'
    doc += '      </thead>\n'
    for code, monomer in alphabet.monomers.items():
        doc += '        <tr>\n'
        doc += '          <td>{}</td>\n'.format(code)
        doc += '          <td>{}</td>\n'.format(monomer.get_image(width=width, height=height, include_xml_header=False))

        if monomer.structure and bpform.can_monomer_bond_left(monomer) and bpform.can_monomer_bond_right(monomer):
            dimer = bpform_type()
            dimer.seq.append(monomer)
            dimer.seq.append(monomer)
            doc += '          <td>{}</td>\n'.format(draw_molecule(dimer.export('cml'), 'cml',
                                                                  show_atom_nums=False, width=width, height=height))
        else:
            doc += '          <td>{}</td>\n'.format('')
        doc += '          <td>{}</td>\n'.format(monomer.export('smiles', options=('c',)))
        doc += '          <td>{}</td>\n'.format(', '.join('{}{}'.format(atom.element, atom.position)
                                                          for atom in monomer.backbone_bond_atoms))
        doc += '          <td>{}</td>\n'.format(', '.join('{}{}'.format(atom.element, atom.position)
                                                          for atom in monomer.backbone_displaced_atoms))
        doc += '          <td>{}</td>\n'.format(', '.join('{}{}'.format(atom.element, atom.position) for atom in monomer.r_bond_atoms))
        doc += '          <td>{}</td>\n'.format(', '.join('{}{}'.format(atom.element, atom.position)
                                                          for atom in monomer.r_displaced_atoms))
        doc += '          <td>{}</td>\n'.format(', '.join('{}{}'.format(atom.element, atom.position) for atom in monomer.l_bond_atoms))
        doc += '          <td>{}</td>\n'.format(', '.join('{}{}'.format(atom.element, atom.position)
                                                          for atom in monomer.l_displaced_atoms))
        doc += '        </tr>\n'
    doc += '    </table>\n'
    doc += '  </body>\n'
    doc += '</html>\n'

    with open(filename, 'w') as file:
        file.write(doc)


def validate_bpform_bonds(form_type):
    """ Validate bonds in alphabet

    Args:
        form_type (:obj:`type`): type of BpForm

    Raises:
        :obj:`ValueError`: if any of the bonds are invalid
    """

    form = form_type()

    element_table = openbabel.OBElementTable()

    errors = []

    # validate bonds to backbone
    atom_types = [
        ['backbone', 'monomer_bond_atoms'],
        ['backbone', 'monomer_displaced_atoms'],
        ['bond', 'l_bond_atoms'],
        ['bond', 'r_bond_atoms'],
        ['bond', 'l_displaced_atoms'],
        ['bond', 'r_displaced_atoms'],
    ]
    for molecule_md, atom_type in atom_types:
        molecule = getattr(form, molecule_md)
        selected_hydrogens = []
        for atom_md in getattr(molecule, atom_type):
            if atom_md.molecule == core.Backbone:
                if form.backbone.structure:
                    n_backbone_atoms = form.backbone.structure.NumAtoms()
                else:
                    n_backbone_atoms = 0
                if atom_md.position < 1 or atom_md.position > n_backbone_atoms:
                    errors.append('Invalid position {} for {}.{}'.format(atom_md.position, molecule_md, atom_type))
                    continue

                atom = form.backbone.structure.GetAtom(atom_md.position)
                if atom_md.element == 'H' and atom.GetAtomicNum() != 1:
                    atom = core.get_hydrogen_atom(atom, selected_hydrogens, None)
                    if atom is None:
                        continue

                if element_table.GetSymbol(atom.GetAtomicNum()) != atom_md.element:
                    errors.append('Invalid element {} != {} at position {} for {}.{}'.format(
                        element_table.GetSymbol(atom.GetAtomicNum()), atom_md.element,
                        atom_md.position, molecule_md, atom_type))

    # validate bonds to monomer
    atom_types = [
        'backbone_bond_atoms',
        'backbone_displaced_atoms',
        'r_bond_atoms',
        'l_bond_atoms',
        'r_displaced_atoms',
        'l_displaced_atoms',
    ]
    for i_monomer, monomer in enumerate(form.alphabet.monomers.values()):
        for atom_type in atom_types:
            selected_hydrogens = []
            for atom_md in getattr(monomer, atom_type):
                if atom_md.molecule == core.Monomer:
                    if atom_md.position < 1 or atom_md.position > monomer.structure.NumAtoms():
                        errors.append('Invalid position {} for monomeric form:{} {}'.format(atom_md.position, monomer.id, atom_type))
                        continue

                    atom = monomer.structure.GetAtom(atom_md.position)
                    if atom_md.element == 'H' and atom.GetAtomicNum() != 1:
                        atom = core.get_hydrogen_atom(atom, selected_hydrogens, i_monomer)
                        if atom is None:
                            continue

                    if element_table.GetSymbol(atom.GetAtomicNum()) != atom_md.element:
                        errors.append('Invalid element {} != {} at position {} for monomeric form:{} {}'.format(
                            element_table.GetSymbol(atom.GetAtomicNum()), atom_md.element,
                            atom_md.position, monomer.id, atom_type))

    # validate monomeric forms and dimers
    for monomer in form.alphabet.monomers.values():
        monomer_form = form_type(seq=[monomer])
        try:
            monomer_structure = monomer_form.get_structure()[0]
            if monomer_form.get_formula() != OpenBabelUtils.get_formula(monomer_structure):
                errors.append('Monomeric form of {} has incorrect formula: {} != {}'.format(
                    monomer.id, str(monomer_form.get_formula()), str(OpenBabelUtils.get_formula(monomer_structure))))
                continue
            if monomer_form.get_charge() != monomer_structure.GetTotalCharge():
                errors.append('Monomeric form of {} has incorrect charge: {} != {}'.format(
                    monomer.id, monomer_form.get_charge(), monomer_structure.GetTotalCharge()))
                continue
            OpenBabelUtils.export(monomer_structure, 'smiles')
            OpenBabelUtils.export(monomer_structure, 'inchi')
        except Exception as error:
            errors.append('Unable to create monomeric form of {}:\n    {}'.format(monomer.id, str(error)))

        if form.can_monomer_bond_left(monomer) and form.can_monomer_bond_right(monomer):
            dimer_form = form_type(seq=[monomer, monomer])
            try:
                dimer_structure = dimer_form.get_structure()[0]
                if dimer_form.get_formula() != OpenBabelUtils.get_formula(dimer_structure):
                    errors.append('Dimer of {} has incorrect formula: {} != {}'.format(
                        monomer.id, str(dimer_form.get_formula()), str(OpenBabelUtils.get_formula(dimer_structure))))
                    continue
                if dimer_form.get_charge() != dimer_structure.GetTotalCharge():
                    errors.append('Dimer of {} has incorrect charge: {} != {}'.format(
                        monomer.id, dimer_form.get_charge(), dimer_structure.GetTotalCharge()))
                    continue
                OpenBabelUtils.export(dimer_structure, 'smiles')
                OpenBabelUtils.export(dimer_structure, 'inchi')
            except Exception as error:
                errors.append('Unable to form dimer of {}:\n    {}'.format(monomer.id, str(error)))

    # report errors
    if errors:
        raise ValueError('BpForm {} is invalid:\n  {}'.format(form_type.__name__, '\n  '.join(errors)))


def write_to_fasta(forms, filename):
    """ Write BpForms to a FASTA-formatted file

    Args:
        forms (:obj:`dict`): dictionary which maps the ids of molecules to their BpForms-encoded
            sequences
        filename (:obj:`str`): path to FASTA-formatted file
    """
    seqs = (SeqRecord(id=id, seq=Seq(str(form))) for id, form in forms.items())
    SeqIO.write(seqs, filename, "fasta")


def read_from_fasta(filename, alphabet):
    """ Read BpForms from a FASTA-formatted file

    Args:
        filename (:obj:`str`): path to FASTA-formatted file
        alphabet (:obj:`str`): alphabet of BpForms in file

    Returns:
        :obj:`dict`: dictionary which maps the ids of molecules to their BpForms-encoded
            sequences
    """
    Form = get_form(alphabet)

    forms = {}
    for record in SeqIO.parse(filename, "fasta"):
        forms[record.id] = Form().from_str(str(record.seq))

    return forms
