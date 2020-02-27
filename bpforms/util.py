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
from wc_utils.util.chem import OpenBabelUtils
from wc_utils.util.chem.marvin import draw_molecule
import importlib
import jinja2
import math
import openbabel
import os
import pkg_resources
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


def gen_genomic_viz(polymers, inter_crosslinks=None, polymer_labels=None, seq_features=None,
                    width=800, cols=1, polymer_margin=25,
                    nt_per_track=100, track_sep=10,
                    polymer_label_font_size=15, seq_font_size=13, tick_label_font_size=10,
                    legend_font_size=13, tooltip_font_size=13,
                    x_link_stroke_width=2, x_link_radius=4,
                    nick_stroke_width=2,
                    axis_stroke_width=0.5,
                    seq_color='#000000', non_canonical_color='#e74624',
                    intra_x_link_color='#2daae1', inter_x_link_color='#90e227',
                    nick_color='#dabe2e',
                    axis_color='#000000', polymer_label_color='#000000'):
    """ Get a genomic visualization of the :obj:`BpForm`

    Args:
        polymers (:obj:`list` of :obj:`core.BpForm`): polymers
        inter_crosslinks (:obj:`list`): list of inter-polymer crosslinks
        polymer_labels (:obj:`dict`, optional): dictionary that maps polymers to their labels
        seq_features (:obj:`list` of :obj:`dict`, optional): list of features each
            represented by a dictionary with three keys

            * label (:obj:`str`): description of the type of feature
            * color (:obj:`str`): color
            * positions (:obj:`list` of :obj:`dict`): dictionary which maps 
              indices of polymers to a list of position ranges of the type 
              of feature
        width (:obj:`int`, optional): width
        cols (:obj:`int`, optional): number of columns of polymers
        polymer_margin (:obj:`int`, optional): horizontal and vertical spacing between polymers
        nt_per_track (:obj:`int`, optional): number of nucleotides per track
        track_sep (:obj:`int`, optional): vertical separation between tracks in pixels
        polymer_label_font_size (:obj:`float`, optional): font size of polymer label
        seq_font_size (:obj:`float`, optional): font size of sequence
        tick_label_font_size (:obj:`float`, optional): font size of tick labels
        legend_font_size (:obj:`float`, optional): font size of legend
        tooltip_font_size (:obj:`float`, optional): font size of tooltip
        x_link_stroke_width (:obj:`float`, optional): stroke width of crosslinks
        x_link_radius (:obj:`float`, optional): radius of crosslinks line
        nick_stroke_width (:obj:`float`, optional): stroke width of nicks
        axis_stroke_width (:obj:`float`, optional): stroke width of axis
        seq_color (:obj:`str`, optional): color of canonical monomers
        non_canonical_color (:obj:`str`, optional): color of non-canonical monomers
        intra_x_link_color (:obj:`str`, optional): colors of intrastrand crosslinks
        inter_x_link_color (:obj:`str`, optional): colors of interstrand crosslinks
        nick_color (:obj:`str`, optional): colors of nicks
        axis_color (:obj:`str`, optional): color of axis
        polymer_label_color (:obj:`str`, optional): color of polymer labels

    Returns:
        :obj:`str`: SVG image
    """
    import bpforms
    from bpforms.xlink.core import onto_crosslink_to_id

    inter_crosslinks = inter_crosslinks or []
    polymer_labels = polymer_labels or {}
    seq_features = seq_features or []

    nc_label_sep = 0.1 * seq_font_size
    axis_sep = 0.2 * tick_label_font_size
    tick_len = 0.4 * tick_label_font_size
    tick_label_sep = 0.6 * tick_label_font_size
    legend_label_sep = 1/3 * legend_font_size
    if not polymer_labels:
        polymer_label_font_size = 0
    polymer_label_sep = 0.5 * polymer_label_font_size

    has_nicks = False

    polymers_context = []
    for i_polymer, polymer in enumerate(polymers):
        if isinstance(polymer, dna.DnaForm):
            canonical_monomers = [dna.dna_alphabet.monomers.get(code) for code in 'ACGT']
        elif isinstance(polymer, rna.RnaForm):
            canonical_monomers = [rna.rna_alphabet.monomers.get(code) for code in 'ACGU']
        elif isinstance(polymer, protein.ProteinForm):
            canonical_monomers = [protein.protein_alphabet.monomers.get(code)
                                  for code in bpforms.canonical_protein_alphabet.monomers.keys()]
        else:
            raise ValueError('BpForm must be an instance of `DnaForm`, `RnaForm`, or `ProteinForm`')

        monomer_codes = {monomer: code for code, monomer in polymer.alphabet.monomers.items()}

        seq_tracks = []
        monomer_seq = polymer.seq
        canonical_seq = polymer.get_canonical_seq()

        n_tracks = math.ceil(len(monomer_seq) / nt_per_track)
        max_nc_label_len = 0
        has_non_canonical = False
        for i_track in range(n_tracks):
            seq_track = {
                'i_track': i_track,
                'min': i_track * nt_per_track + 1,
                'max': min((i_track + 1) * nt_per_track, len(monomer_seq)),
                'len': min((i_track + 1) * nt_per_track, len(monomer_seq)) - (i_track * nt_per_track + 1) + 1,
                'seq': [],
                'ticks': [],
            }
            seq_tracks.append(seq_track)

            # seq
            for i_monomer, (monomer, canonical_code) in enumerate(zip(monomer_seq[i_track * nt_per_track:(i_track + 1) * nt_per_track],
                                                                      canonical_seq[i_track * nt_per_track:(i_track + 1) * nt_per_track])):
                non_canonical_label = monomer_codes.get(monomer, None)
                if not non_canonical_label:
                    if monomer.id:
                        non_canonical_label = monomer.id
                    elif monomer.name:
                        non_canonical_label = monomer.name
                    elif monomer.synonyms:
                        non_canonical_label = list(monomer.synonyms)[0]
                    elif monomer.identifiers:
                        non_canonical_label = list(monomer.identifiers)[0].id
                    else:
                        non_canonical_label = 'Non-canonical'

                if monomer.name:
                    tooltip = monomer.name
                elif monomer.synonyms:
                    tooltip = list(monomer.synonyms)[0]
                elif monomer.identifiers:
                    tooltip = list(monomer.identifiers)[0].id
                elif monomer.id:
                    tooltip = monomer.id
                else:
                    tooltip = 'Non-canonical monomer'

                track_monomer = {
                    'i_monomer': i_monomer,
                    'canonical_code': canonical_code,
                    'canonical': monomer in canonical_monomers,
                    'non_canonical_label': non_canonical_label,
                    'tooltip': tooltip,
                }
                has_non_canonical = has_non_canonical or track_monomer['canonical']

                pos = i_track * nt_per_track + i_monomer + 1
                for feature in seq_features:
                    for position in feature['positions'].get(i_polymer, []):
                        if pos >= position[0] and pos <= position[1]:
                            track_monomer['color'] = feature['color']
                            break

                seq_track['seq'].append(track_monomer)

                max_nc_label_len = max(max_nc_label_len, len(non_canonical_label))

            # ticks
            seq_track['ticks'].append({'x': 1, 'label': i_track * nt_per_track + 1})
            for i_tick in range(math.floor((seq_track['max'] - seq_track['min'] + 1) / 10)):
                seq_track['ticks'].append({'x': (i_tick + 1) * 10, 'label': i_track * nt_per_track + (i_tick + 1) * 10})

        x_links = []
        for crosslink in polymer.crosslinks:
            l = crosslink.get_l_bond_atoms()[0].monomer
            r = crosslink.get_r_bond_atoms()[0].monomer
            l_track = math.floor((l - 1) / nt_per_track)
            r_track = math.floor((r - 1) / nt_per_track)
            l_pos = (l - 1) % nt_per_track + 1
            r_pos = (r - 1) % nt_per_track + 1

            if l_pos < r_pos:
                x_link = {
                    'l_track': l_track,
                    'r_track': r_track,
                    'l_pos': l_pos,
                    'r_pos': r_pos,
                }
            else:
                x_link = {
                    'l_track': r_track,
                    'r_track': l_track,
                    'l_pos': r_pos,
                    'r_pos': l_pos,
                }
            if isinstance(crosslink, core.OntoBond):
                x_link['tooltip'] = onto_crosslink_to_id[crosslink.type]
            else:
                x_link['tooltip'] = None
            x_links.append(x_link)

        if polymer.nicks:
            has_nicks = True
        nicks = []
        for nick in polymer.nicks:
            nicks.append({
                'pos': (nick.position - 1) % nt_per_track + 1,
                'track': math.floor((nick.position - 1) / nt_per_track),
                'tooltip': None,
            })

        polymers_context.append({
            'label': polymer_labels.get(i_polymer, None),
            'seq_tracks': seq_tracks,
            'x_links': x_links,
            'nicks': nicks,
        })

    inter_x_links_context = []
    for x_link in inter_crosslinks:
        l = int(float(x_link.get_l_bond_atoms()[0].subunit))
        r = int(float(x_link.get_r_bond_atoms()[0].subunit))
        l_polymer_row = math.floor(l / cols)
        l_polymer_col = l % cols
        r_polymer_row = math.floor(r / cols)
        r_polymer_col = r % cols

        l = x_link.get_l_bond_atoms()[0].monomer
        r = x_link.get_r_bond_atoms()[0].monomer
        l_track = math.floor((l - 1) / nt_per_track)
        r_track = math.floor((r - 1) / nt_per_track)
        l_pos = (l - 1) % nt_per_track + 1
        r_pos = (r - 1) % nt_per_track + 1

        if l_polymer_col > r_polymer_col \
            or (l_polymer_col == r_polymer_col and
                l_pos > r_pos):
            inter_x_link_context = {
                'r_polymer_row': l_polymer_row,
                'r_polymer_col': l_polymer_col,
                'r_track': l_track,
                'r_pos': l_pos,
                'l_polymer_row': r_polymer_row,
                'l_polymer_col': r_polymer_col,
                'l_track': r_track,
                'l_pos': r_pos,
            }
        else:
            inter_x_link_context = {
                'l_polymer_row': l_polymer_row,
                'l_polymer_col': l_polymer_col,
                'l_track': l_track,
                'l_pos': l_pos,
                'r_polymer_row': r_polymer_row,
                'r_polymer_col': r_polymer_col,
                'r_track': r_track,
                'r_pos': r_pos,
            }
        inter_x_link_context['tooltip'] = getattr(x_link, 'comments', None)
        inter_x_links_context.append(inter_x_link_context)

    all_x_links = list(inter_x_links_context)
    for i_polymer, polymer in enumerate(polymers_context):
        for x_link in polymer['x_links']:
            x_link['l_polymer_row'] = x_link['r_polymer_row'] = math.floor(i_polymer / cols)
            x_link['l_polymer_col'] = x_link['r_polymer_col'] = i_polymer % cols
            all_x_links.append(x_link)

    sorted_x_links = sorted(all_x_links, key=lambda x: (
        x['l_polymer_row'] == x['r_polymer_row'] and x['l_track'] == x['r_track'],
        x['l_polymer_row'], x['r_polymer_row'],
        x['l_track'], x['r_track'],
        x['l_polymer_col'], x['l_pos'],
        x['r_polymer_col'], x['r_pos']))
    offset = 0
    prev_horz = False
    prev_row = -1
    prev_col = -1
    prev_track = -1
    prev_pos = -1
    for x_link in sorted_x_links:
        if prev_horz and x_link['l_polymer_row'] == x_link['r_polymer_row'] and \
                x_link['l_track'] == x_link['r_track'] and \
                x_link['l_polymer_row'] == prev_row and \
                x_link['l_track'] == prev_track and \
                (x_link['l_polymer_col'] < prev_col or
                 (x_link['l_polymer_col'] == prev_col and x_link['l_pos'] <= prev_pos)):
            offset += 1
        else:
            offset = 0
            prev_horz = x_link['l_polymer_row'] == x_link['r_polymer_row'] and \
                x_link['l_track'] == x_link['r_track']
            prev_row = x_link['r_polymer_row']
            prev_col = x_link['r_polymer_col']
            prev_track = x_link['r_track']
            prev_pos = x_link['r_pos']
        x_link['v_offset'] = offset

    sorted_x_links = sorted(all_x_links, key=lambda x: (
        x['l_polymer_col'] == x['r_polymer_col'] and x['l_pos'] == x['r_pos'],
        x['l_polymer_col'], x['r_polymer_col'],
        x['l_pos'], x['r_pos'],
        x['l_polymer_row'], x['l_track'],
        x['r_polymer_row'], x['r_track'],))
    offset = 0
    prev_vert = False
    prev_row = -1
    prev_col = -1
    prev_track = -1
    prev_pos = -1
    for x_link in sorted_x_links:
        if prev_vert and x_link['l_polymer_col'] == x_link['r_polymer_col'] and \
                x_link['l_pos'] == x_link['r_pos'] and \
                x_link['l_polymer_col'] == prev_col and \
                x_link['l_pos'] == prev_pos and \
                (x_link['l_polymer_row'] < prev_row or
                 (x_link['l_polymer_row'] == prev_row and x_link['l_track'] <= prev_track)):
            offset += 1
        else:
            offset = 0
            prev_vert = x_link['l_polymer_col'] == x_link['r_polymer_col'] and \
                x_link['l_pos'] == x_link['r_pos']
            prev_row = x_link['r_polymer_row']
            prev_col = x_link['r_polymer_col']
            prev_track = x_link['r_track']
            prev_pos = x_link['r_pos']
        x_link['h_offset'] = offset

    # read template
    with open(pkg_resources.resource_filename('bpforms', 'genomic_viz.template.svg')) as file:
        template = jinja2.Template(file.read())

    # render template
    max_polymer_len = max(len(polymer.seq) for polymer in polymers)
    h_padding = tick_label_font_size * (math.floor(math.log10(max_polymer_len)) + 1) * 0.5 * 0.6
    code_h = seq_font_size
    nc_label_h = max_nc_label_len * math.sin(math.pi / 3) * seq_font_size * 0.65
    track_h = nc_label_h + nc_label_sep + code_h \
        + axis_sep + tick_len + tick_label_sep + tick_label_font_size * 0.25
    legend_sep = polymer_margin

    polymer_w = (width - (cols - 1) * polymer_margin) / cols
    max_n_tracks = math.ceil(max_polymer_len / nt_per_track)
    polymer_hs = [0] * math.ceil(len(polymers) / cols)
    for i_polymer, polymer in enumerate(polymers):
        n_tracks = math.ceil(len(polymer.seq) / nt_per_track)
        polymer_h = (polymer_label_font_size + polymer_label_sep) \
            + track_h * n_tracks + track_sep * (n_tracks - 1)
        i_row = math.floor(i_polymer / cols)
        polymer_hs[i_row] = max(polymer_h, polymer_hs[i_row])

    cum_polymer_h = []
    for i_row, polymer_h in enumerate(polymer_hs):
        cum_polymer_h.append(sum(polymer_hs[0:i_row]))

    legend_rows = []

    if has_non_canonical:
        legend_rows.append({
            'label': 'Non-canonical monomeric form',
            'color': non_canonical_color,
            'symbol': 'X',
        })
    for seq_feature in seq_features:
        legend_rows.append({
            'label': seq_feature['label'],
            'color': seq_feature['color'],
            'symbol': 'X',
        })
    if any(len(polymer.crosslinks) >= 1 for polymer in polymers):
        legend_rows.append({
            'label': 'Intrachain crosslink',
            'color': intra_x_link_color,
            'symbol': None,
            'stroke_width': x_link_stroke_width,
        })
    if inter_crosslinks:
        legend_rows.append({
            'label': 'Interchain crosslink',
            'color': inter_x_link_color,
            'symbol': None,
            'stroke_width': x_link_stroke_width,
        })
    if has_nicks:
        legend_rows.append({
            'label': 'Nick',
            'color': nick_color,
            'symbol': None,
            'stroke_width': nick_stroke_width,
        })

    context = {
        # BpForm
        'polymers': polymers_context,
        'inter_x_links': inter_x_links_context,

        # legend
        'legend_rows': legend_rows,

        # size
        'width': width,
        'height': sum(polymer_hs) + \
        polymer_margin * (math.ceil(len(polymers) / cols)-1) + \
        (len(legend_rows) >= 1) * (\
            legend_sep + \
            len(legend_rows) * legend_font_size + \
            (len(legend_rows) - 1) * legend_label_sep),
        'h_padding': h_padding,
        'cols': cols,
        'polymer_w': polymer_w,
        'cum_polymer_h': cum_polymer_h,
        'polymer_margin': polymer_margin,

        # track
        'nt_per_track': nt_per_track,
        'track_w': polymer_w - 2 * h_padding,
        'track_h': track_h,
        'px_per_nt': (polymer_w - 2 * h_padding) / (nt_per_track - 1),
        'nc_label_h': nc_label_h,
        'nc_label_sep': nc_label_sep,
        'code_h': code_h,
        'axis_sep': axis_sep,
        'tick_len': tick_len,
        'tick_label_sep': tick_label_sep,
        'track_sep': track_sep,
        'polymer_label_sep': polymer_label_sep,

        # crosslinks
        'x_link_radius': x_link_radius,

        # legend
        'legend_sep': legend_sep,
        'legend_label_sep': legend_label_sep,

        # font sizes
        'polymer_label_font_size': polymer_label_font_size,
        'seq_font_size': seq_font_size,
        'tick_label_font_size': tick_label_font_size,
        'legend_font_size': legend_font_size,

        # stroke widths
        'x_link_stroke_width': x_link_stroke_width,
        'nick_stroke_width': nick_stroke_width,
        'axis_stroke_width': axis_stroke_width,

        # colors
        'seq_color': seq_color,
        'non_canonical_color': non_canonical_color,
        'intra_x_link_color': intra_x_link_color,
        'inter_x_link_color': inter_x_link_color,
        'nick_color': nick_color,
        'axis_color': axis_color,
        'polymer_label_color': polymer_label_color,
    }

    return template.render(**context)


def export_ontos_to_obo(alphabets=None, filename=None, _max_monomers=None, _max_xlinks=None):
    """ Exports alphabets of residues and ontology of crosslinks to OBO format

    Args:
        alphabets (:obj:`list` of :obj:`core.Alphabet`, optional): alphabets to export        
        filename (:obj:`str`, optional): path to export alphabets
        _max_monomers (:obj:`int`, optional): maximum number of monomers to export
        _max_xlinks (:obj:`int`, optional): maximum number of crosslinks to export
    """
    from .xlink.core import crosslinks_onto
    import pronto

    if alphabets is None:
        alphabets = [dna.dna_alphabet, rna.rna_alphabet, protein.protein_alphabet]
    if filename is None:
        filename = pkg_resources.resource_filename('bpforms',
                                                   os.path.join('alphabet', 'onto.obo'))

    # create ontology
    onto = pronto.Ontology()

    onto.metadata.synonymtypedefs.add(pronto.SynonymType('structure', 'SMILES-encoded structure', 'EXACT'))
    onto.metadata.synonymtypedefs.add(pronto.SynonymType('backbone_bond_atoms', 'Backbone bond atoms', 'EXACT'))
    onto.metadata.synonymtypedefs.add(pronto.SynonymType('backbone_displaced_atoms', 'Backbone displaced atoms', 'EXACT'))
    onto.metadata.synonymtypedefs.add(pronto.SynonymType('l_bond_atoms', 'Left bond atoms', 'EXACT'))
    onto.metadata.synonymtypedefs.add(pronto.SynonymType('l_displaced_atoms', 'Left displaced atoms', 'EXACT'))
    onto.metadata.synonymtypedefs.add(pronto.SynonymType('r_bond_atoms', 'Right bond atoms', 'EXACT'))
    onto.metadata.synonymtypedefs.add(pronto.SynonymType('r_displaced_atoms', 'Right displaced atoms', 'EXACT'))
    onto.metadata.synonymtypedefs.add(pronto.SynonymType('delta_mass', 'Delta mass', 'EXACT'))
    onto.metadata.synonymtypedefs.add(pronto.SynonymType('delta_charge', 'Delta charge', 'EXACT'))
    onto.metadata.synonymtypedefs.add(pronto.SynonymType('bond_order', 'Bond order', 'EXACT'))
    onto.metadata.synonymtypedefs.add(pronto.SynonymType('bond_stereo', 'Bond stereochemistry', 'EXACT'))

    is_a_relationship = onto.get_relationship('is_a')
    part_of_relationship = onto.create_relationship('part_of')
    derives_from_relationship = onto.create_relationship('derives_from')
    l_monomer_relationship = onto.create_relationship('l_monomer')
    r_monomer_relationship = onto.create_relationship('r_monomer')

    # add terms for monomers to ontology
    alphabet_type_term = onto.create_term(id='BpForms:alphabet')
    alphabet_type_term.name = 'alphabet'

    monomer_type_term = onto.create_term(id='BpForms:monomer')
    monomer_type_term.name = 'monomeric form'
    monomer_type_term.add_synonym('residue', scope='BROAD')
    monomer_type_term.add_synonym('monomer', scope='BROAD')

    xlink_type_term = onto.create_term(id='BpForms:crosslink')
    xlink_type_term.name = 'crosslink'

    monomer_to_term = {}
    for alphabet in alphabets:
        alphabet_term = onto.create_term(id='BpForms:' + alphabet.id)
        alphabet_term.name = alphabet.name
        alphabet_term.desc = alphabet.description
        alphabet_term.relationships = {
            is_a_relationship: [alphabet_type_term],
        }

        for code, monomer in list(alphabet.monomers.items())[0:_max_monomers]:
            synonyms = []
            for syn in monomer.synonyms:
                synonyms.append(pronto.SynonymData(syn, scope='BROAD'))
            for identifier in monomer.identifiers:
                synonyms.append(pronto.SynonymData(
                    '{}: {}'.format(identifier.ns, identifier.id),
                    scope='EXACT'))

            relationships = {
                is_a_relationship: [monomer_type_term],
                part_of_relationship: [alphabet_term],
            }

            other = {'structure': [monomer.export('smiles')]}
            for atom_type in ['backbone_bond_atoms', 'backbone_displaced_atoms',
                              'l_bond_atoms', 'l_displaced_atoms',
                              'r_bond_atoms', 'r_displaced_atoms']:
                atoms = getattr(monomer, atom_type)
                if atoms:
                    other[atom_type] = [_atom_to_str(atom) for atom in atoms]
            if monomer.delta_mass:
                other['delta_mass'] = [monomer.delta_mass]
            if monomer.delta_charge:
                other['delta_charge'] = [monomer.delta_charge]

            term = onto.create_term(id='BpForms:{}_{}'.format(alphabet.id, code))
            term.name = monomer.name or None
            term.desc = monomer.comments or ''
            term.relationships = relationships
            term.other = other
            for syn in synonyms:
                term.add_synonym(syn.description, scope=syn.scope)

            monomer_to_term[monomer] = term

    for alphabet in alphabets:
        monomer_codes = {monomer: code for code, monomer in alphabet.monomers.items()}
        for code, monomer in list(alphabet.monomers.items())[0:_max_monomers]:
            monomer_term = onto['BpForms:{}_{}'.format(alphabet.id, code)]
            relationships = {}
            for base_monomer in monomer.base_monomers:
                base_monomer_term = onto.get('BpForms:{}_{}'.format(
                    alphabet.id, monomer_codes[base_monomer]), None)
                if base_monomer_term is not None:
                    relationships[derives_from_relationship] = [base_monomer_term]
            monomer_term.relationships = relationships

    # add terms for crosslinks to ontology
    for xlink in list(crosslinks_onto.values())[0:_max_xlinks]:
        relationships = {
            is_a_relationship: [xlink_type_term],
        }

        l_monomer_term = monomer_to_term.get(xlink.l_monomer, None)
        r_monomer_term = monomer_to_term.get(xlink.r_monomer, None)
        if l_monomer_term:
            relationships[l_monomer_relationship] = [l_monomer_term]
        if r_monomer_term:
            relationships[r_monomer_relationship] = [r_monomer_term]

        other = {}
        for atom_type in ['l_bond_atoms', 'r_bond_atoms', 'l_displaced_atoms', 'r_displaced_atoms']:
            other[atom_type] = [_atom_to_str(atom) for atom in getattr(xlink, atom_type)]

        other['bond_order'] = [xlink.order.name]
        if xlink.stereo:
            other['bond_stereo'] = [xlink.stereo.name]

        term = onto.create_term(id='BpForms:crosslink_{}'.format(xlink.id))
        term.name = xlink.name or None
        for syn in xlink.synonyms:
            term.add_synonym(syn, scope='BROAD')
        term.desc = xlink.comments or ''
        term.relationships = relationships
        term.other = other

    # export ontology
    with open(filename, 'wb') as file:
        onto.dump(file)


def _atom_to_str(atom):
    """ Get the string representation of an atom in a monomeric form in an
    alphabet

    Args:
        atom (:obj:`core.Atom`): atom

    Returns:
        :obj:`str`: string representation
    """
    if atom.charge:
        return '{0}{1}{2:+d}'.format(atom.element, atom.position, atom.charge)
    else:
        return '{0}{1}'.format(atom.element, atom.position)
