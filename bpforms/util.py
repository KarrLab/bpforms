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
        protein.curated_protein_alphabet,
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
    if alphabet == 'curated_protein':
        return protein.CuratedProteinForm

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
        url = monomer.get_image_url()
        if url:
            doc += '          <td><img src="{}"/></td>'.format(url)
        else:
            doc += '          <td></td>'.format(url)
        doc += '        </tr>'
    doc += '    </table>'
    doc += '  </body>'
    doc += '</html>'

    with open(filename, 'w') as file:
        file.write(doc)
