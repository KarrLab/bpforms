""" Utilities for BpForms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from .alphabet import dna
from .alphabet import rna
from .alphabet import protein


def get_form(alphabet):
    """ Get a subclass of BpFrom

    Args:
        alphabet (:obj:`str`): alphabet

    Returns:
        :obj:`type`: subclass of BpForm
    """
    if alphabet == 'dna':
        return dna.DnaForm
    if alphabet == 'rna':
        return rna.RnaForm
    if alphabet == 'protein':
        return protein.ProteinForm

    raise ValueError('Alphabet "{}" must be "dna", "rna", or "protein"'.format(alphabet))


def build_alphabets():
    """ Build DNA, RNA, and protein alphabets """
    dna.DnaAlphabetBuilder().run()
    rna.RnaAlphabetBuilder().run()
    protein.ProteinAlphabetBuilder().run()
