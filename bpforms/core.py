""" bpforms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

import Bio
import openbabel


class BpForm(object):
    """ Biopolymer form 

    Attributes:
        structure (:obj:`str`): structure of the biopolymer
    """

    def __init__(self, structure, ph=7.):
        self.structure = structure
        self.ph = ph

    def get_formula(self):
        """ Get the chemical formula

        Returns:
            :obj:`str`: chemical formula
        """
        pass

    def protonate(self):
        """ Get the major protonation state of the modifications of the biopolymer at a pH

        Args:
            ph (obj:`float`): pH
        """
        pass

    def gen_visualization(self, filename):
        """ Generate a visual depiction of the biopolymer and save it to a file

        Args:
            filename (:obj:`str`): path to save the visual depiction of the biopolymer
        """
        pass


class DnaForm(BpForm):
    """ DNA form """
    pass


class RnaForm(BpForm):
    """ RNA form """
    pass


class ProteinForm(BpForm):
    """ Protein form """
    pass


def get_form(type):
    if type == 'dna':
        return DnaForm
    if type == 'rna':
        return RnaForm
    if type == 'protein':
        return ProteinForm

    raise ValueError('Type "{}" must be "dna", "rna", or "protein"'.format(type))
