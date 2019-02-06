""" Alphabet and BpForm to represent modified RNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from .core import Alphabet, Base, BpForm, Identifier
from wc_utils.util.chem import EmpiricalFormula
import pkg_resources

filename = pkg_resources.resource_filename('bpforms', 'rna.yml')
rna_alphabet = Alphabet.from_yaml(filename)
# :obj:`Alphabet`: Alphabet for RNA nucleotides


class RnaForm(BpForm):
    """ RNA form """
    DEFAULT_ALPHABET = rna_alphabet
    DEFAULT_BOND_FORMULA = EmpiricalFormula('H') * -1
    DEFAULT_BOND_CHARGE = 1

    def __init__(self, base_seq=None):
        """
        Args:
            base_seq (:obj:`BaseSequence`, optional): bases of the DNA form
        """
        super(RnaForm, self).__init__(base_seq=base_seq)
