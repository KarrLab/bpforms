""" Alphabet and BpForm to represent modified proteins

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from .core import Alphabet, Base, BpForm, Identifier
from wc_utils.util.chem import EmpiricalFormula
import pkg_resources

filename = pkg_resources.resource_filename('bpforms', 'protein.yml')
protein_alphabet = Alphabet.from_yaml(filename)
# :obj:`Alphabet`: Alphabet for protein amino acids


class ProteinForm(BpForm):
    """ Protein form """
    DEFAULT_ALPHABET = protein_alphabet
    DEFAULT_BOND_FORMULA = EmpiricalFormula('H2O') * -1
    DEFAULT_BOND_CHARGE = 0

    def __init__(self, base_seq=None):
        """
        Args:
            base_seq (:obj:`BaseSequence`, optional): bases of the DNA form
        """
        super(ProteinForm, self).__init__(base_seq=base_seq)
