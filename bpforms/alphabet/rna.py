""" Alphabet and BpForm to represent modified RNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Alphabet, AlphabetBuilder, Base, BpForm, Identifier, IdentifierSet, SynonymSet
from wc_utils.util.chem import EmpiricalFormula
import os.path
import pkg_resources

filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'rna.yml'))
rna_alphabet = Alphabet().from_yaml(filename)
# :obj:`Alphabet`: Alphabet for RNA nucleotides


class RnaAlphabetBuilder(AlphabetBuilder):
    """ Build RNA alphabet from MODOMICS """

    def run(self, path=filename):
        """ Build alphabet and, optionally, save to YAML file

        Args:
            path (:obj:`str`, optional): path to save alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        return super(RnaAlphabetBuilder, self).run(path)

    def build(self):
        """ Build alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        # initialize alphabet
        alphabet = Alphabet()

        # create bases
        alphabet.A = Base(
            id='AMP',
            name="2'-adenosine 5'-monophosphate(2−)",
            identifiers=IdentifierSet([
                Identifier('pubchem.compound', '15938965'),
                Identifier('metacyc.compound', 'AMP'),
                Identifier('chebi', 'CHEBI:456215'),
            ]),
            structure=('InChI=1S/C10H14N5O7P/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(22-10)1-21-23(18,19)20'
                       '/h2-4,6-7,10,16-17H,1H2,(H2,11,12,13)(H2,18,19,20)/p-2/t4-,6-,7-,10-/m1/s1')
        )

        # return alphabet
        return alphabet


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
