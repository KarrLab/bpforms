""" Alphabet and BpForm to represent modified proteins

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Alphabet, AlphabetBuilder, Base, BpForm, Identifier, IdentifierSet, SynonymSet
from wc_utils.util.chem import EmpiricalFormula
import os.path
import pkg_resources

filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.yml'))
protein_alphabet = Alphabet().from_yaml(filename)
# :obj:`Alphabet`: Alphabet for protein amino acids


class ProteinAlphabetBuilder(AlphabetBuilder):
    """ Build protein alphabet from RESID """

    def run(self, path=filename):
        """ Build alphabet and, optionally, save to YAML file

        Args:
            path (:obj:`str`, optional): path to save alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        return super(ProteinAlphabetBuilder, self).run(path)

    def build(self):
        """ Build alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        # initialize alphabet
        alphabet = Alphabet()

        # create bases
        alphabet.A = Base(
            id='ALA',
            name="alanine",
            synonyms=SynonymSet([
                'L-alpha-alanine',
            ]),
            identifiers=IdentifierSet([
                Identifier('pubchem.compound', '7311724'),
                Identifier('metacyc.compound', 'L-ALPHA-ALANINE'),
                Identifier('chebi', 'CHEBI:57972'),
            ]),
            structure='InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1',
        )

        # return alphabet
        return alphabet


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
