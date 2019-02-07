""" Alphabet and BpForm to represent modified DNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Alphabet, AlphabetBuilder, Base, BpForm, Identifier, IdentifierSet, SynonymSet
from wc_utils.util.chem import EmpiricalFormula
import os.path
import pkg_resources

filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'dna.yml'))
dna_alphabet = Alphabet().from_yaml(filename)
# :obj:`Alphabet`: Alphabet for DNA nucleotides


class DnaAlphabetBuilder(AlphabetBuilder):
    """ Build DNA alphabet from DNAmod """

    def run(self, path=filename):
        """ Build alphabet and, optionally, save to YAML file

        Args:
            path (:obj:`str`, optional): path to save alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        return super(DnaAlphabetBuilder, self).run(path)

    def build(self):
        """ Build alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        # initialize alphabet
        alphabet = Alphabet()

        # load canonical bases
        alphabet.from_yaml(pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'dna.canonical.yml')))

        # create bases
        alphabet.A = Base(
            id='dAMP',
            name="2'-deoxyadenosine 5'-monophosphate(2−)",
            synonyms=SynonymSet([
                'dAMP(2−)',
                'dAMP dianion',
            ]),
            identifiers=IdentifierSet([
                Identifier('pubchem.compound', '22848660'),
                Identifier('metacyc.compound', 'DAMP'),
                Identifier('chebi', 'CHEBI:58245'),
            ]),
            structure=('InChI=1S/C10H14N5O6P/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(16)6(21-7)2-20-22(17,18)19'
                       '/h3-7,16H,1-2H2,(H2,11,12,13)(H2,17,18,19)/p-2/t5-,6+,7+/m0/s1')
        )

        # return alphabet
        return alphabet


class DnaForm(BpForm):
    """ DNA form """

    def __init__(self, base_seq=None):
        """
        Args:
            base_seq (:obj:`BaseSequence`, optional): bases of the DNA form
        """
        super(DnaForm, self).__init__(base_seq=base_seq, alphabet=dna_alphabet,
                                      bond_formula=EmpiricalFormula('H') * -1, bond_charge=1)
