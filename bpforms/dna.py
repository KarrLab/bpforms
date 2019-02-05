""" Alphabet and BpForm to represent modified DNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from .core import Alphabet, Base, BpForm, Identifier
from wc_utils.util.chem import EmpiricalFormula

dna_alphabet = Alphabet(dict(
    A=Base(
        id="dAMP", name="2'-deoxyadenosine 5'-monophosphate(2−)",
        identifiers=(
            Identifier('chebi', 'CHEBI:58245'),
            Identifier('pubchem.compound', '22848660'),
            Identifier('metacyc.compound', 'DAMP'),
        ),
        structure=('InChI=1S/C10H14N5O6P/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(16)6(21-7)'
                   '2-20-22(17,18)19/h3-7,16H,1-2H2,(H2,11,12,13)(H2,17,18,19)/p-2/t5-,6+,7+/m0/s1')),
    C=Base(
        id="dCMP", name="2'-deoxycytidine 5'-monophosphate(2−)",
        identifiers=(
            Identifier('chebi', 'CHEBI:57566'),
            Identifier('pubchem.compound', '7058169'),
            Identifier('metacyc.compound', 'DCMP'),
        ),
        structure=('InChI=1S/C9H14N3O7P/c10-7-1-2-12(9(14)11-7)8-3-5(13)6(19-8)4-18-20(15,16)17/'
                   'h1-2,5-6,8,13H,3-4H2,(H2,10,11,14)(H2,15,16,17)/p-2/t5-,6+,8+/m0/s1')),
    G=Base(
        id="dGMP", name="2'-deoxyguanosine 5'-monophosphate(2−)",
        identifiers=(
            Identifier('chebi', 'CHEBI:57673'),
            Identifier('pubchem.compound', '6994968'),
            Identifier('metacyc.compound', 'DGMP'),
        ),
        structure=('InChI=1S/C10H14N5O7P/c11-10-13-8-7(9(17)14-10)12-3-15(8)6-1-4(16)5(22-6)'
                   '2-21-23(18,19)20/h3-6,16H,1-2H2,(H2,18,19,20)(H3,11,13,14,17)/p-2/t4-,5+,6+/m0/s1')),
    T=Base(
        id="dTMP", name="2'-deoxythymidine 5'-monophosphate(2−)",
        identifiers=(
            Identifier('chebi', 'CHEBI:63528'),
            Identifier('pubchem.compound', '16755631'),
            Identifier('metacyc.compound', 'TMP'),
        ),
        structure=('InChI=1S/C10H15N2O8P/c1-5-3-12(10(15)11-9(5)14)8-2-6(13)7(20-8)4-19-21'
                   '(16,17)18/h3,6-8,13H,2,4H2,1H3,(H,11,14,15)(H2,16,17,18)/p-2/t6-,7+,8+/m0/s1')),
))
# :obj:`Alphabet`: Alphabet for DNA nucleotides


class DnaForm(BpForm):
    """ DNA form """
    DEFAULT_ALPHABET = dna_alphabet
    DEFAULT_BOND_FORMULA = EmpiricalFormula('H') * -1
    DEFAULT_BOND_CHARGE = 1

    def __init__(self, base_seq=None):
        """
        Args:
            base_seq (:obj:`BaseSequence`, optional): bases of the DNA form
        """
        super(DnaForm, self).__init__(base_seq=base_seq)
