""" Alphabet and BpForm to represent modified RNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from .core import Alphabet, Base, BpForm, Identifier
from wc_utils.util.chem import EmpiricalFormula

rna_alphabet = Alphabet(dict(
    A=Base(
        id="AMP", name="2'-adenosine 5'-monophosphate(2−)",
        identifiers=(
            Identifier('chebi', 'CHEBI:456215'),
            Identifier('pubchem.compound', '15938965'),
            Identifier('metacyc.compound', 'AMP'),
        ),
        structure=('InChI=1S/C10H14N5O7P/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(22-10)'
                   '1-21-23(18,19)20/h2-4,6-7,10,16-17H,1H2,(H2,11,12,13)(H2,18,19,20)/p-2/t4-,6-,7-,10-/m1/s1')),
    C=Base(
        id="CMP", name="2'-cytidine 5'-monophosphate(2−)",
        identifiers=(
            Identifier('chebi', 'CHEBI:60377'),
            Identifier('pubchem.compound', '7058165'),
            Identifier('metacyc.compound', 'CMP'),
        ),
        structure=('InChI=1S/C9H14N3O8P/c10-5-1-2-12(9(15)11-5)8-7(14)6(13)4(20-8)3-19-21(16,17)18/'
                   'h1-2,4,6-8,13-14H,3H2,(H2,10,11,15)(H2,16,17,18)/p-2/t4-,6-,7-,8-/m1/s1')),
    G=Base(
        id="GMP", name="2'-guanosine 5'-monophosphate(2−)",
        identifiers=(
            Identifier('chebi', 'CHEBI:58115'),
            Identifier('pubchem.compound', '1807035'),
            Identifier('metacyc.compound', 'GMP'),
        ),
        structure=('InChI=1S/C10H14N5O8P/c11-10-13-7-4(8(18)14-10)12-2-15(7)9-6(17)5(16)3(23-9)1-22-24(19,20)21/'
                   'h2-3,5-6,9,16-17H,1H2,(H2,19,20,21)(H3,11,13,14,18)/p-2/t3-,5-,6-,9-/m1/s1')),
    U=Base(
        id="UMP", name="2'-uridine 5'-monophosphate(2−)",
        identifiers=(
            Identifier('chebi', 'CHEBI:57865'),
            Identifier('pubchem.compound', '1778309'),
            Identifier('metacyc.compound', 'UMP'),
        ),
        structure=('InChI=1S/C9H13N2O9P/c12-5-1-2-11(9(15)10-5)8-7(14)6(13)4(20-8)3-19-21(16,17)18/'
                   'h1-2,4,6-8,13-14H,3H2,(H,10,12,15)(H2,16,17,18)/p-2/t4-,6-,7-,8-/m1/s1')),
))
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
