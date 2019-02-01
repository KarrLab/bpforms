""" bpforms

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_utils.util.chem import EmpiricalFormula
import Bio
import itertools
import lark
import openbabel
import re


class BpForm(object):
    """ Biopolymer form 

    Attributes:
        bases (:obj:`list` of :obj:`Base`): bases of the biopolymer
    """

    def __init__(self, bases):
        """
        Args:
            bases (:obj:`list` of :obj:`Base`): bases of the biopolymer
        """
        self.bases = bases

    def get_base_counts(self):
        """ Get the frequency of each base within the biopolymer

        Returns:
            :obj:`dict`: dictionary that maps bases to their counts
        """
        counts = {}
        for base, group in itertools.groupby(self.bases):
            counts[base] = len(list(group))
        return counts

    def protonate(self, ph=7.):
        """ Update to the major protonation state of each modification at the pH 

        Args:
            ph (:obj:`float`, optional): pH
        """
        for base in set(self.bases):
            base.protonate(ph=ph)

    def get_length(self):
        """ Get the length

        Returns:
            :obj:`int`: length
        """
        len = 0
        for base, count in self.get_base_counts().items():
            len += count * base.get_length()
        return len

    def get_formula(self):
        """ Get the chemical formula

        Returns:
            :obj:`EmpiricalFormula`: chemical formula
        """
        formula = EmpiricalFormula()
        for base, count in self.get_base_counts().items():
            formula += count * base.get_formula()
        return formula

    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        mol_wt = 0.
        for base, count in self.get_base_counts().items():
            mol_wt += count * base.mol_wt()
        return mol_wt

    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        charge = 0
        for base, count in self.get_base_counts().items():
            charge += count * base.get_charge()
        return charge

    def is_equal(self, other):
        """ Check if two biopolymer forms are semantically equal

        Args:
            other (:obj:`BpForm'): another biopolymer form

        Returns:
            :obj:`bool`: :obj:`True`, if the objects have the same structure
        """
        return self is other or (self.__class__ == other.__class__ and self.structure == other.structure)

    def __str__(self):
        """ Get a string representation of the biopolymer form

        Returns:
            :obj:`str`: string representation of the biopolymer form
        """
        return ''.join(str(base) for base in bases)


class DnaForm(BpForm):
    """ DNA form """

    def __init__(self, structure):
        """
        Args:
            structure (:obj:`str`): structure of the biopolymer
        """
        super(DnaForm, self).__init__(structure, unmodified_alphabet=UnmodifiedDnaAlphabet)


class RnaForm(BpForm):
    """ RNA form """

    def __init__(self, structure):
        """
        Args:
            structure (:obj:`str`): structure of the biopolymer
        """
        super(RnaForm, self).__init__(structure, unmodified_alphabet=UnmodifiedRnaAlphabet)


class ProteinForm(BpForm):
    """ Protein form """

    def __init__(self, structure):
        """
        Args:
            structure (:obj:`str`): structure of the biopolymer
        """
        super(ProteinForm, self).__init__(structure, unmodified_alphabet=UnmodifiedProteinAlphabet)


def get_form(type):
    if type == 'dna':
        return DnaForm
    if type == 'rna':
        return RnaForm
    if type == 'protein':
        return ProteinForm

    raise ValueError('Type "{}" must be "dna", "rna", or "protein"'.format(type))


class Identifier(object):
    """ A identifier in a namespace for an external database

    Attributes:
        ns (:obj:`str`): namespace
        id (:obj:`str`): id in namespace
    """

    def __init__(self, ns, id):
        """
        Args:
            ns (:obj:`str`): namespace
            id (:obj:`str`): id in namespace
        """
        self.ns = ns
        self.id = id


class Base(object):
    """ A base in a biopolymer

    Attributes:
        chars (:obj:`str`): string representation which must be one uppercase character, optionally followed by one or more lowercase characters
        id (:obj:`str`): id
        name (:obj:`str`): name
        synonyms (:obj:`list` of :obj:`str`): synonyms
        structure (:obj:`openbabel.OBMol`): chemical structure
        identifiers (:obj:`tuple` of :obj:`Identifier`, optional): identifiers in namespaces for external databases
    """

    def __init__(self, chars='', id='', name='', synonyms=None, structure=None, identifiers=None):
        """
        Attributes:
            chars (:obj:`str`, optional): string representation which must be one uppercase character, optionally followed by one or more lowercase characters
            id (:obj:`str`, optional): id
            name (:obj:`str`, optional): name
            synonyms (:obj:`list` of :obj:`str`, optional): synonyms
            structure (:obj:`openbabel.OBMol` or :obj:`str`, optional): chemical structure
            identifiers (:obj:`tuple` of :obj:`Identifier`, optional): identifiers in namespaces for external databases
        """
        if chars and not re.match(r'^[A-Z][a-z0-9_]*$', chars):
            raise ValueError('String representation must be one uppercase character, '
                             'optionally followed by one or more lowercase characters')
        self.chars = chars
        self.id = id
        self.name = name
        self.synonyms = synonyms or []

        if not isinstance(structure, openbabel.OBMol):
            ob_mol = openbabel.OBMol()
            conversion = openbabel.OBConversion()
            assert conversion.SetInFormat('inchi'), 'Unable to set format to inchi'
            if not conversion.ReadString(ob_mol, structure):
                raise ValueError('Unable to parse structure for {}'.format(id))
            structure = ob_mol
        self.structure = structure

        self.identifiers = identifiers or ()

    def protonate(self, ph=7.):
        """ Update to the major protonation state at the pH 

        Args:
            ph (:obj:`float`, optional): pH
        """
        self.structure.CorrectForPH(ph)

    def get_length(self):
        """ Get the length

        Returns:
            :obj:`int`: length
        """
        return 1

    def get_inchi(self):
        """ Get InChI representration of structure

        Returns:
            :obj:`str`: InChI representration of structure
        """
        conversion = openbabel.OBConversion()
        assert conversion.SetOutFormat('inchi'), 'Unable to set format to inchi'
        return conversion.WriteString(ob_mol)

    def get_formula(self):
        """ Get the chemical formula

        Returns:
            :obj:`EmpiricalFormula`: chemical formula
        """
        return EmpiricalFormula(self.structure.GetFormula().rstrip('+-'))

    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        return self.structure.GetTotalCharge()

    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        return self.structure.GetMolWt()

    def __str__(self):
        """ Get a string representation of the base

        Returns:
            :obj:`str`: string representation of the base
        """
        if self.chars:
            return self.chars
        else:
            els = []
            if self.id:
                els.append('id: ' + self.id)
            if self.name:
                els.append('name: ' + self.name)
            for synonym in self.synonyms:
                els.append('synonym: ' + synonym)
            for identifier in self.identifiers:
                els.append('identifier: ' + identifier.ns + '/' + identifier.id)
            return '[' + ' | '.join(els) + ']'

    def is_equal(self, other):
        """ Check if two bases are semantically equal

        Args:
            other (:obj:`Base'): another base

        Returns:
            :obj:`bool`: :obj:`True`, if the objects have the same structure
        """
        return self is other or (self.__class__ == other.__class__ and self.get_inchi() == other.get_inchi())


class Alphabet(object):
    """ Alphabet for unmodified bases 

    Attributes:
        parser (:obj:`lark.Lark`): parser
    """
    parser = lark.Lark(
        ''' ?start: seq

            seq: base+
            ?base: alphabet_base | inline_base
            alphabet_base: CHARS
            inline_base: "[" WS* inline_base_attr (ATTR_SEP inline_base_attr)* WS* "]"
            inline_base_attr: chars | id | name | synonym | identifier | structure | position | mass_shift
            chars: "chars" FIELD_SEP WS* CHARS
            id: "id" FIELD_SEP WS* ESCAPED_STRING
            name: "name" FIELD_SEP WS* ESCAPED_STRING
            synonym: "synonym" FIELD_SEP WS* ESCAPED_STRING
            identifier: "identifier" FIELD_SEP WS* identifier_ns "/" identifier_id
            identifier_ns: /[A-Za-z0-9_\.]+/
            identifier_id: /[A-Za-z0-9_\-]+/
            structure: "structure" FIELD_SEP WS* INCHI
            position: "position" FIELD_SEP WS* LOCATION WS* "-" WS* LOCATION
            mass_shift: "mass-shift" FIELD_SEP WS* DALTON
            ATTR_SEP: WS* "|" WS*
            FIELD_SEP: ":"
            CHARS: /[A-Z][a-z0-9_]*/
            INCHI: /InChI=1S\/[A-Za-z0-9\(\)\-\+,\/]+/
            DALTON: /[\-\+]?[0-9]+(\.[0-9]*)?/
            LOCATION: /[0-9]+/
            %import common.WS
            %import common.ESCAPED_STRING
            ''')

    @classmethod
    def get_bases(cls):
        """ Get bases of alphabet

        Returns:
            :obj:`list` of :obj:`Bases`: bases of alphabet
        """
        bases = []
        for attr_name in dir(cls):
            if not attr_name.startswith('_'):
                attr = getattr(cls, attr_name)
                if isinstance(attr, Base):
                    bases.append(attr)
        return bases

    @classmethod
    def create_bpform(cls, str):
        """ Create biopolymer form from string representation

        Args:
            str (:obj:`str`): string representation of biopolymer

        Returns:
            :obj:`BpForm`: biopolymer form
        """
        cls.parser.parse(str)
        return cls(bases)


class UnmodifiedDnaAlphabet(Alphabet):
    """ Alphabet for unmodified DNA nucleotides """
    A = Base(
        chars="A", id="dAMP", name="2'-deoxyadenosine 5'-monophosphate(2−)",
        structure='InChI=1S/C10H14N5O6P/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(16)6(21-7)2-20-22(17,18)19/h3-7,16H,1-2H2,(H2,11,12,13)(H2,17,18,19)/p-2/t5-,6+,7+/m0/s1',
        identifiers=(
            Identifier('chebi', 'CHEBI:58245'),
            Identifier('pubchem.compound', '22848660'),
            Identifier('metacyc.compound', 'DAMP'),
        ))
    C = Base(
        chars="C", id="dCMP", name="2'-deoxyadenosine 5'-monophosphate(2−)",
        structure='InChI=1S/C9H14N3O7P/c10-7-1-2-12(9(14)11-7)8-3-5(13)6(19-8)4-18-20(15,16)17/h1-2,5-6,8,13H,3-4H2,(H2,10,11,14)(H2,15,16,17)/p-2/t5-,6+,8+/m0/s1',
        identifiers=(
            Identifier('chebi', 'CHEBI:57566'),
            Identifier('pubchem.compound', '7058169'),
            Identifier('metacyc.compound', 'DCMP'),
        ))
    G = Base(
        chars="G", id="dGMP", name="2'-deoxyadenosine 5'-monophosphate(2−)",
        structure='InChI=1S/C10H14N5O7P/c11-10-13-8-7(9(17)14-10)12-3-15(8)6-1-4(16)5(22-6)2-21-23(18,19)20/h3-6,16H,1-2H2,(H2,18,19,20)(H3,11,13,14,17)/p-2/t4-,5+,6+/m0/s1',
        identifiers=(
            Identifier('chebi', 'CHEBI:57673'),
            Identifier('pubchem.compound', '6994968'),
            Identifier('metacyc.compound', 'DGMP'),
        ))
    T = Base(
        chars="T", id="dTMP", name="2'-deoxyadenosine 5'-monophosphate(2−)",
        structure='InChI=1S/C10H15N2O8P/c1-5-3-12(10(15)11-9(5)14)8-2-6(13)7(20-8)4-19-21(16,17)18/h3,6-8,13H,2,4H2,1H3,(H,11,14,15)(H2,16,17,18)/p-2/t6-,7+,8+/m0/s1',
        identifiers=(
            Identifier('chebi', 'CHEBI:63528'),
            Identifier('pubchem.compound', '16755631'),
            Identifier('metacyc.compound', 'TMP'),
        ))


class UnmodifiedRnaAlphabet(Alphabet):
    """ Alphabet for unmodified RNA nucleotides """
    pass


class UnmodifiedProteinAlphabet(Alphabet):
    """ Alphabet for unmodified amino acids """
    pass
