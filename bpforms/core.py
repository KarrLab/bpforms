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
        bond_formula (:obj:`EmpiricalFormula`): empirical formula for bonds between bases
        bond_mol_wt (:obj:`float`): molecular weight of bonds between bases
        bond_charge (:obj:`int`): charge of bonds between bases
    """

    def __init__(self, bases, bond_formula=None, bond_mol_wt=0., bond_charge=0):
        """
        Args:
            bases (:obj:`list` of :obj:`Base`): bases of the biopolymer
            bond_formula (:obj:`EmpiricalFormula`, optional): empirical formula for bonds between bases
            bond_mol_wt (:obj:`float`, optional): molecular weight of bonds between bases
            bond_charge (:obj:`int`, optional): charge of bonds between bases
        """
        self.bases = bases
        self.bond_formula = bond_formula or EmpiricalFormula()
        self.bond_mol_wt = bond_mol_wt
        self.bond_charge = bond_charge

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
            len += count
        return len

    def get_formula(self):
        """ Get the chemical formula

        Returns:
            :obj:`EmpiricalFormula`: chemical formula
        """
        formula = EmpiricalFormula()
        for base, count in self.get_base_counts().items():
            formula += count * base.get_formula()
        return formula + (self.get_length() - 1) * self.bond_formula

    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        mol_wt = 0.
        for base, count in self.get_base_counts().items():
            mol_wt += count * base.mol_wt()
        return mol_wt + (self.get_length() - 1) * self.bond_mol_wt

    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        charge = 0
        for base, count in self.get_base_counts().items():
            charge += count * base.get_charge()
        return charge + (self.get_length() - 1) * self.bond_charge

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
        super(DnaForm, self).__init__(structure, alphabet=DnaAlphabet,
                                      bond_formula=-EmpiricalFormula('HP2O7'),
                                      bond_mol_wt=-EmpiricalFormula('HP2O7').get_molecular_weight(),
                                      bond_charge=3)


class RnaForm(BpForm):
    """ RNA form """

    def __init__(self, structure):
        """
        Args:
            structure (:obj:`str`): structure of the biopolymer
        """
        super(RnaForm, self).__init__(structure, alphabet=RnaAlphabet,
                                      bond_formula=-EmpiricalFormula('HP2O7'),
                                      bond_mol_wt=-EmpiricalFormula('HP2O7').get_molecular_weight(),
                                      bond_charge=3)


class ProteinForm(BpForm):
    """ Protein form """

    def __init__(self, structure):
        """
        Args:
            structure (:obj:`str`): structure of the biopolymer
        """
        super(ProteinForm, self).__init__(structure, alphabet=ProteinAlphabet,
                                          bond_formula=-EmpiricalFormula('H2O'),
                                          bond_mol_wt=-EmpiricalFormula('H2O').get_molecular_weight(),
                                          bond_charge=0)


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

    def __eq__(self, other):
        """ Check if two identifiers are semantically equal

        Args:
            other (:obj:`Identifier`): another identifier

        Returns:
            :obj:`bool`: True, if the identifiers are semantically equal
        """
        return self is other or (self.__class__ == other.__class__ and self.ns == other.ns and self.id == other.id)

    def __hash__(self):
        """ Generate a hash

        Returns:
            :obj:`int`: hash
        """
        return hash((self.ns, self.id))


class Base(object):
    """ A base in a biopolymer

    Attributes:
        chars (:obj:`str`): string representation which must be one uppercase character, optionally followed by one or more lowercase characters
        id (:obj:`str`): id
        name (:obj:`str`): name
        synonyms (:obj:`set` of :obj:`str`): synonyms        
        identifiers (:obj:`set` of :obj:`Identifier`, optional): identifiers in namespaces for external databases
        structure (:obj:`openbabel.OBMol`): chemical structure
        mass_shift (:obj:`float`): additional mass (Dalton) relative to structure
        charge_shift (:obj:`float`): additional charge relative to structure
        start_position (:obj:`tuple`): uncertainty in location of base
        end_position (:obj:`tuple`): uncertainty in location of base
    """

    def __init__(self, chars='', id='', name='', synonyms=None, identifiers=None, structure=None,
                 mass_shift=0., charge_shift=0, start_position=None, end_position=None):
        """
        Attributes:
            chars (:obj:`str`, optional): string representation which must be one uppercase character, optionally followed by one or more lowercase characters
            id (:obj:`str`, optional): id
            name (:obj:`str`, optional): name
            synonyms (:obj:`set` of :obj:`str`, optional): synonyms
            identifiers (:obj:`set` of :obj:`Identifier`, optional): identifiers in namespaces for external databases
            structure (:obj:`openbabel.OBMol` or :obj:`str`, optional): chemical structure
            mass_shift (:obj:`float`, optional): additional mass (Dalton) relative to structure
            charge_shift (:obj:`float`, optional): additional charge relative to structure
            start_position (:obj:`tuple, optional`): uncertainty in location of base
            end_position (:obj:`tuple, optional`): uncertainty in location of base
        """
        if chars and not re.match(r'^[A-Z][a-z0-9_]*$', chars):
            raise ValueError('String representation must be one uppercase character, '
                             'optionally followed by one or more lowercase characters')
        self.chars = chars
        self.id = id
        self.name = name
        self.synonyms = synonyms or set()
        self.identifiers = identifiers or set()

        if not isinstance(structure, openbabel.OBMol):
            ob_mol = openbabel.OBMol()
            conversion = openbabel.OBConversion()
            assert conversion.SetInFormat('inchi'), 'Unable to set format to inchi'
            if not conversion.ReadString(ob_mol, structure):
                raise ValueError('Unable to parse structure for {}'.format(id))
            structure = ob_mol
        self.structure = structure

        self.mass_shift = mass_shift
        self.charge_shift = charge_shift
        self.start_position = start_position
        self.end_position = end_position

    def protonate(self, ph=7.):
        """ Update to the major protonation state at the pH 

        Args:
            ph (:obj:`float`, optional): pH
        """
        self.structure.CorrectForPH(ph)

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
        return self.structure.GetTotalCharge() + self.mass_shift

    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        return self.structure.GetMolWt() + self.charge_shift

    def __str__(self):
        """ Get a string representation of the base

        Returns:
            :obj:`str`: string representation of the base
        """
        els = []
        if self.chars:
            els.append('chars: ' + self.chars)
        if self.id:
            els.append('id: ' + self.id)
        if self.name:
            els.append('name: "' + self.name + '"')
        for synonym in self.synonyms:
            els.append('synonym: "' + synonym + '"')
        for identifier in self.identifiers:
            els.append('identifier: ' + identifier.ns + '/' + identifier.id)
        if self.mass_shift:
            els.append('mass-shift: ' + str(self.mass_shift))
        if self.charge_shift:
            els.append('charge-shift: ' + str(self.charge_shift))
        if self.start_position or self.end_position:
            els.append('position: {}-{}'.format(self.start_position or '', self.end_position or ''))

        return '[' + ' | '.join(els) + ']'

    def is_equal(self, other):
        """ Check if two bases are semantically equal

        Args:
            other (:obj:`Base'): another base

        Returns:
            :obj:`bool`: :obj:`True`, if the objects have the same structure
        """
        if self is other:
            return True
        if not self.__class__ == other.__class__:
            return False

        attrs = ['chars', 'id', 'name', 'synonyms', 'identifiers', 'mass_shift', 'charge_shift', 'start_position', 'end_position']
        for attr in attrs:
            if getattr(self, attr) != getattr(other, attr):
                return False

        if self.get_inchi() != other.get_inchi():
            return False


class Alphabet(object):
    """ Alphabet for bases 

    Attributes:
        parser (:obj:`lark.Lark`): parser
    """
    parser = lark.Lark(
        ''' ?start: seq

            seq: base+
            ?base: alphabet_base | inline_base
            alphabet_base: CHARS
            inline_base: "[" WS* inline_base_attr (ATTR_SEP inline_base_attr)* WS* "]"
            inline_base_attr: chars | id | name | synonym | identifier | structure | mass_shift | charge_shift | position
            chars: "chars" FIELD_SEP WS* CHARS
            id: "id" FIELD_SEP WS* ESCAPED_STRING
            name: "name" FIELD_SEP WS* ESCAPED_STRING
            synonym: "synonym" FIELD_SEP WS* ESCAPED_STRING
            identifier: "identifier" FIELD_SEP WS* identifier_ns "/" identifier_id
            identifier_ns: /[A-Za-z0-9_\.]+/
            identifier_id: /[A-Za-z0-9_\-]+/
            structure: "structure" FIELD_SEP WS* INCHI            
            mass_shift: "mass-shift" FIELD_SEP WS* DALTON
            charge_shift: "charge-shift" FIELD_SEP WS* CHARGE
            position: "position" FIELD_SEP WS* START_LOCATION WS* "-" WS* END_LOCATION
            ATTR_SEP: WS* "|" WS*
            FIELD_SEP: ":"
            CHARS: /[A-Z][a-z0-9_]*/
            INCHI: /InChI=1S\/[A-Za-z0-9\(\)\-\+,\/]+/
            DALTON: /[\-\+]?[0-9]+(\.[0-9]*)?/
            CHARGE: /[\-\+]?[0-9]+/
            START_LOCATION: INT
            END_LOCATION: INT
            %import common.WS
            %import common.ESCAPED_STRING
            %import common.INT
            ''')

    class ParseTreeTransformer(lark.Transformer):
        @lark.v_args(inline=True)
        def alphabet_base(self, s):
            return Base(s[0].__str__())

        @lark.v_args(inline=True)
        def inline_base(self, s):
            return Base(s[0].__str__())
    parse_tree_transformer = ParseTreeTransformer()

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
        tree = cls.parser.parse(str)
        return cls.parse_tree_transformer.transform(tree)


class DnaAlphabet(Alphabet):
    """ Alphabet for DNA nucleotides """
    A = Base(
        chars="A", id="dAMP", name="2'-deoxyadenosine 5'-monophosphate(2−)",
        identifiers=(
            Identifier('chebi', 'CHEBI:58245'),
            Identifier('pubchem.compound', '22848660'),
            Identifier('metacyc.compound', 'DAMP'),
        ),
        structure='InChI=1S/C10H14N5O6P/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(16)6(21-7)2-20-22(17,18)19/h3-7,16H,1-2H2,(H2,11,12,13)(H2,17,18,19)/p-2/t5-,6+,7+/m0/s1')
    C = Base(
        chars="C", id="dCMP", name="2'-deoxyadenosine 5'-monophosphate(2−)",
        identifiers=(
            Identifier('chebi', 'CHEBI:57566'),
            Identifier('pubchem.compound', '7058169'),
            Identifier('metacyc.compound', 'DCMP'),
        ),
        structure='InChI=1S/C9H14N3O7P/c10-7-1-2-12(9(14)11-7)8-3-5(13)6(19-8)4-18-20(15,16)17/h1-2,5-6,8,13H,3-4H2,(H2,10,11,14)(H2,15,16,17)/p-2/t5-,6+,8+/m0/s1')
    G = Base(
        chars="G", id="dGMP", name="2'-deoxyadenosine 5'-monophosphate(2−)",
        identifiers=(
            Identifier('chebi', 'CHEBI:57673'),
            Identifier('pubchem.compound', '6994968'),
            Identifier('metacyc.compound', 'DGMP'),
        ),
        structure='InChI=1S/C10H14N5O7P/c11-10-13-8-7(9(17)14-10)12-3-15(8)6-1-4(16)5(22-6)2-21-23(18,19)20/h3-6,16H,1-2H2,(H2,18,19,20)(H3,11,13,14,17)/p-2/t4-,5+,6+/m0/s1')
    T = Base(
        chars="T", id="dTMP", name="2'-deoxyadenosine 5'-monophosphate(2−)",
        identifiers=(
            Identifier('chebi', 'CHEBI:63528'),
            Identifier('pubchem.compound', '16755631'),
            Identifier('metacyc.compound', 'TMP'),
        ),
        structure='InChI=1S/C10H15N2O8P/c1-5-3-12(10(15)11-9(5)14)8-2-6(13)7(20-8)4-19-21(16,17)18/h3,6-8,13H,2,4H2,1H3,(H,11,14,15)(H2,16,17,18)/p-2/t6-,7+,8+/m0/s1')


class RnaAlphabet(Alphabet):
    """ Alphabet for RNA nucleotides """
    pass


class ProteinAlphabet(Alphabet):
    """ Alphabet for amino acids """
    pass
