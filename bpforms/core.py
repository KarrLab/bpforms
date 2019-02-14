""" Classes to represent modified forms of DNA, RNA, and proteins

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from ruamel import yaml
from wc_utils.util.chem import EmpiricalFormula, Protonator
import abc
import attrdict
import lark
import openbabel
import pkg_resources
import re
import subprocess
import time


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

    @property
    def ns(self):
        """ Get the namespace

        Returns:
            :obj:`str`: namespace
        """
        return self._ns

    @ns.setter
    def ns(self, value):
        """ Set the namespace

        Args:
            value (:obj:`str`): namespace

        Raises:
            :obj:`ValueError`: if the namespace is empty
        """
        if not value:
            raise ValueError('`ns` cannot be empty')
        self._ns = value

    @property
    def id(self):
        """ Get the id

        Returns:
            :obj:`str`: id
        """
        return self._id

    @id.setter
    def id(self, value):
        """ Set the id

        Args:
            value (:obj:`str`): id

        Raises:
            :obj:`ValueError`: if the id is empty
        """
        if not value:
            raise ValueError('`id` cannot be empty')
        self._id = value

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


class IdentifierSet(set):
    """ Set of identifiers """

    def __init__(self, identifiers=None):
        """
        Args:
            identifiers (:obj:iterable of :obj:`Identifier`): iterable of identifiers
        """
        super(IdentifierSet, self).__init__()
        if identifiers is not None:
            for identifier in identifiers:
                self.add(identifier)

    def add(self, identifier):
        """ Add an identifier

        Args:
            identifier (:obj:`Identifier`): identifier

        Raises:
            :obj:`ValueError`: if the `identifier` is not an instance of `Indentifier`
        """
        if not isinstance(identifier, Identifier):
            raise ValueError('`identifier` must be an instance of `Indentifier`')
        super(IdentifierSet, self).add(identifier)

    def update(self, identifiers):
        """ Add a set of identifiers

        Args:
            identifiers (iterable of :obj:`Identifier`): identifiers
        """
        for identifier in identifiers:
            self.add(identifier)

    def symmetric_difference_update(self, other):
        """ Remove common elements with other and add elements from other not in self

        Args:
            other (:obj:`IdentifierSet`): other set of identifiers
        """
        if not isinstance(other, IdentifierSet):
            other = IdentifierSet(other)
        super(IdentifierSet, self).symmetric_difference_update(other)


class SynonymSet(set):
    """ Set of synonyms """

    def __init__(self, synonyms=None):
        """
        Args:
            synonyms (:obj:iterable of :obj:`str`): iterable of synonyms

        Raises:
            :obj:`ValueError`: if synonyms is not an iterable of string
        """
        super(SynonymSet, self).__init__()
        if isinstance(synonyms, str):
            raise ValueError('synonyms must be an iterable of strings')
        if synonyms is not None:
            for synonym in synonyms:
                self.add(synonym)

    def add(self, synonym):
        """ Add an synonym

        Args:
            synonym (:obj:`str`): synonym

        Raises:
            :obj:`ValueError`: if the `synonym` is not an instance of `Indentifier`
        """
        if not synonym or not isinstance(synonym, str):
            raise ValueError('`synonyms` must be a non-empty string')
        super(SynonymSet, self).add(synonym)

    def update(self, synonyms):
        """ Add a set of synonyms

        Args:
            synonyms (iterable of :obj:`SynonymSet`): synonyms
        """
        for synonym in synonyms:
            self.add(synonym)

    def symmetric_difference_update(self, other):
        """ Remove common synonyms with other and add synonyms from other not in self

        Args:
            other (:obj:`SynonymSet`): other set of synonyms
        """
        if not isinstance(other, SynonymSet):
            other = SynonymSet(other)
        super(SynonymSet, self).symmetric_difference_update(other)


class Monomer(object):
    """ A monomer in a biopolymer

    Attributes:
        id (:obj:`str`): id
        name (:obj:`str`): name
        synonyms (:obj:`set` of :obj:`str`): synonyms
        identifiers (:obj:`set` of :obj:`Identifier`, optional): identifiers in namespaces for external databases
        structure (:obj:`openbabel.OBMol`): chemical structure
        delta_mass (:obj:`float`): additional mass (Dalton) relative to structure
        delta_charge (:obj:`int`): additional charge relative to structure
        start_position (:obj:`tuple`): uncertainty in the location of the monomer
        end_position (:obj:`tuple`): uncertainty in the location of the monomer
        base_monomers (:obj:`set` of :obj:`Monomer`): monomers which this monomer is derived from
        comments (:obj:`str`): comments
    """

    def __init__(self, id=None, name=None, synonyms=None, identifiers=None, structure=None,
                 delta_mass=None, delta_charge=None, start_position=None, end_position=None,
                 base_monomers=None, comments=None):
        """
        Attributes:
            id (:obj:`str`, optional): id
            name (:obj:`str`, optional): name
            synonyms (:obj:`set` of :obj:`str`, optional): synonyms
            identifiers (:obj:`set` of :obj:`Identifier`, optional): identifiers in namespaces for external databases
            structure (:obj:`openbabel.OBMol` or :obj:`str`, optional): chemical structure
            delta_mass (:obj:`float`, optional): additional mass (Dalton) relative to structure
            delta_charge (:obj:`float`, optional): additional charge relative to structure
            start_position (:obj:`int`, optional): uncertainty in the location of the monomer
            end_position (:obj:`int`, optional): uncertainty in the location of the monomer
            base_monomers (:obj:`set` of :obj:`Monomer`, optional): monomers which this monomer is derived from
            comments (:obj:`str`, optional): comments
        """
        self.id = id
        self.name = name
        self.synonyms = synonyms or SynonymSet()
        self.identifiers = identifiers or IdentifierSet()
        self.structure = structure
        self.delta_mass = delta_mass
        self.delta_charge = delta_charge
        self.start_position = start_position
        self.end_position = end_position
        self.base_monomers = base_monomers or set()
        self.comments = comments

    @property
    def id(self):
        """ Get id

        Returns:
            :obj:`str`: id
        """
        return self._id

    @id.setter
    def id(self, value):
        """ Set id

        Args:
            value (:obj:`str`): id

        Raises:
            :obj:`ValueError`: if `value` is not a a string or None
        """
        if value and not isinstance(value, str):
            raise ValueError('`id` must be a string or None')
        self._id = value

    @property
    def name(self):
        """ Get name

        Returns:
            :obj:`str`: name
        """
        return self._name

    @name.setter
    def name(self, value):
        """ Set name

        Args:
            value (:obj:`str`): name

        Raises:
            :obj:`ValueError`: if `value` is not a a string or None
        """
        if value and not isinstance(value, str):
            raise ValueError('`name` must be a string or None')
        self._name = value

    @property
    def synonyms(self):
        """ Get synonyms

        Returns:
            :obj:`SynonymSet`: synonyms
        """
        return self._synonyms

    @synonyms.setter
    def synonyms(self, value):
        """ Set synonyms

        Args:
            value (:obj:`SynonymSet`): synonyms

        Raises:
            :obj:`ValueError`: if `synonyms` is not an instance of `SynonymSet`
        """
        if value is None:
            raise ValueError('`synonyms` must be an instance `SynonymSet`')
        if not isinstance(value, SynonymSet):
            value = SynonymSet(value)
        self._synonyms = value

    @property
    def identifiers(self):
        """ Get identifiers

        Returns:
            :obj:`IdentifierSet`: identifiers
        """
        return self._identifiers

    @identifiers.setter
    def identifiers(self, value):
        """ Set identifiers

        Args:
            value (:obj:`IdentifierSet`): identifiers

        Raises:
            :obj:`ValueError`: if `identifiers` is not an instance of `Indentifiers`
        """
        if value is None:
            raise ValueError('`identifiers` must be an instance `Indentifiers`')
        if not isinstance(value, IdentifierSet):
            value = IdentifierSet(value)
        self._identifiers = value

    @property
    def structure(self):
        """ Get structure

        Returns:
            :obj:`openbabel.OBMol`: structure
        """
        return self._structure

    @structure.setter
    def structure(self, value):
        """ Set structure

        Args:
            value (:obj:`openbabel.OBMol` or :obj:`str`): OpenBabel molecule, InChI-encoded structure, or None

        Raises:
            :obj:`ValueError`: if value is not an OpenBabel molecule, InChI-encoded structure, or None
        """
        if value and not isinstance(value, openbabel.OBMol):
            ob_mol = openbabel.OBMol()
            conversion = openbabel.OBConversion()
            assert conversion.SetInFormat('inchi'), 'Unable to set format to InChI'
            if not conversion.ReadString(ob_mol, value):
                raise ValueError('`structure` must be an OpenBabel molecule, InChI-encoded structure, or None')
            value = ob_mol

        self._structure = value or None

    @property
    def delta_mass(self):
        """ Get extra mass

        Returns:
            :obj:`float`: extra mass
        """
        return self._delta_mass

    @delta_mass.setter
    def delta_mass(self, value):
        """ Set extra mass

        Args:
            value (:obj:`float`): extra mass

        Raises:
            :obj:`ValueError`: if value is not a float or None
        """
        if value is not None:
            if not isinstance(value, (int, float)):
                raise ValueError('`delta_mass` must be a float or None')
            value = float(value)
        self._delta_mass = value

    @property
    def delta_charge(self):
        """ Get extra charge

        Returns:
            :obj:`int`: extra charge
        """
        return self._delta_charge

    @delta_charge.setter
    def delta_charge(self, value):
        """ Set extra charge

        Args:
            value (:obj:`int`): extra charge

        Raises:
            :obj:`ValueError`: if value is not an int or None
        """
        if value is not None:
            if not isinstance(value, (int, float)):
                raise ValueError('`delta_charge` must be an integer or None')
            if value != int(value):
                raise ValueError('`delta_charge` must be an integer or None')
            value = int(value)
        self._delta_charge = value

    @property
    def start_position(self):
        """ Get start position

        Returns:
            :obj:`int`: start position
        """
        return self._start_position

    @start_position.setter
    def start_position(self, value):
        """ Set start position

        Args:
            value (:obj:`float`): start position

        Raises:
            :obj:`ValueError`: if value is not an int or None
        """
        if value is not None:
            if not isinstance(value, (int, float)):
                raise ValueError('`start_position` must be a positive integer or None')
            if value != int(value) or value < 1:
                raise ValueError('`start_position` must be a positive integer or None')
            value = int(value)
        self._start_position = value

    @property
    def end_position(self):
        """ Get end position

        Returns:
            :obj:`int`: end position
        """
        return self._end_position

    @end_position.setter
    def end_position(self, value):
        """ Set end position

        Args:
            value (:obj:`float`): end position

        Raises:
            :obj:`ValueError`: if value is not an int or None
        """
        if value is not None:
            if not isinstance(value, (int, float)):
                raise ValueError('`end_position` must be a positive integer or None')
            if value != int(value) or value < 1:
                raise ValueError('`end_position` must be a positive integer or None')
            value = int(value)
        self._end_position = value

    @property
    def base_monomers(self):
        """ Get base monomers

        Returns:
            :obj:`set` of :obj:`Monomer`: base monomers
        """
        return self._base_monomers

    @base_monomers.setter
    def base_monomers(self, value):
        """ Set base monomers

        Args:
            value (:obj:`set` of :obj:`Monomer`): base monomers

        Raises:
            :obj:`ValueError`: if value is not an instance of :obj:`set`
        """
        if isinstance(value, list):
            value = set(value)
        if not isinstance(value, set):
            raise ValueError('`base_monomers` must be an instance of `set`')
        self._base_monomers = value

    @property
    def comments(self):
        """ Get comments

        Returns:
            :obj:`str`: comments
        """
        return self._comments

    @comments.setter
    def comments(self, value):
        """ Set comments

        Args:
            value (:obj:`str`): comments

        Raises:
            :obj:`ValueError`: if value is not a str or None
        """
        if value and not isinstance(value, str):
            raise ValueError('`comments` must be a string or None')
        self._comments = value

    def get_root_monomers(self):
        """ Get root monomers

        Returns:
            :obj:`set` of :obj:`Monomer`: root monomers            
        """
        if not self.base_monomers:
            return set([self])

        roots = set()
        for base_monomer in self.base_monomers:
            roots.update(base_monomer.get_root_monomers())

        return roots

    def protonate(self, ph):
        """ Update to the major protonation state at the pH

        Args:
            ph (:obj:`float`): pH
        """
        if self.structure:
            self.structure = Protonator.run(self.get_inchi(), ph=ph)

    def get_inchi(self):
        """ Get InChI representration of structure

        Returns:
            :obj:`str`: InChI representration of structure
        """
        if self.structure:
            conversion = openbabel.OBConversion()
            assert conversion.SetOutFormat('inchi'), 'Unable to set format to InChI'
            conversion.SetOptions('r', conversion.OUTOPTIONS)
            conversion.SetOptions('F', conversion.OUTOPTIONS)
            inchi = conversion.WriteString(self.structure).strip()
            inchi = inchi.replace('InChI=1/', 'InChI=1S/')
            i_fixed_h = inchi.find('/f')
            if i_fixed_h >= 0:
                inchi = inchi[0:i_fixed_h]
            return inchi
        else:
            return None

    IMAGE_URL_PATTERN = ('https://cactus.nci.nih.gov/chemical/structure/{}/image'
                         '?format=png'
                         '&bgcolor=transparent'
                         '&antialiasing=0')

    def get_image_url(self):
        """ Get URL for image of structure

        Returns:
            :obj:`str`: URL for image of structure
        """
        inchi = self.get_inchi()
        if inchi:
            inchi, _, _ = inchi.partition('/t')
            return self.IMAGE_URL_PATTERN.format(inchi)
        return None

    def get_formula(self):
        """ Get the chemical formula

        Returns:
            :obj:`EmpiricalFormula`: chemical formula
        """
        if not self.structure:
            raise ValueError('A structure must be defined to calculate the formula')

        el_table = openbabel.OBElementTable()
        formula = {}
        mass = 0
        for i_atom in range(self.structure.NumAtoms()):
            atom = self.structure.GetAtom(i_atom + 1)
            el = el_table.GetSymbol(atom.GetAtomicNum())
            if el in formula:
                formula[el] += 1
            else:
                formula[el] = 1
            mass += el_table.GetMass(atom.GetAtomicNum())
        formula = EmpiricalFormula(formula)

        # calc hydrogens because OpenBabel doesn't output this
        formula['H'] = round((self.structure.GetMolWt() - mass) / el_table.GetMass(1))
        return formula

    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        if self.structure:
            return self.get_formula().get_molecular_weight() + (self.delta_mass or 0.)
        return None

    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        if not self.structure:
            raise ValueError('A structure must be defined to calculate the charge')
        return self.structure.GetTotalCharge() + (self.delta_charge or 0)

    def to_dict(self, alphabet=None):
        """ Get a dictionary representation of the monomer

        Args:
            alphabet (:obj:`Alphabet`, optional): alphabet

        Returns:
            :obj:`dict`: dictionary representation of the monomer
        """
        dict = {}

        attrs = ['id', 'name', 'delta_mass', 'delta_charge', 'start_position', 'end_position', 'comments']
        for attr in attrs:
            val = getattr(self, attr)
            if val is not None:
                dict[attr] = val

        if self.synonyms:
            dict['synonyms'] = list(self.synonyms)

        if self.identifiers:
            dict['identifiers'] = [{'ns': i.ns, 'id': i.id} for i in self.identifiers]

        if self.structure:
            dict['structure'] = self.get_inchi()

        if self.base_monomers and alphabet:
            dict['base_monomers'] = []
            for monomer in self.base_monomers:
                monomer_code = alphabet.get_monomer_code(monomer)
                dict['base_monomers'].append(monomer_code)

        return dict

    def from_dict(self, dict, alphabet=None):
        """ Get a dictionary representation of the monomer

        Args:
            dict (:obj:`dict`): dictionary representation of the monomer
            alphabet (:obj:`Alphabet`, optional): alphabet

        Returns:
            :obj:`Monomer`: monomer
        """
        self.id = None
        self.name = None
        self.synonyms.clear()
        self.identifiers.clear()
        self.structure = None
        self.delta_mass = None
        self.delta_charge = None
        self.start_position = None
        self.end_position = None
        self.base_monomers.clear()
        self.comments = None

        attrs = ['id', 'name', 'delta_mass', 'delta_charge', 'start_position', 'end_position', 'comments']
        for attr in attrs:
            val = dict.get(attr, None)
            if val is not None:
                setattr(self, attr, val)

        synonyms = dict.get('synonyms', [])
        if synonyms:
            self.synonyms = SynonymSet(synonyms)

        identifiers = dict.get('identifiers', [])
        if identifiers:
            self.identifiers = IdentifierSet([Identifier(i['ns'], i['id']) for i in identifiers])

        structure = dict.get('structure', None)
        if structure:
            self.structure = structure

        base_monomer_ids = dict.get('base_monomers', [])
        if base_monomer_ids and alphabet:
            self.base_monomers = set([alphabet.monomers.get(monomer_id) for monomer_id in base_monomer_ids])

        return self

    def __str__(self, alphabet=None):
        """ Get a string representation of the monomer

        Args:
            alphabet (:obj:`Alphabet`, optional): alphabet

        Returns:
            :obj:`str`: string representation of the monomer
        """
        els = []
        if self.id:
            els.append('id: "' + self.id + '"')
        if self.name:
            els.append('name: "' + self.name.replace('"', '\\"') + '"')
        for synonym in self.synonyms:
            els.append('synonym: "' + synonym.replace('"', '\\"') + '"')
        for identifier in self.identifiers:
            els.append('identifier: "' + identifier.ns + '" / "' + identifier.id + '"')
        if self.structure:
            els.append('structure: ' + self.get_inchi())
        if self.delta_mass is not None:
            els.append('delta-mass: ' + str(self.delta_mass))
        if self.delta_charge is not None:
            els.append('delta-charge: ' + str(self.delta_charge))
        if self.start_position is not None or self.end_position is not None:
            els.append('position: {}-{}'.format(self.start_position or '', self.end_position or ''))
        if alphabet:
            for base_monomer in self.base_monomers:
                els.append('base-monomer: "{}"'.format(alphabet.get_monomer_code(base_monomer)))
        if self.comments:
            els.append('comments: "' + self.comments.replace('"', '\\"') + '"')

        return '[' + ' | '.join(els) + ']'

    def is_equal(self, other):
        """ Check if two monomers are semantically equal

        Args:
            other (:obj:`Monomer`): another monomer

        Returns:
            :obj:`bool`: :obj:`True`, if the objects have the same structure
        """
        if self is other:
            return True
        if self.__class__ != other.__class__:
            return False

        attrs = ['id', 'name', 'synonyms', 'identifiers',
                 'delta_mass', 'delta_charge', 'start_position', 'end_position',
                 'comments']
        for attr in attrs:
            if getattr(self, attr) != getattr(other, attr):
                return False

        if self.get_inchi() != other.get_inchi():
            return False

        if len(self.base_monomers) != len(other.base_monomers):
            return False
        for base_monomer in self.base_monomers:
            has_equal = False
            for other_base_monomer in other.base_monomers:
                if base_monomer.is_equal(other_base_monomer):
                    has_equal = True
                    break
            if not has_equal:
                return False

        return True


class MonomerSequence(list):
    """ Sequence of monomers """

    def __init__(self, monomers=None):
        """
        Args:
            monomers (:obj:iterable of :obj:`Monomer`): iterable of monomers
        """
        super(MonomerSequence, self).__init__()
        if monomers is not None:
            for monomer in monomers:
                self.append(monomer)

    def append(self, monomer):
        """ Add a monomer

        Args:
            monomer (:obj:`Monomer`): monomer

        Raises:
            :obj:`ValueError`: if the `monomer` is not an instance of `Monomer`
        """
        if not isinstance(monomer, Monomer):
            raise ValueError('`monomer` must be an instance of `Monomer`')
        super(MonomerSequence, self).append(monomer)

    def extend(self, monomers):
        """ Add a list of monomers

        Args:
            monomers (iterable of :obj:`Monomer`): iterable of monomers
        """
        for monomer in monomers:
            self.append(monomer)

    def insert(self, i, monomer):
        """ Insert a monomer at a position

        Args:
            i (:obj:`int`): position to insert monomer
            monomer (:obj:`Monomer`): monomer
        """
        if not isinstance(monomer, Monomer):
            raise ValueError('`monomer` must be an instance of `Monomer`')
        super(MonomerSequence, self).insert(i, monomer)

    def __setitem__(self, slice, monomer):
        """ Set monomer(s) at slice

        Args:
            slice (:obj:`int` or :obj:`slice`): position(s) to set monomer
            monomer (:obj:`Monomer` or :obj:`list` of :obj:`Monomer`): monomer or monomers
        """
        if isinstance(slice, int):
            if not isinstance(monomer, Monomer):
                raise ValueError('`monomer` must be a `Monomer`')
        else:
            for b in monomer:
                if not isinstance(b, Monomer):
                    raise ValueError('`monomer` must be an iterable of `Monomer`')

        super(MonomerSequence, self).__setitem__(slice, monomer)

    def get_monomer_counts(self):
        """ Get the frequency of each monomer within the sequence

        Returns:
            :obj:`dict`: dictionary that maps monomers to their counts
        """
        counts = {}
        for monomer in self:
            if monomer in counts:
                counts[monomer] += 1
            else:
                counts[monomer] = 1
        return counts

    def is_equal(self, other):
        """ Determine if two monomer sequences are semantically equal

        Args:
            other (:obj:`MonomerSequence`): other monomer sequence

        Returns:
            :obj:`bool`: True, of the monomer sequences are semantically equal
        """
        if self is other:
            return True
        if self.__class__ != other.__class__ or len(self) != len(other):
            return False
        for self_monomer, other_monomer in zip(self, other):
            if not self_monomer.is_equal(other_monomer):
                return False
        return True


class MonomerDict(attrdict.AttrDict):
    """ Dictionary for monomers """

    def __setitem__(self, chars, monomer):
        """ Set monomer with chars

        Args:
            chars (:obj:`str`): characters for monomer
            monomer (:obj:`Monomer`): monomer
        """
        if not re.match(r'^[^\[\]\{\}]+$', chars):
            raise ValueError(f'`chars` "{chars}" must be at least one character, excluding '
                             'square brackets and curly brackets')
        super(MonomerDict, self).__setitem__(chars, monomer)


class Alphabet(object):
    """ Alphabet for monomers 

    Attributes:
        id (:obj:`str`): id
        name (:obj:`str`): name
        description (:obj:`str`): description
        monomers (:obj:`dict`): monomers
    """

    def __init__(self, id=None, name=None, description=None, monomers=None):
        """
        Args:
            id (:obj:`str`, optional): id
            name (:obj:`str`, optional): name
            description (:obj:`str`, optional): description
            monomers (:obj:`dict`, optional): monomers
        """
        self.id = id
        self.name = name
        self.description = description
        self.monomers = monomers or MonomerDict()

    @property
    def monomers(self):
        """ Get the monomers

        Returns:
            :obj:`MonomerDict`: monomers
        """
        return self._monomers

    @monomers.setter
    def monomers(self, value):
        """ Set the monomers

        Args:
            value (:obj:`MonomerDict`): monomers

        Raises:
            :obj:`ValueError`: if `monomers` is not an instance of `MonomerDict`
        """
        if value is None:
            raise ValueError('`monomers` must be an instance of `MonomerDict`')
        if not isinstance(value, MonomerDict):
            value = MonomerDict(value)
        self._monomers = value

    def get_monomer_code(self, monomer):
        """ Get the code for a monomer in the alphabet

        Args:
            monomer (:obj:`Monomer`): monomer

        Returns:
            :obj:`str`: code for monomer

        Raises:
            :obj:`ValueError`: if monomer is not in alphabet
        """
        for code, alph_monomer in self.monomers.items():
            if monomer == alph_monomer:
                return code
        raise ValueError('Monomer {} is not in alphabet'.format(monomer.id))

    def protonate(self, ph):
        """ Protonate monomers

        Args:
            ph (:obj:`float`): pH
        """
        monomers = list(filter(lambda monomer: monomer.structure is not None, self.monomers.values()))

        inchis = []
        for monomer in monomers:
            inchis.append(monomer.get_inchi())

        new_inchis = Protonator.run(inchis, ph=ph)

        for monomer, new_inchi in zip(monomers, new_inchis):
            monomer.structure = new_inchi

    def is_equal(self, other):
        """ Determine two alphabets are semantically equal

        Args:
            other (:obj:`type`): other alphabet

        Returns:
            :obj:`bool`: True, if the alphabets are semantically equal
        """
        if self is other:
            return True
        if self.__class__ != other.__class__:
            return False
        for attr in ['id', 'name', 'description']:
            if getattr(self, attr) != getattr(other, attr):
                return False
        if len(self.monomers) != len(other.monomers):
            return False
        for chars, self_monomer in self.monomers.items():
            if not self_monomer.is_equal(other.monomers.get(chars, None)):
                return False
        return True

    def to_dict(self):
        """ Get dictionary representation of alphabet

        Returns:
            :obj:`dict`: dictionary representation of alphabet
        """
        dict = {}

        for attr in ['id', 'name', 'description']:
            val = getattr(self, attr)
            if val:
                dict[attr] = val

        dict['monomers'] = {}
        for chars, monomer in self.monomers.items():
            dict['monomers'][chars] = monomer.to_dict(alphabet=self)

        return dict

    def from_dict(self, dict):
        """ Create alphabet from a dictionary representation

        Args:
            dict (:obj:`dict`): dictionary representation of alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        for attr in ['id', 'name', 'description']:
            val = dict.get(attr, None)
            setattr(self, attr, val)

        self.monomers.clear()
        for chars, monomer in dict['monomers'].items():
            self.monomers[chars] = Monomer().from_dict(monomer)
        for chars, monomer in dict['monomers'].items():
            self.monomers[chars].from_dict(monomer, alphabet=self)

        return self

    def to_yaml(self, path):
        """ Save alphabet to YAML file

        Args:
            path (:obj:`str`): path to save alphabet in YAML format
        """
        yaml_writer = yaml.YAML()
        yaml_writer.default_flow_style = False
        with open(path, 'wb') as file:
            yaml_writer.dump(self.to_dict(), file)

    def from_yaml(self, path):
        """ Read alphabet from YAML file

        Args:
            path (:obj:`str`): path to YAML file which defines alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        yaml_reader = yaml.YAML()
        with open(path, 'rb') as file:
            self.from_dict(yaml_reader.load(file))
        return self


class AlphabetBuilder(abc.ABC):
    """ Builder for alphabets 

    Attributes:
        _max_monomers (:obj:`float`): maximum number of monomers to build; used to limit length of tests
    """

    def __init__(self, _max_monomers=float('inf')):
        """
        Args:
            _max_monomers (:obj:`float`, optional): maximum number of monomers to build; used to limit length of tests
        """
        self._max_monomers = _max_monomers

    def run(self, path=None):
        """ Build alphabet and, optionally, save to YAML file

        Args:
            path (:obj:`str`, optional): path to save alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        alphabet = self.build()
        if path:
            self.save(alphabet, path)
        return alphabet

    @abc.abstractmethod
    def build(self):
        """ Build alphabet 

        Returns:
            :obj:`Alphabet`: alphabet
        """
        pass  # pragma: no cover

    def save(self, alphabet, path):
        """ Save alphabet to YAML file

        Args:
            alphabet (:obj:`Alphabet`): alphabet
            path (:obj:`str`): path to save alphabet
        """
        alphabet.to_yaml(path)


class BpForm(object):
    """ Biopolymer form

    Attributes:
        monomer_seq (:obj:`MonomerSequence`): monomers of the biopolymer
        alphabet (:obj:`Alphabet`): monomer alphabet
        backbone_formula (:obj:`EmpiricalFormula`): empirical formula for backbone that connects monomers
        backbone_charge (:obj:`int`): charge for backbone that connects monomers
        bond_formula (:obj:`EmpiricalFormula`): empirical formula for bonds between monomers
        bond_charge (:obj:`int`): charge of bonds between monomers

        _parser (:obj:`lark.Lark`): parser
    """

    def __init__(self, monomer_seq=None, alphabet=None, backbone_formula=None, backbone_charge=0, bond_formula=None, bond_charge=0):
        """
        Args:
            monomer_seq (:obj:`MonomerSequence`, optional): monomers of the biopolymer
            alphabet (:obj:`Alphabet`, optional): monomer alphabet
            backbone_formula (:obj:`EmpiricalFormula`, optional): empirical formula for backbone that connects monomers
            backbone_charge (:obj:`int`, optional): charge for backbone that connects monomers
            bond_formula (:obj:`EmpiricalFormula`, optional): empirical formula for bonds between monomers
            bond_charge (:obj:`int`, optional): charge of bonds between monomers
        """
        if alphabet is None:
            alphabet = Alphabet()
        if backbone_formula is None:
            backbone_formula = EmpiricalFormula()
        if bond_formula is None:
            bond_formula = EmpiricalFormula()

        self.monomer_seq = monomer_seq or MonomerSequence()
        self.alphabet = alphabet
        self.backbone_formula = backbone_formula
        self.backbone_charge = backbone_charge
        self.bond_formula = bond_formula
        self.bond_charge = bond_charge

    @property
    def monomer_seq(self):
        """ Get the monomer sequence

        Returns:
            :obj:`MonomerSequence`: monomer sequence
        """
        return self._monomer_seq

    @monomer_seq.setter
    def monomer_seq(self, value):
        """ Set the monomer sequence

        Args:
            value (:obj:`MonomerSequence`): monomer sequence

        Raises:
            :obj:`ValueError`: if `monomer_seq` is not an instance of `MonomerSequence`
        """
        if value is None:
            raise ValueError('`monomer_seq` must be an instance of `MonomerSequence`')
        if not isinstance(value, MonomerSequence):
            value = MonomerSequence(value)
        self._monomer_seq = value

    @property
    def alphabet(self):
        """ Get the alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        return self._alphabet

    @alphabet.setter
    def alphabet(self, value):
        """ Set the monomer sequence

        Args:
            value (:obj:`Alphabet`): alphabet

        Raises:
            :obj:`ValueError`: if `alphabet` is not an instance of `Alphabet`
        """
        if not isinstance(value, Alphabet):
            raise ValueError('`alphabet` must be an instance of `Alphabet`')
        self._alphabet = value

    @property
    def backbone_formula(self):
        """ Get the formula of the backbones

        Returns:
            :obj:`EmpiricalFormula`: formula of the backbones
        """
        return self._backbone_formula

    @backbone_formula.setter
    def backbone_formula(self, value):
        """ Set the formula of the backbones

        Args:
            value (:obj:`EmpiricalFormula` or :obj:`str`): formula of the backbones
        """
        if not isinstance(value, EmpiricalFormula):
            value = EmpiricalFormula(value)
        self._backbone_mol_wt = value.get_molecular_weight()
        self._backbone_formula = value

    @property
    def backbone_charge(self):
        """ Get the backbone charge

        Returns:
            :obj:`str`: backbone charge
        """
        return self._backbone_charge

    @backbone_charge.setter
    def backbone_charge(self, value):
        """ Set the backbone charge

        Args:
            value (:obj:`str`): backbone charge

        Raises:
            :obj:`ValueError`: if the backbone charge is not an integer
        """
        if not isinstance(value, (int, float)) or value != int(value):
            raise ValueError('`backbone_charge` must be an integer')
        self._backbone_charge = int(value)

    @property
    def bond_formula(self):
        """ Get the formula of the bonds

        Returns:
            :obj:`EmpiricalFormula`: formula of the bonds
        """
        return self._bond_formula

    @bond_formula.setter
    def bond_formula(self, value):
        """ Set the formula of the bonds

        Args:
            value (:obj:`EmpiricalFormula` or :obj:`str`): formula of the bonds
        """
        if not isinstance(value, EmpiricalFormula):
            value = EmpiricalFormula(value)
        self._bond_mol_wt = value.get_molecular_weight()
        self._bond_formula = value

    @property
    def bond_charge(self):
        """ Get the bond charge

        Returns:
            :obj:`str`: bond charge
        """
        return self._bond_charge

    @bond_charge.setter
    def bond_charge(self, value):
        """ Set the bond charge

        Args:
            value (:obj:`str`): bond charge

        Raises:
            :obj:`ValueError`: if the bond charge is not an integer
        """
        if not isinstance(value, (int, float)) or value != int(value):
            raise ValueError('`bond_charge` must be an integer')
        self._bond_charge = int(value)

    def is_equal(self, other):
        """ Check if two biopolymer forms are semantically equal

        Args:
            other (:obj:`BpForm`): another biopolymer form

        Returns:
            :obj:`bool`: :obj:`True`, if the objects have the same structure
        """
        return self is other or (self.__class__ == other.__class__
                                 and self.monomer_seq.is_equal(other.monomer_seq)
                                 and self.alphabet.is_equal(other.alphabet)
                                 and self.backbone_formula == other.backbone_formula
                                 and self.backbone_charge == other.backbone_charge
                                 and self.bond_formula == other.bond_formula
                                 and self.bond_charge == other.bond_charge)

    def __getitem__(self, slice):
        """ Get monomer(s) at slice

        Args:
            slice (:obj:`int` or :obj:`slice`): position(s)

        Returns:
            :obj:`Monomer` or :obj:`Monomers`: monomer or monomers
        """
        return self.monomer_seq.__getitem__(slice)

    def __setitem__(self, slice, monomer):
        """ Set monomer(s) at slice

        Args:
            slice (:obj:`int` or :obj:`slice`): position(s)
            monomer (:obj:`Monomer` or :obj:`Monomers`): monomer or monomers
        """
        self.monomer_seq.__setitem__(slice, monomer)

    def __delitem__(self, slice):
        """ Delete monomer(s) at slice

        Args:
            slice (:obj:`int` or :obj:`slice`): position(s)
        """
        self.monomer_seq.__delitem__(slice)

    def __iter__(self):
        """ Get iterator over monomer sequence

        Returns:
            :obj:`iterator` of :obj:`Monomer`: iterator of monomers
        """
        return self.monomer_seq.__iter__()

    def __reversed__(self):
        """ Get reverse iterator over monomer sequence

        Returns:
            :obj:`iterator` of :obj:`Monomer`: iterator of monomers
        """
        return self.monomer_seq.__reversed__()

    def __contains__(self, monomer):
        """ Determine if a monomer is in the form

        Args:
            monomer (:obj:`Monomer`): monomer

        Returns:
            :obj:`bool`: true if the monomer is in the sequence
        """
        return self.monomer_seq.__contains__(monomer)

    def __len__(self):
        """ Get the length of the sequence of the form

        Returns:
            :obj:`int`: length
        """
        return len(self.monomer_seq)

    def get_monomer_counts(self):
        """ Get the frequency of each monomer within the biopolymer

        Returns:
            :obj:`dict`: dictionary that maps monomers to their counts
        """
        return self.monomer_seq.get_monomer_counts()

    def protonate(self, ph):
        """ Update to the major protonation state of each monomer at the pH

        Args:
            ph (:obj:`float`): pH
        """
        monomers = list(filter(lambda monomer: monomer.structure is not None, set(self.monomer_seq)))

        inchis = []
        for monomer in monomers:
            inchis.append(monomer.get_inchi())

        new_inchis = Protonator.run(inchis, ph=ph)

        for monomer, new_inchi in zip(monomers, new_inchis):
            monomer.structure = new_inchi

    def get_formula(self):
        """ Get the chemical formula

        Returns:
            :obj:`EmpiricalFormula`: chemical formula
        """
        formula = EmpiricalFormula()
        for monomer, count in self.get_monomer_counts().items():
            formula += monomer.get_formula() * count
        return formula + self.backbone_formula * len(self) + self.bond_formula * (len(self) - 1)

    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        mol_wt = 0.
        for monomer, count in self.get_monomer_counts().items():
            mol_wt += monomer.get_mol_wt() * count
        return mol_wt + self._backbone_mol_wt * len(self) + self._bond_mol_wt * (len(self) - 1)

    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        charge = 0
        for monomer, count in self.get_monomer_counts().items():
            charge += monomer.get_charge() * count
        return charge + self.backbone_charge * len(self) + self.bond_charge * (len(self) - 1)

    def __str__(self):
        """ Get a string representation of the biopolymer form

        Returns:
            :obj:`str`: string representation of the biopolymer form
        """
        alphabet_monomers = {monomer: chars for chars, monomer in self.alphabet.monomers.items()}
        val = ''
        for monomer in self.monomer_seq:
            chars = alphabet_monomers.get(monomer, None)
            if chars:
                if len(chars) == 1:
                    val += chars
                else:
                    val += '{' + chars + '}'
            else:
                val += monomer.__str__(alphabet=self.alphabet)
        return val

    _grammar_filename = pkg_resources.resource_filename('bpforms', 'grammar.lark')
    with open(_grammar_filename, 'r') as file:
        _parser = lark.Lark(file.read())

    def from_str(self, str):
        """ Create biopolymer form its string representation

        Args:
            str (:obj:`str`): string representation of the biopolymer

        Returns:
            :obj:`BpForm`: biopolymer form
        """
        class ParseTreeTransformer(lark.Transformer):
            def __init__(self, bp_form):
                super(ParseTreeTransformer, self).__init__()
                self.bp_form = bp_form

            @lark.v_args(inline=True)
            def seq(self, *monomer_seq):
                self.bp_form.monomer_seq.clear()
                self.bp_form.monomer_seq = monomer_seq
                return self.bp_form

            @lark.v_args(inline=True)
            def alphabet_monomer(self, chars):
                if chars[0] == '{' and chars[-1] == '}':
                    chars = chars[1:-1]
                monomer = self.bp_form.alphabet.monomers.get(chars, None)
                if monomer is None:
                    raise ValueError('"{}" not in alphabet'.format(chars))
                return monomer

            @lark.v_args(inline=True)
            def inline_monomer(self, *args):
                kwargs = {
                    'synonyms': SynonymSet(),
                    'identifiers': IdentifierSet(),
                    'base_monomers': set(),
                }
                for arg in args:
                    if isinstance(arg, tuple):
                        arg_name, arg_val = arg
                        if arg_name in ['id', 'name', 'structure', 'delta_mass', 'delta_charge', 'position', 'comments']:
                            if arg_name in kwargs:
                                raise ValueError('{} attribute cannot be repeated'.format(arg_name))
                            if arg_name == 'position':
                                kwargs['start_position'], kwargs['end_position'] = arg_val
                            else:
                                kwargs[arg_name] = arg_val
                        elif arg_name in ['synonyms', 'identifiers', 'base_monomers']:
                            kwargs[arg_name].add(arg_val)
                        else:  # pragma: no cover # the grammar ensures this will never be reached
                            raise ValueError('Invalid attribute {}'.format(arg_name))
                return Monomer(**kwargs)

            @lark.v_args(inline=True)
            def id(self, *args):
                return ('id', args[-1].value[1:-1])

            @lark.v_args(inline=True)
            def name(self, *args):
                return ('name', args[-1].value[1:-1])

            @lark.v_args(inline=True)
            def synonym(self, *args):
                return ('synonyms', args[-1].value[1:-1])

            @lark.v_args(inline=True)
            def identifier(self, *args):
                return ('identifiers', Identifier(args[-3].value[1:-1], args[-1].value[1:-1]))

            @lark.v_args(inline=True)
            def structure(self, *args):
                return ('structure', args[-1].value)

            @lark.v_args(inline=True)
            def delta_mass(self, *args):
                return ('delta_mass', float(args[-1].value))

            @lark.v_args(inline=True)
            def delta_charge(self, *args):
                return ('delta_charge', int(float(args[-1].value)))

            @lark.v_args(inline=True)
            def position(self, *args):
                start_position = None
                end_position = None
                for arg in args:
                    if arg.type == 'START_POSITION':
                        start_position = int(float(arg.value))
                    elif arg.type == 'END_POSITION':
                        end_position = int(float(arg.value))

                return ('position', (start_position, end_position))

            @lark.v_args(inline=True)
            def base_monomer(self, separator, chars):
                if chars[0] == '"' and chars[-1] == '"':
                    chars = chars[1:-1]
                monomer = self.bp_form.alphabet.monomers.get(chars, None)
                if monomer is None:
                    raise ValueError('"{}" not in alphabet'.format(chars))
                return ('base_monomers', monomer)

            @lark.v_args(inline=True)
            def comments(self, *args):
                return ('comments', args[-1].value[1:-1])

        tree = self._parser.parse(str)
        parse_tree_transformer = ParseTreeTransformer(self)
        return parse_tree_transformer.transform(tree)

    def to_fasta(self):
        """ Get FASTA representation of a monomer with bases represented by the character codes
        of their parent monomers (e.g. methyl-2-adenosine is represented by 'A')

        Returns:
            :obj:`str`: FASTA representation of a monomer
        """

        monomer_codes = {monomer: code for code, monomer in self.alphabet.monomers.items()}

        seq = ''
        for monomer in self.monomer_seq:
            root_monomers = monomer.get_root_monomers()
            root_monomers_codes = list(set(monomer_codes[root_monomer] for root_monomer in root_monomers))

            if len(root_monomers_codes) == 1:
                seq += root_monomers_codes[0]
            else:
                seq += 'N'

        return seq
