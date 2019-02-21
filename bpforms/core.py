""" Classes to represent modified forms of DNA, RNA, and proteins

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from ruamel import yaml
from wc_utils.util.chem import EmpiricalFormula, Protonator, OpenBabelUtils
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

    def protonate(self, ph, major_tautomer=False):
        """ Update to the major protonation and tautomerization state at the pH

        Args:
            ph (:obj:`float`): pH
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
        """
        if self.structure:
            self.structure = Protonator.run(self.get_inchi(), ph=ph, major_tautomer=major_tautomer)

    def get_inchi(self):
        """ Get InChI representration of structure

        Returns:
            :obj:`str`: InChI representration of structure
        """
        if self.structure:
            return OpenBabelUtils.get_inchi(self.structure)
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
        return OpenBabelUtils.get_formula(self.structure)

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

    def protonate(self, ph, major_tautomer=False):
        """ Calculate the major protonation and tautomerization of each monomer

        Args:
            ph (:obj:`float`): pH
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
        """
        monomers = list(filter(lambda monomer: monomer.structure is not None, self.monomers.values()))

        inchis = []
        for monomer in monomers:
            inchis.append(monomer.get_inchi())

        new_inchis = Protonator.run(inchis, ph=ph, major_tautomer=major_tautomer)

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

    def run(self, ph=None, major_tautomer=False, path=None):
        """ Build alphabet and, optionally, save to YAML file

        Args:
            ph (:obj:`float`, optional): pH at which to calculate the major protonation state of each monomer
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            path (:obj:`str`, optional): path to save alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        alphabet = self.build()
        if ph is not None:
            alphabet.protonate(ph, major_tautomer=major_tautomer)
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


class Atom(object):
    """ An atom in a compound or bond

    Attributes:
        element (:obj:`str`): code for the element (e.g. 'H') 
        position (:obj:`int`): IUPAC position of the atom within the compound
        charge (:obj:`int`): charge of the atom
    """

    def __init__(self, element, position=None, charge=0):
        """
        Args:
            element (:obj:`str`, optional): code for the element (e.g. 'H') 
            position (:obj:`int`, optional): IUPAC position of the atom within the compound
            charge (:obj:`int`, optional): charge of the atom
        """
        self.element = element
        self.position = position
        self.charge = charge

    @property
    def element(self):
        """ Get the element

        Returns:
            :obj:`str`: element
        """
        return self._element

    @element.setter
    def element(self, value):
        """ Set the element

        Args:
            value (:obj:`str`): element
        """
        if not isinstance(value, str):
            raise ValueError('`element` must be a string')
        self._element = value

    @property
    def position(self):
        """ Get the position

        Returns:
            :obj:`int`: position
        """
        return self._position

    @position.setter
    def position(self, value):
        """ Set the position

        Args:
            value (:obj:`int`): position
        """
        if value is not None:
            if (not isinstance(value, (int, float)) or value != int(value) or value < 1):
                raise ValueError('`position` must be a positive integer or None')
            value = int(value)
        self._position = value

    @property
    def charge(self):
        """ Get the charge

        Returns:
            :obj:`str`: charge
        """
        return self._charge

    @charge.setter
    def charge(self, value):
        """ Set the charge

        Args:
            value (:obj:`str`): charge
        """
        if not isinstance(value, (float, int)) or int(value) != value:
            raise ValueError('`charge` must be an integer')
        value = int(value)
        self._charge = value

    def is_equal(self, other):
        """ Determine if two atoms are semantically equal

        Args:
            other (:obj:`Atom`): other atom

        Returns:
            :obj:`bool`: obj:`True` if the atoms are semantically equal
        """
        if self is other:
            return True
        if self.__class__ != other.__class__:
            return False
        return self is other or (self.__class__ == other.__class__
                                 and self.element == other.element
                                 and self.charge == other.charge
                                 and self.position == other.position)


class AtomList(list):
    """ List of atoms """

    def __init__(self, atoms=None):
        """
        Args:
            atoms (:obj:iterable of :obj:`Atom`): iterable of atoms
        """
        super(AtomList, self).__init__()
        if atoms is not None:
            for atom in atoms:
                self.append(atom)

    def append(self, atom):
        """ Add a atom

        Args:
            atom (:obj:`Atom`): atom

        Raises:
            :obj:`ValueError`: if the `atom` is not an instance of `Atom`
        """
        if not isinstance(atom, Atom):
            raise ValueError('`atom` must be an instance of `Atom`')
        super(AtomList, self).append(atom)

    def extend(self, atoms):
        """ Add a list of atoms

        Args:
            atoms (iterable of :obj:`Atom`): iterable of atoms
        """
        for atom in atoms:
            self.append(atom)

    def insert(self, i, atom):
        """ Insert an atom at a position

        Args:
            i (:obj:`int`): position to insert atom
            atom (:obj:`Atom`): atom
        """
        if not isinstance(atom, Atom):
            raise ValueError('`atom` must be an instance of `Atom`')
        super(AtomList, self).insert(i, atom)

    def __setitem__(self, slice, atom):
        """ Set atom(s) at slice

        Args:
            slice (:obj:`int` or :obj:`slice`): position(s) to set atom
            atom (:obj:`Atom` or :obj:`AtomList`): atom or atoms
        """
        if isinstance(slice, int):
            if not isinstance(atom, Atom):
                raise ValueError('`atom` must be a `Atom`')
        else:
            for b in atom:
                if not isinstance(b, Atom):
                    raise ValueError('`atom` must be an iterable of `Atom`')

        super(AtomList, self).__setitem__(slice, atom)

    def is_equal(self, other):
        """ Determine if two lists of atoms are semantically equal

        Args:
            other (:obj:`AtomList`): other list of atoms

        Returns:
            :obj:`bool`: True, of the lists of atoms are semantically equal
        """
        if self is other:
            return True
        if self.__class__ != other.__class__ or len(self) != len(other):
            return False
        for self_atom, other_atom in zip(self, other):
            if not self_atom.is_equal(other_atom):
                return False
        return True


class Backbone(object):
    """ Backbone of a monomer

    Attributes:
        structure (:obj:`openbabel.OBMol`): chemical structure
        backbone_bond_atoms (:obj:`AtomList`): atoms from backbone that bonds to monomer
        monomer_bond_atoms (:obj:`AtomList`): atoms from monomer that bonds to backbone
        backbone_displaced_atoms (:obj:`AtomList`): atoms from backbone displaced by bond to monomer
        monomer_displaced_atoms (:obj:`AtomList`): atoms from monomer displaced by bond to backbone
    """

    def __init__(self, structure=None, backbone_bond_atoms=None, monomer_bond_atoms=None,
                 backbone_displaced_atoms=None, monomer_displaced_atoms=None):
        """
        Args:
            structure (:obj:`str` or :obj:`openbabel.OBMol`, optional): chemical structure as InChI-encoded string or OpenBabel molecule
            backbone_bond_atoms (:obj:`AtomList`, optional): atoms from backbone that bonds to monomer
            monomer_bond_atoms (:obj:`AtomList`, optional): atoms from monomer that bonds to backbone
            backbone_displaced_atoms (:obj:`AtomList`, optional): atoms from backbone displaced by bond to monomer
            monomer_displaced_atoms (:obj:`AtomList`, optional): atoms from monomer displaced by bond to backbone
        """
        self.structure = structure
        self.backbone_bond_atoms = backbone_bond_atoms or AtomList()
        self.monomer_bond_atoms = monomer_bond_atoms or AtomList()
        self.backbone_displaced_atoms = backbone_displaced_atoms or AtomList()
        self.monomer_displaced_atoms = monomer_displaced_atoms or AtomList()

    @property
    def structure(self):
        """ Get the structure

        Returns:
            :obj:`openbabel.OBMol`: structure
        """
        return self._structure

    @structure.setter
    def structure(self, value):
        """ Set the structure

        Args:
            value (:obj:`str` or :obj:`openbabel.OBMol`): structure as InChI-encoded string or OpenBabel molecule
        """
        if value is not None and not isinstance(value, openbabel.OBMol):
            ob_mol = openbabel.OBMol()
            conversion = openbabel.OBConversion()
            assert conversion.SetInFormat('inchi'), 'Unable to set format to InChI'
            if not conversion.ReadString(ob_mol, value):
                raise ValueError('`structure` must be an OpenBabel molecule, InChI-encoded structure, or None')
            value = ob_mol
        self._structure = value

    @property
    def monomer_bond_atoms(self):
        """ Get the monomer bond atoms

        Returns:
            :obj:`AtomList`: monomer bond atoms
        """
        return self._monomer_bond_atoms

    @monomer_bond_atoms.setter
    def monomer_bond_atoms(self, value):
        """ Set the monomer bond atoms

        Args:
            value (:obj:`AtomList`): monomer bond atoms

        Raises:
            :obj:`ValueError`: if `monomer_bond_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`monomer_bond_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._monomer_bond_atoms = value

    @property
    def backbone_bond_atoms(self):
        """ Get the backbone bond atoms

        Returns:
            :obj:`AtomList`: backbone bond atoms
        """
        return self._backbone_bond_atoms

    @backbone_bond_atoms.setter
    def backbone_bond_atoms(self, value):
        """ Set the backbone bond atoms

        Args:
            value (:obj:`AtomList`): backbone bond atoms

        Raises:
            :obj:`ValueError`: if `backbone_bond_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`backbone_bond_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._backbone_bond_atoms = value

    @property
    def monomer_displaced_atoms(self):
        """ Get the monomer displaced atoms

        Returns:
            :obj:`AtomList`: monomer displaced atoms
        """
        return self._monomer_displaced_atoms

    @monomer_displaced_atoms.setter
    def monomer_displaced_atoms(self, value):
        """ Set the monomer displaced atoms

        Args:
            value (:obj:`AtomList`): monomer displaced atoms

        Raises:
            :obj:`ValueError`: if `monomer_displaced_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`monomer_displaced_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._monomer_displaced_atoms = value

    @property
    def backbone_displaced_atoms(self):
        """ Get the backbone displaced atoms

        Returns:
            :obj:`AtomList`: backbone displaced atoms
        """
        return self._backbone_displaced_atoms

    @backbone_displaced_atoms.setter
    def backbone_displaced_atoms(self, value):
        """ Set the backbone displaced atoms

        Args:
            value (:obj:`AtomList`): backbone displaced atoms

        Raises:
            :obj:`ValueError`: if `backbone_displaced_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`backbone_displaced_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._backbone_displaced_atoms = value

    def get_inchi(self):
        """ Get InChI representration of structure

        Returns:
            :obj:`str`: InChI representration of structure
        """
        if self.structure:
            return OpenBabelUtils.get_inchi(self.structure)
        return None

    def get_formula(self):
        """ Get the formula

        Returns:
            :obj:`EmpiricalFormula`: formula
        """
        if self.structure:
            formula = OpenBabelUtils.get_formula(self.structure)
        else:
            formula = EmpiricalFormula()
        for atom in self.backbone_displaced_atoms:
            formula[atom.element] -= 1
        for atom in self.monomer_displaced_atoms:
            formula[atom.element] -= 1
        return formula

    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        return self.get_formula().get_molecular_weight()

    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        if self.structure:
            charge = self.structure.GetTotalCharge()
        else:
            charge = 0
        for atom in self.backbone_displaced_atoms:
            charge -= atom.charge
        for atom in self.monomer_displaced_atoms:
            charge -= atom.charge
        return charge

    def is_equal(self, other):
        """ Determine if two backbones are semantically equal

        Args:
            other (:obj:`Backbone`): other backbone

        Returns:
            :obj:`bool`: :obj:`True` if the backbones are semantically equal
        """
        if self is other:
            return True
        if self.__class__ != other.__class__:
            return False
        if self.get_inchi() != other.get_inchi():
            return False
        if not self.monomer_bond_atoms.is_equal(other.monomer_bond_atoms)\
                or not self.backbone_bond_atoms.is_equal(other.backbone_bond_atoms)\
                or not self.monomer_displaced_atoms.is_equal(other.monomer_displaced_atoms)\
                or not self.backbone_displaced_atoms.is_equal(other.backbone_displaced_atoms):
            return False
        return True


class Bond(object):
    """ Bond between monomers

    Attributes:
        left_participant (:obj:`type`): type of left participant (monomer or backbone)
        right_participant (:obj:`type`): type of right participant (monomer or backbone)
        left_bond_atoms (:obj:`AtomList`): atoms from left monomer that bonds with right monomer
        right_bond_atoms (:obj:`AtomList`): atoms from right monomer that bonds with left monomer
        left_displaced_atoms (:obj:`AtomList`): atoms from left monomer displaced by bond
        right_displaced_atoms (:obj:`AtomList`): atoms from right monomer displaced by bond
    """

    def __init__(self, left_participant=None, right_participant=None, left_bond_atoms=None, right_bond_atoms=None,
                 left_displaced_atoms=None, right_displaced_atoms=None):
        """
        Args:
            left_participant (:obj:`type`, optional): type of left participant (monomer or backbone)
            right_participant (:obj:`type`, optional): type of right participant (monomer or backbone)
            left_bond_atoms (:obj:`AtomList`, optional): atoms from left monomer that bonds with right monomer
            right_bond_atoms (:obj:`AtomList`, optional): atoms from right monomer that bonds with left monomer
            left_displaced_atoms (:obj:`AtomList`, optional): atoms from left monomer displaced by bond
            right_displaced_atoms (:obj:`AtomList`, optional): atoms from right monomer displaced by bond
        """
        self.left_participant = left_participant
        self.right_participant = right_participant
        self.left_bond_atoms = left_bond_atoms or AtomList()
        self.right_bond_atoms = right_bond_atoms or AtomList()
        self.left_displaced_atoms = left_displaced_atoms or AtomList()
        self.right_displaced_atoms = right_displaced_atoms or AtomList()

    @property
    def left_participant(self):
        """ Get type of the left participant

        Returns:
            :obj:`type`: type of the left participant
        """
        return self._left_participant

    @left_participant.setter
    def left_participant(self, value):
        """ Set the type of the left participant

        Args:
            value (:obj:`type`): type of the left participant

        Raises:
            :obj:`ValueError`: if `left_participant` is not :obj:`None`, :obj:`Monomer`, or :obj:`Backbone`
        """
        if value not in [None, Monomer, Backbone]:
            raise ValueError('`left_participant` must be `None`, `Monomer`, or `Backbone`')
        self._left_participant = value

    @property
    def right_participant(self):
        """ Get type of the right participant

        Returns:
            :obj:`type`: type of the right participant
        """
        return self._right_participant

    @right_participant.setter
    def right_participant(self, value):
        """ Set the type of the right participant

        Args:
            value (:obj:`type`): type of the right participant

        Raises:
            :obj:`ValueError`: if `right_participant` is not :obj:`None`, :obj:`Monomer`, or :obj:`Backbone`
        """
        if value not in [None, Monomer, Backbone]:
            raise ValueError('`right_participant` must be `None`, `Monomer`, or `Backbone`')
        self._right_participant = value

    @property
    def left_bond_atoms(self):
        """ Get the left bond atoms

        Returns:
            :obj:`AtomList`: left bond atoms
        """
        return self._left_bond_atoms

    @left_bond_atoms.setter
    def left_bond_atoms(self, value):
        """ Set the left bond atoms

        Args:
            value (:obj:`AtomList`): left bond atoms

        Raises:
            :obj:`ValueError`: if `left_bond_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`left_bond_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._left_bond_atoms = value

    @property
    def right_bond_atoms(self):
        """ Get the right bond atoms

        Returns:
            :obj:`AtomList`: right bond atoms
        """
        return self._right_bond_atoms

    @right_bond_atoms.setter
    def right_bond_atoms(self, value):
        """ Set the right bond atoms

        Args:
            value (:obj:`AtomList`): right bond atoms

        Raises:
            :obj:`ValueError`: if `right_bond_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`right_bond_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._right_bond_atoms = value

    @property
    def left_displaced_atoms(self):
        """ Get the left displaced atoms

        Returns:
            :obj:`AtomList`: left displaced atoms
        """
        return self._left_displaced_atoms

    @left_displaced_atoms.setter
    def left_displaced_atoms(self, value):
        """ Set the left displaced atoms

        Args:
            value (:obj:`AtomList`): left displaced atoms

        Raises:
            :obj:`ValueError`: if `left_displaced_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`left_displaced_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._left_displaced_atoms = value

    @property
    def right_displaced_atoms(self):
        """ Get the right displaced atoms

        Returns:
            :obj:`AtomList`: right displaced atoms
        """
        return self._right_displaced_atoms

    @right_displaced_atoms.setter
    def right_displaced_atoms(self, value):
        """ Set the right displaced atoms

        Args:
            value (:obj:`AtomList`): right displaced atoms

        Raises:
            :obj:`ValueError`: if `right_displaced_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`right_displaced_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._right_displaced_atoms = value

    def get_formula(self):
        """ Get the formula

        Returns:
            :obj:`EmpiricalFormula`: formula
        """
        formula = EmpiricalFormula()
        for atom in self.left_displaced_atoms:
            formula[atom.element] -= 1
        for atom in self.right_displaced_atoms:
            formula[atom.element] -= 1
        return formula

    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        return self.get_formula().get_molecular_weight()

    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        charge = 0
        for atom in self.left_displaced_atoms:
            charge -= atom.charge
        for atom in self.right_displaced_atoms:
            charge -= atom.charge
        return charge

    def is_equal(self, other):
        """ Determine if two bonds are semantically equal

        Args:
            other (:obj:`Bond`): other bond

        Returns:
            :obj:`bool`: :obj:`True` if the bond are semantically equal
        """
        if self is other:
            return True
        if self.__class__ != other.__class__:
            return False
        if not self.left_bond_atoms.is_equal(other.left_bond_atoms) \
                or not self.right_bond_atoms.is_equal(other.right_bond_atoms) \
                or not self.left_displaced_atoms.is_equal(other.left_displaced_atoms) \
                or not self.right_displaced_atoms.is_equal(other.right_displaced_atoms):
            return False
        return True


class BpForm(object):
    """ Biopolymer form

    Attributes:
        monomer_seq (:obj:`MonomerSequence`): monomers of the biopolymer
        alphabet (:obj:`Alphabet`): monomer alphabet
        backbone (:obj:`Backbone`): backbone that connects monomers
        bond (:obj:`Bond`): bonds between (backbones of) monomers
        circular (:obj:`bool`): if :obj:`True`, indicates that the biopolymer is circular

        _parser (:obj:`lark.Lark`): parser
    """

    DEFAULT_FASTA_CODE = '?'

    def __init__(self, monomer_seq=None, alphabet=None, backbone=None, bond=None, circular=False):
        """
        Args:
            monomer_seq (:obj:`MonomerSequence`, optional): monomers of the biopolymer
            alphabet (:obj:`Alphabet`, optional): monomer alphabet
            backbone (:obj:`Backbone`, optional): backbone that connects monomers
            bond (:obj:`Bond`, optional): bonds between (backbones of) monomers
            circular (:obj:`bool`, optional): if :obj:`True`, indicates that the biopolymer is circular
        """
        if alphabet is None:
            alphabet = Alphabet()
        if backbone is None:
            backbone = Backbone()
        if bond is None:
            bond = Bond()

        self.monomer_seq = monomer_seq or MonomerSequence()
        self.alphabet = alphabet
        self.backbone = backbone
        self.bond = bond
        self.circular = circular

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
    def backbone(self):
        """ Get the backbones

        Returns:
            :obj:`Backbone`: backbones
        """
        return self._backbone

    @backbone.setter
    def backbone(self, value):
        """ Set the backbones

        Args:
            value (:obj:`Backbone`): backbones
        """
        if not isinstance(value, Backbone):
            raise ValueError('`backbone` must be an instance of `Backbone`')
        self._backbone = value

    @property
    def bond(self):
        """ Get the bonds

        Returns:
            :obj:`Bond`: bonds
        """
        return self._bond

    @bond.setter
    def bond(self, value):
        """ Set the bonds

        Args:
            value (:obj:`Bond`): bonds
        """
        if not isinstance(value, Bond):
            raise ValueError('`bond` must be an instance of `Bond`')
        self._bond = value

    @property
    def circular(self):
        """ Get the circularity

        Returns:
            :obj:`bool`: circularity
        """
        return self._circular

    @circular.setter
    def circular(self, value):
        """ Set the circularity

        Args:
            value (:obj:`bool`): circularity
        """
        if not isinstance(value, bool):
            raise ValueError('`circular` must be an instance of `bool`')
        self._circular = value

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
                                 and self.backbone.is_equal(other.backbone)
                                 and self.bond.is_equal(other.bond)
                                 and self.circular == other.circular)

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

    def protonate(self, ph, major_tautomer=True):
        """ Update to the major protonation and tautomerization state of each monomer at the pHf

        Args:
            ph (:obj:`float`): pH
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
        """
        monomers = list(filter(lambda monomer: monomer.structure is not None, set(self.monomer_seq)))

        inchis = []
        for monomer in monomers:
            inchis.append(monomer.get_inchi())

        new_inchis = Protonator.run(inchis, ph=ph, major_tautomer=major_tautomer)

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
        return formula + self.backbone.get_formula() * len(self) + self.bond.get_formula() * (len(self) - (1 - self.circular))

    def get_mol_wt(self):
        """ Get the molecular weight

        Returns:
            :obj:`float`: molecular weight
        """
        mol_wt = 0.
        for monomer, count in self.get_monomer_counts().items():
            mol_wt += monomer.get_mol_wt() * count
        return mol_wt + self.backbone.get_mol_wt() * len(self) + self.bond.get_mol_wt() * (len(self) - (1 - self.circular))

    def get_charge(self):
        """ Get the charge

        Returns:
            :obj:`int`: charge
        """
        charge = 0
        for monomer, count in self.get_monomer_counts().items():
            charge += monomer.get_charge() * count
        return charge + self.backbone.get_charge() * len(self) + self.bond.get_charge() * (len(self) - (1 - self.circular))

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
                return ('structure', re.sub(r'[ \t\f\r\n]', '', args[-1].value))

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
                seq += self.DEFAULT_FASTA_CODE

        return seq
