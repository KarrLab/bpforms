""" Classes to represent modified forms of DNA, RNA, and proteins

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from ruamel import yaml
from wc_utils.util.chem import EmpiricalFormula, get_major_micro_species, draw_molecule, OpenBabelUtils
import abc
import attrdict
import itertools
import lark
import openbabel
import os
import pkg_resources
import re
import urllib.parse
import warnings
import wc_utils.cache

# setup cache
cache_dir = os.path.expanduser('~/.cache/bpforms')
if not os.path.isdir(cache_dir):
    os.makedirs(cache_dir)
cache = wc_utils.cache.Cache(directory=cache_dir)


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
    """ A monomeric form in a biopolymer

    Attributes:
        id (:obj:`str`): id
        name (:obj:`str`): name
        synonyms (:obj:`set` of :obj:`str`): synonyms
        identifiers (:obj:`set` of :obj:`Identifier`, optional): identifiers in namespaces for external databases
        structure (:obj:`openbabel.OBMol`): chemical structure
        delta_mass (:obj:`float`): additional mass (Dalton) relative to structure
        delta_charge (:obj:`int`): additional charge relative to structure
        start_position (:obj:`tuple`): uncertainty in the location of the monomeric form
        end_position (:obj:`tuple`): uncertainty in the location of the monomeric form
        monomers_position (:obj:`set` of :obj:`Monomer`): originating monomers within :obj:`start_position` to 
            :obj:`end_position` where the monomeric form may be located
        base_monomers (:obj:`set` of :obj:`Monomer`): monomers which this monomeric form is derived from
        backbone_bond_atoms (:obj:`AtomList`): atoms from monomeric form that bond to backbone
        backbone_displaced_atoms (:obj:`AtomList`): atoms from monomeric form displaced by bond to backbone
        r_bond_atoms (:obj:`AtomList`): atoms that bond with right/suceeding/following/forward monomeric form
        l_bond_atoms (:obj:`AtomList`): atoms that bond with left/preceding/previous/backward monomeric form
        r_displaced_atoms (:obj:`AtomList`): atoms displaced by bond with right/suceeding/following/forward monomeric form
        l_displaced_atoms (:obj:`AtomList`): atoms displaced by bond with left/preceding/previous/backward monomeric form
        comments (:obj:`str`): comments
    """

    def __init__(self, id=None, name=None, synonyms=None, identifiers=None, structure=None,
                 delta_mass=None, delta_charge=None,
                 start_position=None, end_position=None, monomers_position=None,
                 base_monomers=None,
                 backbone_bond_atoms=None, backbone_displaced_atoms=None,
                 r_bond_atoms=None, l_bond_atoms=None,
                 r_displaced_atoms=None, l_displaced_atoms=None,
                 comments=None):
        """
        Attributes:
            id (:obj:`str`, optional): id
            name (:obj:`str`, optional): name
            synonyms (:obj:`set` of :obj:`str`, optional): synonyms
            identifiers (:obj:`set` of :obj:`Identifier`, optional): identifiers in namespaces for external databases
            structure (:obj:`openbabel.OBMol` or :obj:`str`, optional): chemical structure
            delta_mass (:obj:`float`, optional): additional mass (Dalton) relative to structure
            delta_charge (:obj:`float`, optional): additional charge relative to structure
            start_position (:obj:`int`, optional): uncertainty in the location of the monomeric form
            end_position (:obj:`int`, optional): uncertainty in the location of the monomeric form
            monomers_position (:obj:`set` of :obj:`Monomer`, optional): originating monomers within :obj:`start_position` to 
                :obj:`end_position` where the monomeric form may be located
            base_monomers (:obj:`set` of :obj:`Monomer`, optional): monomers which this monomeric form is derived from
            backbone_bond_atoms (:obj:`AtomList`, optional): atoms from monomeric form that bond to backbone
            backbone_displaced_atoms (:obj:`AtomList`, optional): atoms from monomeric form displaced by bond to backbone
            r_bond_atoms (:obj:`AtomList`, optional): atoms that bond with right/suceeding/following/forward monomeric form
            l_bond_atoms (:obj:`AtomList`, optional): atoms that bond with left/preceding/previous/backward monomeric form
            r_displaced_atoms (:obj:`AtomList`, optional): atoms displaced by bond with right/suceeding/following/forward monomeric form
            l_displaced_atoms (:obj:`AtomList`, optional): atoms displaced by bond with left/preceding/previous/backward monomeric form
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
        self.monomers_position = monomers_position or set()
        self.base_monomers = base_monomers or set()
        self.backbone_bond_atoms = backbone_bond_atoms or AtomList()
        self.backbone_displaced_atoms = backbone_displaced_atoms or AtomList()
        self.r_bond_atoms = r_bond_atoms or AtomList()
        self.l_bond_atoms = l_bond_atoms or AtomList()
        self.r_displaced_atoms = r_displaced_atoms or AtomList()
        self.l_displaced_atoms = l_displaced_atoms or AtomList()
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
            value (:obj:`openbabel.OBMol` or :obj:`str`): Open Babel molecule, canonical SMILES-encoded structure, or None

        Raises:
            :obj:`ValueError`: if value is not an Open Babel molecule, canonical SMILES-encoded structure, or None
        """
        if value and not isinstance(value, openbabel.OBMol):
            ob_mol = openbabel.OBMol()
            conversion = openbabel.OBConversion()
            assert conversion.SetInFormat('smi'), 'Unable to set format to SMILES'
            if not conversion.ReadString(ob_mol, value):
                raise ValueError('`structure` must be an Open Babel molecule, canonical SMILES-encoded structure, or None')
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
    def monomers_position(self):
        """ Get the originating monomers within :obj:`start_position` to 
        :obj:`end_position` where the monomeric form may be located

        Returns:
            :obj:`set` of :obj:`Monomer`: originating monomers within :obj:`start_position` to 
                :obj:`end_position` where the monomeric form may be located
        """
        return self._monomers_position

    @monomers_position.setter
    def monomers_position(self, value):
        """ Set the originating monomers within :obj:`start_position` to 
        :obj:`end_position` where the monomeric form may be located

        Args:
            value (:obj:`set` of :obj:`Monomer`): originating monomers within :obj:`start_position` to 
                :obj:`end_position` where the monomeric form may be located

        Raises:
            :obj:`ValueError`: if value is not an instance of :obj:`set`
        """
        if isinstance(value, list):
            value = set(value)
        if not isinstance(value, set):
            raise ValueError('`monomers_position` must be an instance of `set`')
        self._monomers_position = value

    @property
    def base_monomers(self):
        """ Get base monomeric forms

        Returns:
            :obj:`set` of :obj:`Monomer`: base monomeric forms
        """
        return self._base_monomers

    @base_monomers.setter
    def base_monomers(self, value):
        """ Set base monomeric forms

        Args:
            value (:obj:`set` of :obj:`Monomer`): base monomeric forms

        Raises:
            :obj:`ValueError`: if value is not an instance of :obj:`set`
        """
        if isinstance(value, list):
            value = set(value)
        if not isinstance(value, set):
            raise ValueError('`base_monomers` must be an instance of `set`')
        self._base_monomers = value

    @property
    def backbone_bond_atoms(self):
        """ Get the atoms from the monomeric form that bond to backbone

        Returns:
            :obj:`AtomList`: atoms from the monomeric form that bond to backbone
        """
        return self._backbone_bond_atoms

    @backbone_bond_atoms.setter
    def backbone_bond_atoms(self, value):
        """ Set the atoms from the monomeric form that bond to backbone

        Args:
            value (:obj:`AtomList`): atoms from the monomeric form that bond to backbone

        Raises:
            :obj:`ValueError`: if `backbone_bond_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`backbone_bond_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._backbone_bond_atoms = value

    @property
    def backbone_displaced_atoms(self):
        """ Get the atoms from the monomeric form displaced by the bond to the backbone

        Returns:
            :obj:`AtomList`: atoms from the monomeric form displaced by the bond to the backbone
        """
        return self._backbone_displaced_atoms

    @backbone_displaced_atoms.setter
    def backbone_displaced_atoms(self, value):
        """ Set the atoms from the monomeric form displaced by the bond to the backbone

        Args:
            value (:obj:`AtomList`): atoms from the monomeric form displaced by the bond to the backbone

        Raises:
            :obj:`ValueError`: if `backbone_displaced_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`backbone_displaced_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._backbone_displaced_atoms = value

    @property
    def r_bond_atoms(self):
        """ Get the left bond atoms

        Returns:
            :obj:`AtomList`: left bond atoms
        """
        return self._r_bond_atoms

    @r_bond_atoms.setter
    def r_bond_atoms(self, value):
        """ Set the left bond atoms

        Args:
            value (:obj:`AtomList`): left bond atoms

        Raises:
            :obj:`ValueError`: if `r_bond_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`r_bond_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._r_bond_atoms = value

    @property
    def l_bond_atoms(self):
        """ Get the right bond atoms

        Returns:
            :obj:`AtomList`: right bond atoms
        """
        return self._l_bond_atoms

    @l_bond_atoms.setter
    def l_bond_atoms(self, value):
        """ Set the right bond atoms

        Args:
            value (:obj:`AtomList`): right bond atoms

        Raises:
            :obj:`ValueError`: if `l_bond_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`l_bond_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._l_bond_atoms = value

    @property
    def r_displaced_atoms(self):
        """ Get the left displaced atoms

        Returns:
            :obj:`AtomList`: left displaced atoms
        """
        return self._r_displaced_atoms

    @r_displaced_atoms.setter
    def r_displaced_atoms(self, value):
        """ Set the left displaced atoms

        Args:
            value (:obj:`AtomList`): left displaced atoms

        Raises:
            :obj:`ValueError`: if `r_displaced_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`r_displaced_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._r_displaced_atoms = value

    @property
    def l_displaced_atoms(self):
        """ Get the right displaced atoms

        Returns:
            :obj:`AtomList`: right displaced atoms
        """
        return self._l_displaced_atoms

    @l_displaced_atoms.setter
    def l_displaced_atoms(self, value):
        """ Set the right displaced atoms

        Args:
            value (:obj:`AtomList`): right displaced atoms

        Raises:
            :obj:`ValueError`: if `l_displaced_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`l_displaced_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._l_displaced_atoms = value

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
        """ Get root monomeric forms

        Returns:
            :obj:`set` of :obj:`Monomer`: root monomeric forms
        """
        if not self.base_monomers:
            return set([self])

        roots = set()
        for base_monomer in self.base_monomers:
            roots.update(base_monomer.get_root_monomers())

        return roots

    def get_major_micro_species(self, ph, major_tautomer=False, dearomatize=False):
        """ Update to the major protonation and tautomerization state at the pH

        Args:
            ph (:obj:`float`): pH
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule
        """
        if self.structure:
            structure = get_major_micro_species(self.export('smi', options=('c',)),
                                                'smiles', 'smiles', ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)
            conv = openbabel.OBConversion()
            assert conv.SetInFormat('smi')
            assert conv.ReadString(self.structure, structure)

    def export(self, format, options=()):
        """ Export structure to format

        Args:
            format (:obj:`str`): format
            options (:obj:`list` of :obj:`str`, optional): export options

        Returns:
            :obj:`str`: format representation of structure
        """
        if self.structure:
            return OpenBabelUtils.export(self.structure, format, options=options)
        else:
            return None

    IMAGE_URL_PATTERN = ('https://cactus.nci.nih.gov/chemical/structure/{}/image'
                         '?format=gif'
                         '&bgcolor=transparent'
                         '&antialiasing=0')

    def get_image_url(self):
        """ Get URL for image of structure

        Returns:
            :obj:`str`: URL for image of structure
        """
        smiles = self.export('smi', options=('c',))
        if smiles:
            return self.IMAGE_URL_PATTERN.format(urllib.parse.quote(smiles))
        return None

    def get_image(self, bond_label='', displaced_label='', bond_opacity=255, displaced_opacity=63,
                  backbone_bond_color=0xff0000, left_bond_color=0x00ff00, right_bond_color=0x0000ff,
                  include_all_hydrogens=True, show_atom_nums=False,
                  width=200, height=200, image_format='svg', include_xml_header=True):
        """ Get image

        Args:
            bond_label (:obj:`str`, optional): label for atoms involved in bonds
            displaced_label (:obj:`str`, optional): labels for atoms displaced by bond formation
            bond_opacity (:obj:`int`, optional): opacity of atoms involved in bonds
            displaced_opacity (:obj:`int`, optional): opacity of atoms dislaced by bond formation
            backbone_bond_color (:obj:`int`, optional): color to paint atoms involved in bond with backbone
            left_bond_color (:obj:`int`, optional): color to paint atoms involved in bond with monomeric form to left
            right_bond_color (:obj:`int`, optional): color to paint atoms involved in bond with monomeric form to right
            include_all_hydrogens (:obj:`bool`, optional): if :obj:`True`, show all hydrogens
            show_atom_nums (:obj:`bool`, optional): if :obj:`True`, show the numbers of the atoms
            width (:obj:`int`, optional): width in pixels
            height (:obj:`int`, optional): height in pixels
            image_format (:obj:`str`, optional): format of generated image {emf, eps, jpeg, msbmp, pdf, png, or svg}
            include_xml_header (:obj:`bool`, optional): if :obj:`True`, include XML header at the beginning of the SVG

        Returns:
            :obj:`object`: image
        """
        if not self.structure:
            return None

        atom_md_types = [
            (self.backbone_bond_atoms, bond_label, self._blend_color_opacity(backbone_bond_color, bond_opacity)),
            (self.backbone_displaced_atoms, displaced_label, self._blend_color_opacity(backbone_bond_color, displaced_opacity)),
            (self.r_bond_atoms, bond_label, self._blend_color_opacity(left_bond_color, bond_opacity)),
            (self.r_displaced_atoms, displaced_label, self._blend_color_opacity(left_bond_color, displaced_opacity)),
            (self.l_bond_atoms, bond_label, self._blend_color_opacity(right_bond_color, bond_opacity)),
            (self.l_displaced_atoms, displaced_label, self._blend_color_opacity(right_bond_color, displaced_opacity)),
        ]
        atom_labels = []
        atom_sets = {}

        mol = openbabel.OBMol()
        mol += self.structure
        if include_all_hydrogens:
            mol.AddHydrogens()

        for atom_mds, label, color in atom_md_types:
            bonding_hydrogens = []
            for atom_md in atom_mds:
                atom = mol.GetAtom(atom_md.position)
                if atom_md.element == 'H' and atom.GetAtomicNum() != 1:
                    atom = get_hydrogen_atom(atom, bonding_hydrogens, 0)

                if atom:
                    atom_labels.append({'position': atom.GetIdx(), 'element': atom_md.element, 'label': label, 'color': color})
                    if color not in atom_sets:
                        atom_sets[color] = {'positions': [], 'elements': [], 'color': color}
                    atom_sets[color]['positions'].append(atom.GetIdx())
                    atom_sets[color]['elements'].append(atom_md.element)

        cml = OpenBabelUtils.export(mol, 'cml')

        return draw_molecule(cml, 'cml', image_format=image_format,
                             atom_labels=atom_labels, atom_sets=atom_sets.values(),
                             show_atom_nums=show_atom_nums,
                             width=width, height=height, include_xml_header=include_xml_header)

    @staticmethod
    def _blend_color_opacity(color, opacity):
        """ Blend color with white to simulate opacity

        Args:
            color (:obj:`int`): color (0-0xffffff)
            opacity (:obj:`int`): opacity (0-0xff)

        Returns:
            :obj:`int`: blended color
        """
        r = color >> 16
        g = (color - r * 2**16) >> 8
        b = color - r * 2**16 - g * 2**8
        w = 0xff

        r_opaque = round((r * opacity + w * (255 - opacity)) / 255)
        g_opaque = round((g * opacity + w * (255 - opacity)) / 255)
        b_opaque = round((b * opacity + w * (255 - opacity)) / 255)

        return (r_opaque << 16) + (g_opaque << 8) + b_opaque

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
        """ Get a dictionary representation of the monomeric form

        Args:
            alphabet (:obj:`Alphabet`, optional): alphabet

        Returns:
            :obj:`dict`: dictionary representation of the monomeric form
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
            dict['structure'] = self.export('smi', options=('c',))

        if self.monomers_position and alphabet:
            dict['monomers_position'] = []
            for monomer in self.monomers_position:
                monomer_code = alphabet.get_monomer_code(monomer)
                dict['monomers_position'].append(monomer_code)

        if self.base_monomers and alphabet:
            dict['base_monomers'] = []
            for monomer in self.base_monomers:
                monomer_code = alphabet.get_monomer_code(monomer)
                dict['base_monomers'].append(monomer_code)

        attr_names = (
            'backbone_bond_atoms', 'backbone_displaced_atoms',
            'r_bond_atoms', 'l_bond_atoms',
            'r_displaced_atoms', 'l_displaced_atoms')
        for attr_name in attr_names:
            atoms = getattr(self, attr_name)
            if atoms:
                dict[attr_name] = []
                for atom in atoms:
                    dict[attr_name].append(atom.to_dict())

        return dict

    def from_dict(self, dict, alphabet=None):
        """ Get a dictionary representation of the monomeric form

        Args:
            dict (:obj:`dict`): dictionary representation of the monomeric form
            alphabet (:obj:`Alphabet`, optional): alphabet

        Returns:
            :obj:`Monomer`: monomeric form
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
        self.monomers_position.clear()
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

        monomers_position_ids = dict.get('monomers_position', [])
        if monomers_position_ids and alphabet:
            self.monomers_position = set([alphabet.monomers.get(monomer_id) for monomer_id in monomers_position_ids])

        base_monomer_ids = dict.get('base_monomers', [])
        if base_monomer_ids and alphabet:
            self.base_monomers = set([alphabet.monomers.get(monomer_id) for monomer_id in base_monomer_ids])

        attr_names = (
            'backbone_bond_atoms', 'backbone_displaced_atoms',
            'r_bond_atoms', 'l_bond_atoms',
            'r_displaced_atoms', 'l_displaced_atoms')
        for attr_name in attr_names:
            atoms = getattr(self, attr_name)
            atoms.from_list(dict.get(attr_name, []))

        return self

    def __str__(self, alphabet=None):
        """ Get a string representation of the monomeric form

        Args:
            alphabet (:obj:`Alphabet`, optional): alphabet

        Returns:
            :obj:`str`: string representation of the monomeric form
        """
        els = []

        if self.id:
            els.append('id: "' + self.id + '"')

        if self.name:
            els.append('name: "' + self.name.replace('"', '\\"') + '"')

        for synonym in self.synonyms:
            els.append('synonym: "' + synonym.replace('"', '\\"') + '"')

        for identifier in self.identifiers:
            els.append('identifier: "' + identifier.id + '" @ "' + identifier.ns + '"')

        if self.structure:
            els.append('structure: "' + self.export('smi', options=('c',)) + '"')

        atom_types = [
            'backbone_bond_atoms', 'backbone_displaced_atoms',
            'r_bond_atoms', 'r_displaced_atoms',
            'l_bond_atoms', 'l_displaced_atoms',
        ]
        for atom_type in atom_types:
            for atom in getattr(self, atom_type):
                if atom.charge > 0:
                    charge = '+' + str(atom.charge)
                elif atom.charge == 0:
                    charge = ''
                else:
                    charge = str(atom.charge)
                els.append('{}: {}{}{}'.format(
                    atom_type[:-1].replace('_', '-'), atom.element, atom.position, charge))

        if self.delta_mass is not None:
            els.append('delta-mass: ' + str(self.delta_mass))

        if self.delta_charge is not None:
            els.append('delta-charge: ' + str(self.delta_charge))

        if self.start_position is not None or self.end_position is not None:
            el = 'position: {}-{}'.format(self.start_position or '', self.end_position or '')

            if self.monomers_position and alphabet:
                codes = sorted(alphabet.get_monomer_code(monomer) for monomer in self.monomers_position)
                el += ' [{}]'.format(' | '.join(codes))

            els.append(el)

        if alphabet:
            for monomer in self.base_monomers:
                els.append('base-monomer: "{}"'.format(alphabet.get_monomer_code(monomer)))

        if self.comments:
            els.append('comments: "' + self.comments.replace('"', '\\"') + '"')

        return '[' + ' | '.join(els) + ']'

    def get_canonical_code(self, monomer_codes, default_code='?'):
        """ Get IUPAC/IUBMB representation of a monomeric form using the character code
        of its parent monomer (e.g. 'methyl-2-adenosine' is represented by 'A')

        Args:
            monomer_codes (:obj:`dict`): dictionary that maps monomeric forms to codes
            default_code (:obj:`str`): default code

        Returns:
            :obj:`str`: IUPAC/IUBMB representation of monomeric form
        """
        roots = self.get_root_monomers()
        root_codes = list(set(monomer_codes.get(root, default_code) for root in roots))

        if len(root_codes) == 1:
            return root_codes[0]
        else:
            return default_code

    def is_equal(self, other):
        """ Check if two monomeric forms are semantically equal

        Args:
            other (:obj:`Monomer`): another monomeric form

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

        if self.export('inchi') != other.export('inchi'):
            return False

        for attr in ['monomers_position', 'base_monomers']:
            if len(getattr(self, attr)) != len(getattr(other, attr)):
                return False
            for base_monomer in getattr(self, attr):
                has_equal = False
                for other_base_monomer in getattr(other, attr):
                    if base_monomer.is_equal(other_base_monomer):
                        has_equal = True
                        break
                if not has_equal:
                    return False

        attr_names = (
            'backbone_bond_atoms', 'backbone_displaced_atoms',
            'r_bond_atoms', 'l_bond_atoms',
            'r_displaced_atoms', 'l_displaced_atoms')
        for attr_name in attr_names:
            self_atoms = getattr(self, attr_name)
            other_atoms = getattr(other, attr_name)
            if not self_atoms.is_equal(other_atoms):
                return False

        return True


class MonomerSequence(list):
    """ Sequence of monomeric forms """

    def __init__(self, monomers=None):
        """
        Args:
            monomers (:obj:iterable of :obj:`Monomer`): iterable of monomeric forms
        """
        super(MonomerSequence, self).__init__()
        if monomers is not None:
            for monomer in monomers:
                self.append(monomer)

    def append(self, monomer):
        """ Add a monomeric form

        Args:
            monomer (:obj:`Monomer`): monomeric form

        Raises:
            :obj:`ValueError`: if the `monomer` is not an instance of `Monomer`
        """
        if not isinstance(monomer, Monomer):
            raise ValueError('`monomer` must be an instance of `Monomer`')
        super(MonomerSequence, self).append(monomer)

    def extend(self, monomers):
        """ Add a list of monomeric forms

        Args:
            monomers (iterable of :obj:`Monomer`): iterable of monomeric forms
        """
        for monomer in monomers:
            self.append(monomer)

    def insert(self, i, monomer):
        """ Insert a monomeric form at a position

        Args:
            i (:obj:`int`): position to insert monomeric form
            monomer (:obj:`Monomer`): monomeric form
        """
        if not isinstance(monomer, Monomer):
            raise ValueError('`monomer` must be an instance of `Monomer`')
        super(MonomerSequence, self).insert(i, monomer)

    def __setitem__(self, slice, monomer):
        """ Set monomeric form(s) at slice

        Args:
            slice (:obj:`int` or :obj:`slice`): position(s) to set monomeric form
            monomer (:obj:`Monomer` or :obj:`list` of :obj:`Monomer`): monomeric form(s)
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
        """ Get the frequency of each monomeric form within the sequence

        Returns:
            :obj:`dict`: dictionary that maps monomeric forms to their counts
        """
        counts = {}
        for monomer in self:
            if monomer in counts:
                counts[monomer] += 1
            else:
                counts[monomer] = 1
        return counts

    def is_equal(self, other):
        """ Determine if two sequences of monomeric forms are semantically equal

        Args:
            other (:obj:`MonomerSequence`): other sequence

        Returns:
            :obj:`bool`: True, of the sequences are semantically equal
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
    """ Dictionary for monomeric forms """

    def __setitem__(self, code, monomer):
        """ Set monomeric form with code

        Args:
            code (:obj:`str`): characters for monomeric form
            monomer (:obj:`Monomer`): monomeric form
        """
        pattern = (r'^([^\[\]\{\}":\| \t\f\r\n]|'
                   r'[^\[\]\{\}":\| \t\f\r\n][^\[\]\{\}":\|\t\f\r\n]*[^\[\]\{\}":\| \t\f\r\n])$')
        if not re.match(pattern, code):
            raise ValueError(f'`code` "{code}" must be at least one character, excluding '
                             'square brackets, curly brackets, double quotes, colons, and pipes '
                             'and must not begin or end with a white space')
        super(MonomerDict, self).__setitem__(code, monomer)


class Alphabet(object):
    """ Alphabet for monomeric forms

    Attributes:
        id (:obj:`str`): id
        name (:obj:`str`): name
        description (:obj:`str`): description
        monomers (:obj:`dict`): monomeric forms
    """

    def __init__(self, id=None, name=None, description=None, monomers=None):
        """
        Args:
            id (:obj:`str`, optional): id
            name (:obj:`str`, optional): name
            description (:obj:`str`, optional): description
            monomers (:obj:`dict`, optional): monomeric forms
        """
        self.id = id
        self.name = name
        self.description = description
        self.monomers = monomers or MonomerDict()

    @property
    def monomers(self):
        """ Get the monomeric forms

        Returns:
            :obj:`MonomerDict`: monomeric forms
        """
        return self._monomers

    @monomers.setter
    def monomers(self, value):
        """ Set the monomeric forms

        Args:
            value (:obj:`MonomerDict`): monomeric forms

        Raises:
            :obj:`ValueError`: if `monomers` is not an instance of `MonomerDict`
        """
        if value is None:
            raise ValueError('`monomers` must be an instance of `MonomerDict`')
        if not isinstance(value, MonomerDict):
            value = MonomerDict(value)
        self._monomers = value

    def get_monomer_code(self, monomer):
        """ Get the code for a monomeric form in the alphabet

        Args:
            monomer (:obj:`Monomer`): monomeric form

        Returns:
            :obj:`str`: code for monomeric form

        Raises:
            :obj:`ValueError`: if monomeric form is not in alphabet
        """
        for code, alph_monomer in self.monomers.items():
            if monomer == alph_monomer:
                return code
        raise ValueError('Monomer {} is not in alphabet'.format(monomer.id))

    def get_major_micro_species(self, ph, major_tautomer=False, dearomatize=False):
        """ Calculate the major protonation and tautomerization of each monomeric form

        Args:
            ph (:obj:`float`): pH
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule
        """
        monomers = list(filter(lambda monomer: monomer.structure is not None, self.monomers.values()))

        structures = []
        for monomer in monomers:
            structure = monomer.export('smi', options=('c',))
            structures.append(structure)

        new_structures = get_major_micro_species(structures, 'smiles', 'smiles',
                                                 ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)

        for monomer, new_structure in zip(monomers, new_structures):
            conv = openbabel.OBConversion()
            assert conv.SetInFormat('smi')
            assert conv.ReadString(monomer.structure, new_structure)

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
        for code, self_monomer in self.monomers.items():
            if not self_monomer.is_equal(other.monomers.get(code, None)):
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
        for code, monomer in self.monomers.items():
            dict['monomers'][code] = monomer.to_dict(alphabet=self)

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
        for code, monomer in dict['monomers'].items():
            self.monomers[code] = Monomer().from_dict(monomer)
        for code, monomer in dict['monomers'].items():
            self.monomers[code].from_dict(monomer, alphabet=self)

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
        self.from_dict(parse_yaml(path))
        return self


@cache.memoize(typed=False, expire=30 * 24 * 60 * 60, filename_args=[0])
def parse_yaml(path):
    """ Read a YAML file

    Args:
        path (:obj:`str`): path to YAML file which defines alphabet

    Returns:
        :obj:`object`: content of file
    """
    yaml_reader = yaml.YAML()
    with open(path, 'rb') as file:
        return yaml_reader.load(file)


class AlphabetBuilder(abc.ABC):
    """ Builder for alphabets

    Attributes:
        _max_monomers (:obj:`float`): maximum number of monomeric forms to build; used to limit length of tests
    """

    def __init__(self, _max_monomers=float('inf')):
        """
        Args:
            _max_monomers (:obj:`float`, optional): maximum number of monomeric forms to build; used to limit length of tests
        """
        self._max_monomers = _max_monomers

    def run(self, ph=None, major_tautomer=False, dearomatize=False, path=None):
        """ Build alphabet and, optionally, save to YAML file

        Args:
            ph (:obj:`float`, optional): pH at which to calculate the major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule
            path (:obj:`str`, optional): path to save alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        alphabet = self.build(ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)
        if path:
            self.save(alphabet, path)
        return alphabet

    @abc.abstractmethod
    def build(self, ph=None, major_tautomer=False, dearomatize=False):
        """ Build alphabet

        Args:
            ph (:obj:`float`, optional): pH at which to calculate the major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule

        Returns:
            :obj:`Alphabet`: alphabet
        """
        pass  # pragma: no cover

    def get_major_micro_species(self, alphabet, ph=None, major_tautomer=False, dearomatize=False):
        """ Get major microspecies for monomeric forms in alphabet

        Args:
            ph (:obj:`float`, optional): pH at which to calculate the major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule
        """
        if ph is not None:
            alphabet.get_major_micro_species(ph, major_tautomer=major_tautomer, dearomatize=dearomatize)

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
        molecule (:obj:`type`): type of parent molecule
        element (:obj:`str`): code for the element (e.g. 'H')
        position (:obj:`int`): position of the atom within the molecule, which should use canonical SMILES atom numbers
        charge (:obj:`int`): charge of the atom
        monomer (:obj:`int`): index of parent monomeric form within sequence
    """

    def __init__(self, molecule, element, position=None, charge=0, monomer=None):
        """
        Args:
            molecule (:obj:`type`): type of parent molecule
            element (:obj:`str`, optional): code for the element (e.g. 'H')
            position (:obj:`int`, optional): position of the atom within the molecule, which should use
                canonical SMILES atom numbers
            charge (:obj:`int`, optional): charge of the atom
            monomer (:obj:`int`, optional): index of parent monomeric form within sequence
        """
        self.molecule = molecule
        self.element = element
        self.position = position
        self.charge = charge
        self.monomer = monomer

    @property
    def molecule(self):
        """ Get type of parent molecule

        Returns:
            :obj:`type`: type of parent molecule
        """
        return self._molecule

    @molecule.setter
    def molecule(self, value):
        """ Set the type of parent molecule

        Args:
            value (:obj:`type`): type of parent molecule

        Raises:
            :obj:`ValueError`: if `molecule` is not :obj:`None`, :obj:`Monomer`, or :obj:`Backbone`
        """
        if value not in [None, Monomer, Backbone]:
            raise ValueError('`molecule` must be `None`, `Monomer`, or `Backbone`')
        self._molecule = value

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

    @property
    def monomer(self):
        """ Get the index of the parent monomer within the sequence

        Returns:
            :obj:`int`: index of the parent monomer within the sequence
        """
        return self._monomer

    @monomer.setter
    def monomer(self, value):
        """ Set the index of the parent monomer within the sequence

        Args:
            value (:obj:`int`): index of the parent monomer within the sequence
        """
        if value is not None:
            if (not isinstance(value, (int, float)) or value != int(value) or value < 1):
                raise ValueError('`monomer` must be a positive integer or None')
            value = int(value)
        self._monomer = value

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
                                 and self.molecule == other.molecule
                                 and self.element == other.element
                                 and self.position == other.position
                                 and self.charge == other.charge
                                 and self.monomer == other.monomer)

    def to_dict(self):
        """ Get dictionary representation

        Returns:
            :obj:`dict`: dictionary representation
        """
        dict = {}
        if self.molecule:
            dict['molecule'] = self.molecule.__name__
        dict['element'] = self.element
        if self.position is not None:
            dict['position'] = self.position
        if self.charge:
            dict['charge'] = self.charge
        if self.monomer is not None:
            dict['monomer'] = self.monomer
        return dict

    def from_dict(self, dict):
        """ Load from dictionary representation

        Args:
            dict (:obj:`dict`): dictionary representation

        Returns:
            :obj:`Atom`: atom
        """
        molecule = dict.get('molecule', None)
        if molecule == 'Monomer':
            self.molecule = Monomer
        elif molecule == 'Backbone':
            self.molecule = Backbone
        else:
            self.molecule = None
        self.element = dict['element']
        self.position = dict.get('position', None)
        self.charge = dict.get('charge', 0)
        self.monomer = dict.get('monomer', None)
        return self


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

    def to_list(self):
        """ Get list representation

        Returns:
            :obj:`list`: list representation
        """
        list = []
        for atom in self:
            list.append(atom.to_dict())
        return list

    def from_list(self, list):
        """ Load from list representation

        Args:
            list (:obj:`list`): list representation

        Returns:
            :obj:`AtomList`: atom list
        """
        self.clear()
        for atom in list:
            self.append(Atom(None, '').from_dict(atom))
        return self


class Backbone(object):
    """ Backbone of a monomeric form

    Attributes:
        structure (:obj:`openbabel.OBMol`): chemical structure
        monomer_bond_atoms (:obj:`AtomList`): atoms from backbone that bond to monomeric form
        monomer_displaced_atoms (:obj:`AtomList`): atoms from backbone displaced by bond to monomeric form
    """

    def __init__(self, structure=None, monomer_bond_atoms=None, monomer_displaced_atoms=None):
        """
        Args:
            structure (:obj:`str` or :obj:`openbabel.OBMol`, optional): chemical structure as canonical 
                SMILES-encoded string or Open Babel molecule
            monomer_bond_atoms (:obj:`AtomList`, optional): atoms from backbone that bond to monomeric form
            monomer_displaced_atoms (:obj:`AtomList`, optional): atoms from backbone displaced by bond to monomeric form
        """
        self.structure = structure
        self.monomer_bond_atoms = monomer_bond_atoms or AtomList()
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
            value (:obj:`str` or :obj:`openbabel.OBMol`): structure as canonical SMILES-encoded string or Open Babel molecule
        """
        if value is not None and not isinstance(value, openbabel.OBMol):
            ob_mol = openbabel.OBMol()
            conversion = openbabel.OBConversion()
            assert conversion.SetInFormat('smi'), 'Unable to set format to SMILES'
            if not conversion.ReadString(ob_mol, value):
                raise ValueError('`structure` must be an Open Babel molecule, canonical SMILES-encoded structure, or None')
            value = ob_mol
        self._structure = value

    @property
    def monomer_bond_atoms(self):
        """ Get the backbone bond atoms

        Returns:
            :obj:`AtomList`: backbone bond atoms
        """
        return self._monomer_bond_atoms

    @monomer_bond_atoms.setter
    def monomer_bond_atoms(self, value):
        """ Set the backbone bond atoms

        Args:
            value (:obj:`AtomList`): backbone bond atoms

        Raises:
            :obj:`ValueError`: if `monomer_bond_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`monomer_bond_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._monomer_bond_atoms = value

    @property
    def monomer_displaced_atoms(self):
        """ Get the backbone displaced atoms

        Returns:
            :obj:`AtomList`: backbone displaced atoms
        """
        return self._monomer_displaced_atoms

    @monomer_displaced_atoms.setter
    def monomer_displaced_atoms(self, value):
        """ Set the backbone displaced atoms

        Args:
            value (:obj:`AtomList`): backbone displaced atoms

        Raises:
            :obj:`ValueError`: if `monomer_displaced_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`monomer_displaced_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._monomer_displaced_atoms = value

    def export(self, format, options=()):
        """ Export structure to format

        Args:
            format (:obj:`str`): format
            options (:obj:`list` of :obj:`str`, optional): export options

        Returns:
            :obj:`str`: format representation of structure
        """
        if self.structure:
            return OpenBabelUtils.export(self.structure, format, options=options)
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
        for atom in self.monomer_bond_atoms:
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
        if self.export('inchi') != other.export('inchi'):
            return False
        if not self.monomer_bond_atoms.is_equal(other.monomer_bond_atoms)\
                or not self.monomer_displaced_atoms.is_equal(other.monomer_displaced_atoms):
            return False
        return True


class Bond(object):
    """ Bond between monomeric forms

    Attributes:
        l_bond_atoms (:obj:`AtomList`): atoms from left monomeric form that bond with right monomeric form
        r_bond_atoms (:obj:`AtomList`): atoms from right monomeric form that bond with left monomeric form
        l_displaced_atoms (:obj:`AtomList`): atoms from left monomeric form displaced by bond
        r_displaced_atoms (:obj:`AtomList`): atoms from right monomeric form displaced by bond
    """

    def __init__(self, l_bond_atoms=None, r_bond_atoms=None,
                 l_displaced_atoms=None, r_displaced_atoms=None):
        """
        Args:
            l_bond_atoms (:obj:`AtomList`, optional): atoms from left monomeric form that bond with right monomeric form
            r_bond_atoms (:obj:`AtomList`, optional): atoms from right monomeric form that bond with left monomeric form
            l_displaced_atoms (:obj:`AtomList`, optional): atoms from left monomeric form displaced by bond
            r_displaced_atoms (:obj:`AtomList`, optional): atoms from right monomeric form displaced by bond
        """
        self.l_bond_atoms = l_bond_atoms or AtomList()
        self.r_bond_atoms = r_bond_atoms or AtomList()
        self.l_displaced_atoms = l_displaced_atoms or AtomList()
        self.r_displaced_atoms = r_displaced_atoms or AtomList()

    @property
    def l_bond_atoms(self):
        """ Get the left bond atoms

        Returns:
            :obj:`AtomList`: left bond atoms
        """
        return self._l_bond_atoms

    @l_bond_atoms.setter
    def l_bond_atoms(self, value):
        """ Set the left bond atoms

        Args:
            value (:obj:`AtomList`): left bond atoms

        Raises:
            :obj:`ValueError`: if `l_bond_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`l_bond_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._l_bond_atoms = value

    @property
    def r_bond_atoms(self):
        """ Get the right bond atoms

        Returns:
            :obj:`AtomList`: right bond atoms
        """
        return self._r_bond_atoms

    @r_bond_atoms.setter
    def r_bond_atoms(self, value):
        """ Set the right bond atoms

        Args:
            value (:obj:`AtomList`): right bond atoms

        Raises:
            :obj:`ValueError`: if `r_bond_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`r_bond_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._r_bond_atoms = value

    @property
    def l_displaced_atoms(self):
        """ Get the left displaced atoms

        Returns:
            :obj:`AtomList`: left displaced atoms
        """
        return self._l_displaced_atoms

    @l_displaced_atoms.setter
    def l_displaced_atoms(self, value):
        """ Set the left displaced atoms

        Args:
            value (:obj:`AtomList`): left displaced atoms

        Raises:
            :obj:`ValueError`: if `l_displaced_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`l_displaced_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._l_displaced_atoms = value

    @property
    def r_displaced_atoms(self):
        """ Get the right displaced atoms

        Returns:
            :obj:`AtomList`: right displaced atoms
        """
        return self._r_displaced_atoms

    @r_displaced_atoms.setter
    def r_displaced_atoms(self, value):
        """ Set the right displaced atoms

        Args:
            value (:obj:`AtomList`): right displaced atoms

        Raises:
            :obj:`ValueError`: if `r_displaced_atoms` is not an instance of `AtomList`
        """
        if value is None:
            raise ValueError('`r_displaced_atoms` must be an instance of `AtomList`')
        if not isinstance(value, AtomList):
            value = AtomList(value)
        self._r_displaced_atoms = value

    def get_formula(self, none_position=True):
        """ Get the formula

        Args:
            none_position (:obj:`bool`, optional): include atoms whose position is :obj:`None`

        Returns:
            :obj:`EmpiricalFormula`: formula
        """
        formula = EmpiricalFormula()
        for atom in self.l_displaced_atoms:
            if atom.position is not None or none_position:
                formula[atom.element] -= 1
        for atom in self.r_displaced_atoms:
            if atom.position is not None or none_position:
                formula[atom.element] -= 1
        return formula

    def get_mol_wt(self, none_position=True):
        """ Get the molecular weight

        Args:
            none_position (:obj:`bool`, optional): include atoms whose position is :obj:`None`

        Returns:
            :obj:`float`: molecular weight
        """
        return self.get_formula(none_position=none_position).get_molecular_weight()

    def get_charge(self, none_position=True):
        """ Get the charge

        Args:
            none_position (:obj:`bool`, optional): include atoms whose position is :obj:`None`

        Returns:
            :obj:`int`: charge
        """
        charge = 0
        for atom in self.l_bond_atoms:
            if atom.position is not None or none_position:
                charge -= atom.charge
        for atom in self.l_displaced_atoms:
            if atom.position is not None or none_position:
                charge -= atom.charge
        for atom in self.r_bond_atoms:
            if atom.position is not None or none_position:
                charge -= atom.charge
        for atom in self.r_displaced_atoms:
            if atom.position is not None or none_position:
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
        if not self.l_bond_atoms.is_equal(other.l_bond_atoms) \
                or not self.r_bond_atoms.is_equal(other.r_bond_atoms) \
                or not self.l_displaced_atoms.is_equal(other.l_displaced_atoms) \
                or not self.r_displaced_atoms.is_equal(other.r_displaced_atoms):
            return False
        return True

    def __str__(self):
        """ Generate string representation of bond

        Returns:
            :obj:`str`: string representation of bond
        """
        atoms = []
        for atom_type in ['l_bond_atoms', 'r_bond_atoms',
                          'l_displaced_atoms',  'r_displaced_atoms']:
            for atom in getattr(self, atom_type):
                if atom.charge > 0:
                    charge = '+' + str(atom.charge)
                elif atom.charge < 0:
                    charge = str(atom.charge)
                else:
                    charge = ''

                atoms.append('{}: {}{}{}{}'.format(atom_type[0:-1].replace('_', '-'),
                                                   atom.monomer or '',
                                                   atom.element,
                                                   atom.position or '',
                                                   charge))
        return "[{}]".format(' | '.join(atoms))


class BondSet(set):
    """ Set of bonds """

    def add(self, bond):
        """ Add a bond

        Args:
            bond (:obj:`Bond`): bond

        Raises:
            :obj:`ValueError`: if the `bond` is not an instance of `Bond`
        """
        if not isinstance(bond, Bond):
            raise ValueError('`bond` must be an instance of `Bond`')
        super(BondSet, self).add(bond)

    def update(self, bonds):
        """ Add a set of bonds

        Args:
            bonds (iterable of :obj:`Bond`): bonds
        """
        for bond in bonds:
            self.add(bond)

    def symmetric_difference_update(self, other):
        """ Remove common elements with other and add elements from other not in self

        Args:
            other (:obj:`BondSet`): other set of bonds
        """
        for o in other:
            if o in self:
                self.remove(o)
            else:
                self.add(o)

    def is_equal(self, other):
        """ Check if two sets of bonds are semantically equal

        Args:
            other (:obj:`BondSet`): other set of bonds

        Returns:
            :obj:`bool`: :obj:`True`, if the bond sets are semantically equal
        """
        if self is other:
            return True

        if self.__class__ != other.__class__:
            return False

        o_bonds = set(other)
        for s_bond in self:
            match = False
            for o_bond in set(o_bonds):
                if s_bond.is_equal(o_bond):
                    match = True
                    o_bonds.remove(o_bond)
                    break
            if not match:
                return False
        return True


class BpForm(object):
    """ Biopolymer form

    Attributes:
        seq (:obj:`MonomerSequence`): sequence of monomeric forms of the biopolymer
        alphabet (:obj:`Alphabet`): alphabet of monomeric forms
        backbone (:obj:`Backbone`): backbone that connects monomeric forms
        bond (:obj:`Bond`): bonds between (backbones of) monomeric forms
        circular (:obj:`bool`): if :obj:`True`, indicates that the biopolymer is circular
        crosslinks (:obj:`BondSet`): crosslinking intrachain bonds
        features (:obj:`BpFormFeatureSet`): set of features

        _parser (:obj:`lark.Lark`): parser
    """

    DEFAULT_FASTA_CODE = '?'

    def __init__(self, seq=None, alphabet=None, backbone=None, bond=None, circular=False, crosslinks=None):
        """
        Args:
            seq (:obj:`MonomerSequence`, optional): sequence of monomeric forms of the biopolymer
            alphabet (:obj:`Alphabet`, optional): alphabet of monomeric forms
            backbone (:obj:`Backbone`, optional): backbone that connects monomeric forms
            bond (:obj:`Bond`, optional): bonds between (backbones of) monomeric forms
            circular (:obj:`bool`, optional): if :obj:`True`, indicates that the biopolymer is circular
            crosslinks (:obj:`BondSet`, optional): crosslinking intrachain bonds
        """
        if alphabet is None:
            alphabet = Alphabet()
        if backbone is None:
            backbone = Backbone()
        if bond is None:
            bond = Bond()
        if crosslinks is None:
            crosslinks = BondSet()

        self.seq = seq or MonomerSequence()
        self.alphabet = alphabet
        self.backbone = backbone
        self.bond = bond
        self.circular = circular
        self.crosslinks = crosslinks
        self.features = BpFormFeatureSet(self)

    @property
    def seq(self):
        """ Get the sequence of monomeric forms

        Returns:
            :obj:`MonomerSequence`: sequence of monomeric forms
        """
        return self._monomer_seq

    @seq.setter
    def seq(self, value):
        """ Set the sequence of monomeric forms

        Args:
            value (:obj:`MonomerSequence`): sequence of monomeric forms

        Raises:
            :obj:`ValueError`: if `seq` is not an instance of `MonomerSequence`
        """
        if value is None:
            raise ValueError('`seq` must be an instance of `MonomerSequence`')
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
        """ Set the sequence of monomeric forms

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

    @property
    def crosslinks(self):
        """ Get the crosslinking intrachain bonds

        Returns:
            :obj:`BondSet`: crosslinking intrachain bonds
        """
        return self._crosslinks

    @crosslinks.setter
    def crosslinks(self, value):
        """ Set the crosslinking intrachain bonds

        Args:
            value (:obj:`list` of :obj:`BondSet`): crosslinking intrachain bonds
        """
        if not isinstance(value, BondSet):
            raise ValueError('`crosslinks` must be an instance of `BondSet`')
        self._crosslinks = value

    @property
    def features(self):
        """ Get the features

        Returns:
            :obj:`BpFormFeatureSet`: features
        """
        return self._features

    @features.setter
    def features(self, value):
        """ Set the features

        Args:
            value (:obj:`BpFormFeatureSet`): features
        """
        if not isinstance(value, BpFormFeatureSet):
            raise ValueError('`features` must be an instance of `BpFormFeatureSet`')

        if hasattr(self, '_features'):
            raise ValueError('`features` cannot be set outside constructor')

        self._features = value

    def is_equal(self, other):
        """ Check if two biopolymer forms are semantically equal

        Args:
            other (:obj:`BpForm`): another biopolymer form

        Returns:
            :obj:`bool`: :obj:`True`, if the objects have the same structure
        """
        return self is other or (self.__class__ == other.__class__
                                 and self.seq.is_equal(other.seq)
                                 and self.alphabet.is_equal(other.alphabet)
                                 and self.backbone.is_equal(other.backbone)
                                 and self.bond.is_equal(other.bond)
                                 and self.circular == other.circular
                                 and self.crosslinks.is_equal(other.crosslinks))

    def __getitem__(self, slice):
        """ Get monomeric form(s) at slice

        Args:
            slice (:obj:`int` or :obj:`slice`): position(s)

        Returns:
            :obj:`Monomer` or :obj:`Monomers`: monomeric form(s)
        """
        return self.seq.__getitem__(slice)

    def __setitem__(self, slice, monomer):
        """ Set monomeric form(s) at slice

        Args:
            slice (:obj:`int` or :obj:`slice`): position(s)
            monomer (:obj:`Monomer` or :obj:`Monomers`): monomeric forms(s)
        """
        self.seq.__setitem__(slice, monomer)

    def __delitem__(self, slice):
        """ Delete monomeric form(s) at slice

        Args:
            slice (:obj:`int` or :obj:`slice`): position(s)
        """
        self.seq.__delitem__(slice)

    def __iter__(self):
        """ Get iterator over sequence of monomeric forms

        Returns:
            :obj:`iterator` of :obj:`Monomer`: iterator of monomeric forms
        """
        return self.seq.__iter__()

    def __reversed__(self):
        """ Get reverse iterator over sequence of monomeric forms

        Returns:
            :obj:`iterator` of :obj:`Monomer`: iterator of monomeric forms
        """
        return self.seq.__reversed__()

    def __contains__(self, monomer):
        """ Determine if a monomeric form is in the biopolymer form

        Args:
            monomer (:obj:`Monomer`): monomeric form

        Returns:
            :obj:`bool`: true if the monomeric form is in the sequence
        """
        return self.seq.__contains__(monomer)

    def __len__(self):
        """ Get the length of the sequence of the form

        Returns:
            :obj:`int`: length
        """
        return len(self.seq)

    def validate(self):
        """ Check that the biopolymer form is valid and return any errors

        * Check that monomeric forms :math:`1 \ldots L-1` can bond to the right (their right bonding attributes are set)
        * Check that monomeric forms :math:`2 \ldots L` can bond to the left (their left bonding attributes are set)
        * No atom is involved in multiple bonds

        Returns:
            :obj:`list` of :obj:`str`: list of errors, if any
        """
        errors = []

        bonding_backbone_hydrogens = []
        bonding_monomer_hydrogens = []

        el_table = openbabel.OBElementTable()

        def check_atom(molecule, atom_type, i_monomer, i_atom, structure, atom_md, bonding_hydrogens, el_table, errors):
            if structure is None:
                errors.append('Structure cannot be None')
                return

            atom_obj = structure.GetAtom(atom_md.position)
            if atom_md.element == 'H':
                atom_obj = get_hydrogen_atom(atom_obj, bonding_hydrogens, i_monomer)
            if atom_obj:
                if el_table.GetSymbol(atom_obj.GetAtomicNum()) != atom_md.element:
                    errors.append("{} atom '{}[{}]' must be {}".format(molecule, atom_type, i_atom, atom_md.element))

        # check that bond atoms are defined
        atom_types = ['monomer_bond_atoms', 'monomer_displaced_atoms']
        for atom_type in atom_types:
            for i_atom, atom in enumerate(getattr(self.backbone, atom_type)):
                if atom.molecule != Backbone or not atom.element or not atom.position:
                    errors.append("Backbone atom '{}[{}]' must have a defined element and position".format(atom_type, i_atom))
                else:
                    check_atom('Backbone', atom_type, None, i_atom, self.backbone.structure, atom,
                               bonding_backbone_hydrogens, el_table, errors)

        atom_types = ['l_bond_atoms', 'l_displaced_atoms',
                      'r_bond_atoms', 'r_displaced_atoms']
        for atom_type in atom_types:
            if len(set(atom.molecule for atom in getattr(self.bond, atom_type))) > 1:
                errors.append("'{}' must have the same molecule type".format(atom_type))

            for i_atom, atom in enumerate(getattr(self.bond, atom_type)):
                if not atom.element or (atom.molecule == Backbone and not atom.position):
                    errors.append("Bond atom '{}[{}]' must have a defined element and position".format(atom_type, i_atom))

                elif atom.molecule == Backbone:
                    check_atom('Bond', atom_type, None, i_atom, self.backbone.structure, atom,
                               bonding_backbone_hydrogens, el_table, errors)

        for i_monomer, monomer in enumerate(self.seq):
            if not monomer.structure:
                errors.append('Monomer {} must have a defined structure'.format(i_monomer + 1))
                continue

            atom_types = ['backbone_bond_atoms', 'backbone_displaced_atoms']
            for atom_type in atom_types:
                for i_atom, atom in enumerate(getattr(monomer, atom_type)):
                    if atom.molecule != Monomer or not atom.element or not atom.position:
                        errors.append("'{}[{}]' of monomer {} must have a defined element and position".format(
                            atom_type, i_atom, i_monomer + 1))
                    else:
                        check_atom('Monomer {}'.format(i_monomer + 1), atom_type, i_monomer, i_atom,
                                   monomer.structure, atom, bonding_monomer_hydrogens, el_table, errors)

            atom_types = ['r_bond_atoms', 'r_displaced_atoms',
                          'l_bond_atoms', 'l_displaced_atoms']
            for atom_type in atom_types:
                for i_atom, atom in enumerate(getattr(monomer, atom_type)):
                    if atom.molecule != Monomer or not atom.element or not atom.position:
                        errors.append("'{}[{}]' of monomer {} must have a defined element and position".format(
                            atom_type, i_atom, i_monomer + 1))
                    else:
                        check_atom('Monomer {}'.format(i_monomer + 1), atom_type, i_monomer, i_atom,
                                   monomer.structure, atom, bonding_monomer_hydrogens, el_table, errors)

        # crosslinks
        for i_crosslink, crosslink in enumerate(self.crosslinks):
            atom_types = ['r_bond_atoms', 'r_displaced_atoms',
                          'l_bond_atoms', 'l_displaced_atoms']
            for atom_type in atom_types:
                for i_atom, atom in enumerate(getattr(crosslink, atom_type)):
                    if atom.molecule != Monomer or not atom.monomer or not atom.element or not atom.position:
                        errors.append("'{}[{}]' of crosslink {} must have a defined monomer, element, and position".format(
                            atom_type, i_atom, i_crosslink + 1))
                    else:
                        check_atom('Crosslink {} - monomer {}'.format(i_crosslink, atom.monomer), atom_type, atom.monomer - 1, i_atom,
                                   self.seq[atom.monomer - 1].structure, atom, bonding_monomer_hydrogens,
                                   el_table, errors)

        # check that monomers 1 .. L-1 can bond to right
        for i_monomer, monomer in enumerate(self.seq[0:-1]):
            if not self.can_monomer_bond_right(monomer):
                errors.append('Monomeric form {} must be able to bond to the right'.format(i_monomer + 1))

        if self.circular and not self.can_monomer_bond_right(self.seq[-1]):
            errors.append('Monomeric form {} must be able to bond to the right'.format(len(self.seq)))

        # check that monomers 2 .. L can bond to left
        for i_monomer, monomer in enumerate(self.seq[1:]):
            if not self.can_monomer_bond_left(monomer):
                errors.append('Monomeric form {} must be able to bond to the left'.format(i_monomer + 2))

        if self.circular and not self.can_monomer_bond_left(self.seq[0]):
            errors.append('Monomeric form {} must be able to bond to the left'.format(1))

        # left/right backbone/monomer atoms same length
        if len(self.bond.l_bond_atoms) != len(self.bond.r_bond_atoms):
            errors.append('Number or left and right bond atoms must be equal')

        n_bond_left_atoms = sum([1 for atom in self.bond.l_bond_atoms if atom.molecule == Backbone])
        n_bond_right_atoms = sum([1 for atom in self.bond.r_bond_atoms if atom.molecule == Backbone])

        for i_monomer, monomer in enumerate(self.seq):
            if len(self.backbone.monomer_bond_atoms) < len(monomer.backbone_bond_atoms):
                errors.append(('Number of monomer-backbone atoms for monomer {} '
                               'must be less than or equal to the number of backbone-monomer atoms').format(i_monomer + 1))

        for i_monomer, (monomer_l, monomer_r) in enumerate(zip(self.seq[0:-1], self.seq[1:])):
            if len(monomer_l.r_bond_atoms) + n_bond_right_atoms != len(monomer_r.l_bond_atoms) + n_bond_left_atoms:
                errors.append('Number of right and left bond atoms must be equal for monomers {}-{}'.format(i_monomer + 1, i_monomer + 2))

        if self.circular and \
                len(self.seq[-1].r_bond_atoms) + n_bond_right_atoms != len(self.seq[0].l_bond_atoms) + n_bond_left_atoms:
            errors.append('Number of right and left bond atoms must be equal for monomers {}, {}'.format(len(self.seq), 1))

        for i_crosslink, crosslink in enumerate(self.crosslinks):
            if len(crosslink.l_bond_atoms) != len(crosslink.r_bond_atoms):
                errors.append('Number of right and left bond atoms must be equal for crosslink {}'.format(i_crosslink + 1))

        crosslinks = list(self.crosslinks)
        for crosslink in crosslinks:
            for other_crosslink in crosslinks[1:]:
                if crosslink.is_equal(other_crosslink):
                    errors.append('Crosslink {} cannot be repeated'.format(str(crosslink)))

        # return errors
        return errors

    def can_monomer_bond_left(self, monomer):
        """ Check if monomeric form can bond to the left

        Args:
            monomer (:obj:`Monomer`): monomeric form

        Returns:
            :obj:`bool`: :obj:`True`, if the monomeric form can bond to the left
        """
        return len(monomer.l_bond_atoms) > 0 or \
            (len(monomer.backbone_bond_atoms) > 0
             and len(self.backbone.monomer_bond_atoms) > 0
             and len(self.bond.l_bond_atoms) > 0
             and self.bond.l_bond_atoms[0].molecule == Backbone)

    def can_monomer_bond_right(self, monomer):
        """ Check if monomeric form can bond to right

        Args:
            monomer (:obj:`Monomer`): monomeric form

        Returns:
            :obj:`bool`: :obj:`True`, if the monomeric form can bond to the right
        """
        return len(monomer.r_bond_atoms) > 0 or \
            (len(monomer.backbone_bond_atoms) > 0
             and len(self.backbone.monomer_bond_atoms) > 0
             and len(self.bond.r_bond_atoms) > 0
             and self.bond.r_bond_atoms[0].molecule == Backbone)

    def get_monomer_counts(self):
        """ Get the frequency of each monomeric form within the biopolymer

        Returns:
            :obj:`dict`: dictionary that maps monomeric forms to their counts
        """
        return self.seq.get_monomer_counts()

    def get_major_micro_species(self, ph, major_tautomer=False, dearomatize=False):
        """ Get the major protonation and tautomerization state

        Args:
            ph (:obj:`float`): pH
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule

        Returns:
            :obj:`openbabel.OBMol`: major protonation and tautomerization state
        """
        if not self.seq:
            return None

        smiles = self.export('smiles')
        smiles = get_major_micro_species(smiles, 'smiles', 'smiles',
                                         ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)
        mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        assert conv.SetInFormat('smiles')
        assert conv.ReadString(mol, smiles)
        return mol

    def get_structure(self, include_all_hydrogens=False):
        """ Get an Open Babel molecule of the structure

        Args:
            include_all_hydrogens (:obj:`bool`, optional): if :obj:`True`, explicitly include all
                hydrogens

        Returns:
            :obj:`tuple`:
                * :obj:`openbabel.OBMol`: Open Babel molecule of the structure
                * :obj:`dict` of :obj:`dict`: dictionary which maps indices (1-based) of monomeric forms
                    to dictionaries which map types of components of monomeric forms ('monomer' or 'backbone')
                    to dictionaries which map indices (1-based) of atoms to atoms (instances of :obj:`openbabel.OBAtom`)
        """
        if not self.seq:
            return (None, None)

        mol = openbabel.OBMol()
        n_atom = 0

        # join molecules and get bonded and displaced atoms
        atom_map = {}
        atoms = []
        n_atoms = []
        for i_monomer, monomer in enumerate(self.seq):
            if monomer.structure is None:
                monomer_structure = openbabel.OBMol()
            else:
                monomer_structure = openbabel.OBMol(monomer.structure)
            if self.backbone.structure is None:
                backbone_structure = openbabel.OBMol()
            else:
                backbone_structure = openbabel.OBMol(self.backbone.structure)

            n_atoms.append(n_atom)
            atom_map[i_monomer + 1] = {'monomer': {}, 'backbone': {}}

            mol += monomer_structure
            atom_map[i_monomer + 1]['monomer'] = {i_atom + 1: mol.GetAtom(n_atom + i_atom + 1)
                                                  for i_atom in range(monomer_structure.NumAtoms())}

            if monomer.backbone_bond_atoms and backbone_structure:
                mol += backbone_structure
                atom_map[i_monomer + 1]['backbone'] = {i_atom + 1:
                                                       mol.GetAtom(n_atom + monomer_structure.NumAtoms() + i_atom + 1)
                                                       for i_atom in range(backbone_structure.NumAtoms())}

            monomer_atom_attrs = ['monomer', [
                ['backbone_bond_atoms', monomer.backbone_bond_atoms],
                ['backbone_displaced_atoms', monomer.backbone_displaced_atoms],
            ]]

            backbone_atom_attrs = ['backbone', [
                ['monomer_bond_atoms', self.backbone.monomer_bond_atoms],
                ['monomer_displaced_atoms', self.backbone.monomer_displaced_atoms]]]
            if not monomer.backbone_bond_atoms:
                backbone_atom_attrs[1][0][1] = []
                backbone_atom_attrs[1][1][1] = []

            left_atom_attrs = ['left', [
                ['l_bond_atoms', []],
                ['l_displaced_atoms', []]]]
            if all(atom.position for atom in self.bond.l_bond_atoms):
                left_atom_attrs[1][0][1] = self.bond.l_bond_atoms
            if all(atom.position for atom in self.bond.l_displaced_atoms):
                left_atom_attrs[1][1][1] = self.bond.l_displaced_atoms
            if monomer.l_bond_atoms:
                left_atom_attrs[1][0][1] = monomer.l_bond_atoms
            if monomer.l_displaced_atoms:
                left_atom_attrs[1][1][1] = monomer.l_displaced_atoms

            right_atom_attrs = ['right', [
                ['r_bond_atoms', []],
                ['r_displaced_atoms', []]]]
            if all(atom.position for atom in self.bond.r_bond_atoms):
                right_atom_attrs[1][0][1] = self.bond.r_bond_atoms
            if all(atom.position for atom in self.bond.r_displaced_atoms):
                right_atom_attrs[1][1][1] = self.bond.r_displaced_atoms
            if monomer.r_bond_atoms:
                right_atom_attrs[1][0][1] = monomer.r_bond_atoms
            if monomer.r_displaced_atoms:
                right_atom_attrs[1][1][1] = monomer.r_displaced_atoms

            subunit_atoms = {}
            for type, attrs in [monomer_atom_attrs, backbone_atom_attrs,
                                left_atom_attrs, right_atom_attrs]:
                subunit_atoms[type] = {}
                for attr in attrs:
                    subunit_atoms[type][attr[0]] = []
                    for atom_md in attr[1]:
                        if atom_md.molecule == Monomer:
                            atom = mol.GetAtom(n_atom + atom_md.position)
                        else:
                            atom = mol.GetAtom(n_atom + atom_md.position + monomer_structure.NumAtoms())
                        subunit_atoms[type][attr[0]].append([atom, atom_md.molecule, atom_md.position, atom_md.element, atom_md.charge])
            atoms.append(subunit_atoms)

            n_atom += monomer_structure.NumAtoms()
            if backbone_structure:
                n_atom += backbone_structure.NumAtoms()

        bonding_hydrogens = []

        for i_monomer, subunit_atoms in enumerate(atoms):
            for residue_atoms in subunit_atoms.values():
                for type_atoms in residue_atoms.values():
                    for i_atom_el, atom_el in enumerate(type_atoms):
                        if atom_el[-2] == 'H':
                            type_atoms[i_atom_el] = (
                                get_hydrogen_atom(atom_el[0], bonding_hydrogens, i_monomer),
                                i_monomer + 1, atom_el[1], atom_el[2], atom_el[4])
                        else:
                            type_atoms[i_atom_el] = (
                                atom_el[0], i_monomer + 1, atom_el[1],
                                atom_el[2], atom_el[4])

        crosslinks_atoms = []
        for crosslink in self.crosslinks:
            crosslink_atoms = {}
            crosslinks_atoms.append(crosslink_atoms)
            for atom_type in ['l_bond_atoms', 'r_bond_atoms', 'l_displaced_atoms', 'r_displaced_atoms']:
                crosslink_atoms[atom_type] = []
                for atom_md in getattr(crosslink, atom_type):
                    atom = mol.GetAtom(n_atoms[atom_md.monomer - 1] + atom_md.position)
                    if atom_md.element == 'H':
                        atom = get_hydrogen_atom(atom, bonding_hydrogens, atom_md.monomer - 1)
                    crosslink_atoms[atom_type].append((atom, atom_md.monomer, atom_md.position))

        # bond monomeric forms to backbones
        for subunit_atoms in atoms:
            self._bond_monomer_backbone(mol, subunit_atoms, atom_map)

        # bond left/right pairs of subunits
        for left_atoms, right_atoms in zip(atoms[0:-1], atoms[1:]):
            self._bond_subunits(mol, left_atoms, right_atoms, atom_map)

        # circularity
        if self.circular:
            self._bond_subunits(mol, atoms[-1], atoms[0], atom_map)

        # crosslinks
        i_crosslinks_bond_atoms = []
        for crosslink_atoms in crosslinks_atoms:
            i_crosslinks_bond_atoms.append(self._form_crosslink(mol, crosslink_atoms, atom_map))

        # include all hydrogens
        if include_all_hydrogens:
            mol.AddHydrogens()

        # atom map
        for monomer in atom_map.values():
            for motif in monomer.keys():
                for i_atom, atom in list(monomer[motif].items()):
                    monomer[motif][i_atom] = atom.GetIdx()

        # return molecule
        return (mol, atom_map)

    def _bond_monomer_backbone(self, mol, subunit_atoms, atom_map):
        """ Bond a monomeric form to a backbone

        Args:
            mol (:obj:`openbabel.OBMol`): molecule with a monomeric form and backbone
            subunit_atoms (:obj:`dict`): dictionary of atoms in monomeric form and backbone to bond
            atom_map (:obj:`dict`): dictionary which maps indices (1-based) of monomeric forms
                to dictionaries which map types of components of monomeric forms ('monomer' or 'backbone')
                to dictionaries which map indices (1-based) of atoms to atoms (instances of :obj:`openbabel.OBAtom`)
        """
        for atom, i_monomer, monomer_type, i_atom, _ in subunit_atoms['monomer']['backbone_displaced_atoms']:
            if atom is not None:
                atom_2 = atom_map[i_monomer][monomer_type.__name__.lower()].get(i_atom, None)
                if atom_2 and atom_2.GetIdx() == atom.GetIdx():
                    atom_map[i_monomer][monomer_type.__name__.lower()].pop(i_atom)
                assert mol.DeleteAtom(atom, True)

        for atom, i_monomer, monomer_type, i_atom, _ in subunit_atoms['backbone']['monomer_displaced_atoms']:
            if atom is not None:
                atom_2 = atom_map[i_monomer][monomer_type.__name__.lower()].get(i_atom, None)
                if atom_2 and atom_2.GetIdx() == atom.GetIdx():
                    atom_map[i_monomer][monomer_type.__name__.lower()].pop(i_atom)
                assert mol.DeleteAtom(atom, True)

        for (monomer_atom, _, _, _, monomer_atom_charge), (backbone_atom, _, _, _, backbone_atom_charge) in zip(
                subunit_atoms['monomer']['backbone_bond_atoms'],
                subunit_atoms['backbone']['monomer_bond_atoms']):
            bond = openbabel.OBBond()
            bond.SetBegin(monomer_atom)
            bond.SetEnd(backbone_atom)
            bond.SetBondOrder(1)
            assert mol.AddBond(bond)

            if monomer_atom_charge:
                monomer_atom.SetFormalCharge(monomer_atom.GetFormalCharge() + monomer_atom_charge)
            if backbone_atom_charge:
                backbone_atom.SetFormalCharge(backbone_atom.GetFormalCharge() + backbone_atom_charge)

    def _bond_subunits(self, mol, left_atoms, right_atoms, atom_map):
        """  Bond a left/right pair of subunits

        Args:
            mol (:obj:`openbabel.OBMol`): molecule with left and right subunits
            left_atoms (:obj:`dict`): dictionary of atoms in left subunit to bond
            right_atoms (:obj:`dict`): dictionary of atoms in right subunit to bond
            atom_map (:obj:`dict`): dictionary which maps indices (1-based) of monomeric forms
                to dictionaries which map types of components of monomeric forms ('monomer' or 'backbone')
                to dictionaries which map indices (1-based) of atoms to atoms (instances of :obj:`openbabel.OBAtom`)
        """
        for atom, i_monomer, monomer_type, i_atom, _ in left_atoms['right']['r_displaced_atoms']:
            if atom is not None:
                atom_2 = atom_map[i_monomer][monomer_type.__name__.lower()].get(i_atom, None)
                if atom_2 and atom_2.GetIdx() == atom.GetIdx():
                    atom_map[i_monomer][monomer_type.__name__.lower()].pop(i_atom)
                assert mol.DeleteAtom(atom, True)

        for atom, i_monomer, monomer_type, i_atom, _ in right_atoms['left']['l_displaced_atoms']:
            if atom is not None:
                atom_2 = atom_map[i_monomer][monomer_type.__name__.lower()].get(i_atom, None)
                if atom_2 and atom_2.GetIdx() == atom.GetIdx():
                    atom_map[i_monomer][monomer_type.__name__.lower()].pop(i_atom)
                assert mol.DeleteAtom(atom, True)

        for (l_atom, _, _, _, l_atom_charge), (r_atom, _, _, _, r_atom_charge) in \
            zip(left_atoms['right']['r_bond_atoms'],
                right_atoms['left']['l_bond_atoms']):
            bond = openbabel.OBBond()
            bond.SetBegin(l_atom)
            bond.SetEnd(r_atom)
            bond.SetBondOrder(1)
            assert mol.AddBond(bond)

            if l_atom_charge:
                l_atom.SetFormalCharge(l_atom.GetFormalCharge() + l_atom_charge)
            if r_atom_charge:
                r_atom.SetFormalCharge(r_atom.GetFormalCharge() + r_atom_charge)

    def _form_crosslink(self, mol, atoms, atom_map):
        """ Form bond(s) for a crosslink

        Args:
            mol (:obj:`openbabel.OBMol`): molecule
            atoms (:obj:`dict`): dictionary of atoms to bond and delete
            atom_map (:obj:`dict`): dictionary which maps indices (1-based) of monomeric forms
                to dictionaries which map types of components of monomeric forms ('monomer' or 'backbone')
                to dictionaries which map indices (1-based) of atoms to atoms (instances of :obj:`openbabel.OBAtom`)
        """
        for (l_atom, _, _), (r_atom, _, _) in zip(atoms['l_bond_atoms'], atoms['r_bond_atoms']):
            bond = openbabel.OBBond()
            bond.SetBegin(l_atom)
            bond.SetEnd(r_atom)
            bond.SetBondOrder(1)
            assert mol.AddBond(bond)

        for atom, i_monomer, i_atom in itertools.chain(atoms['l_displaced_atoms'], atoms['r_displaced_atoms']):
            if atom is not None:
                atom_2 = atom_map[i_monomer]['monomer'].get(i_atom, None)
                if atom_2 and atom_2.GetIdx() == atom.GetIdx():
                    atom_map[i_monomer]['monomer'].pop(i_atom)
                assert mol.DeleteAtom(atom, True)

    def export(self, format, include_all_hydrogens=False, options=()):
        """ Export structure to format

        Args:
            format (:obj:`str`): format
            include_all_hydrogens (:obj:`bool`, optional): if :obj:`True`, explicitly
                include all hydrogens
            options (:obj:`list` of :obj:`str`, optional): export options

        Returns:
            :obj:`str`: format representation of structure
        """
        structure = self.get_structure(include_all_hydrogens=include_all_hydrogens)[0]
        if structure is None:
            return None
        else:
            return OpenBabelUtils.export(structure, format, options=options)

    def get_formula(self):
        """ Get the chemical formula

        Returns:
            :obj:`EmpiricalFormula`: chemical formula
        """
        backbone_formula = self.backbone.get_formula()
        n_backbone = 0

        formula = EmpiricalFormula()
        for monomer, count in self.get_monomer_counts().items():
            formula += monomer.get_formula() * count
            if monomer.backbone_bond_atoms:
                n_backbone += count
            for atom in monomer.backbone_displaced_atoms:
                formula[atom.element] -= count

        for monomer in self.seq[0:-1]:
            for atom in monomer.r_displaced_atoms:
                formula[atom.element] -= 1
        for monomer in self.seq[1:]:
            for atom in monomer.l_displaced_atoms:
                formula[atom.element] -= 1
        if self.circular:
            for atom in self.seq[-1].r_displaced_atoms:
                formula[atom.element] -= 1
            for atom in self.seq[0].l_displaced_atoms:
                formula[atom.element] -= 1

        # crosslinks
        for crosslink in self.crosslinks:
            for atom in itertools.chain(crosslink.l_displaced_atoms, crosslink.r_displaced_atoms):
                formula[atom.element] -= 1

        return formula \
            + backbone_formula * n_backbone \
            + self.bond.get_formula(none_position=False) * (len(self.seq) - (1 - self.circular))

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
        backbone_charge = self.backbone.get_charge()
        n_backbone = 0

        charge = 0
        for monomer, count in self.get_monomer_counts().items():
            charge += monomer.get_charge() * count
            if monomer.backbone_bond_atoms:
                n_backbone += count
            for atom in monomer.backbone_bond_atoms:
                charge -= atom.charge * count
            for atom in monomer.backbone_displaced_atoms:
                charge -= atom.charge * count

        for monomer in self.seq[0:-1]:
            for atom in monomer.r_displaced_atoms:
                charge -= atom.charge
        for monomer in self.seq[1:]:
            for atom in monomer.l_displaced_atoms:
                charge -= atom.charge
        if self.circular:
            for atom in self.seq[-1].r_displaced_atoms:
                charge -= atom.charge
            for atom in self.seq[0].l_displaced_atoms:
                charge -= atom.charge

        # crosslinks
        for crosslink in self.crosslinks:
            for atom in itertools.chain(crosslink.l_displaced_atoms, crosslink.r_displaced_atoms):
                charge -= atom.charge

        return charge \
            + backbone_charge * n_backbone \
            + self.bond.get_charge(none_position=False) * (len(self.seq) - (1 - self.circular))

    def __str__(self):
        """ Get a string representation of the biopolymer form

        Returns:
            :obj:`str`: string representation of the biopolymer form
        """
        alphabet_monomers = {monomer: code for code, monomer in self.alphabet.monomers.items()}
        val = ''
        for monomer in self.seq:
            code = alphabet_monomers.get(monomer, None)
            if code:
                if len(code) == 1:
                    val += code
                else:
                    val += '{' + code + '}'
            else:
                val += monomer.__str__(alphabet=self.alphabet)

        # circular
        if self.circular:
            val += ' | circular'

        # crosslinks
        for crosslink in self.crosslinks:
            atoms = []
            for atom_type in ['l_bond_atoms', 'r_bond_atoms', 'l_displaced_atoms', 'r_displaced_atoms']:
                for atom in getattr(crosslink, atom_type):
                    if atom.charge > 0:
                        charge = '+' + str(atom.charge)
                    elif atom.charge == 0:
                        charge = ''
                    else:
                        charge = str(atom.charge)
                    atoms.append('{}: {}{}{}{}'.format(atom_type[0:-1].replace('_', '-'),
                                                       atom.monomer, atom.element, atom.position, charge))

            val += ' | x-link: [{}]'.format(' | '.join(atoms))

        return val

    _grammar_filename = pkg_resources.resource_filename('bpforms', 'grammar.lark')
    with open(_grammar_filename, 'r') as file:
        _parser = lark.Lark(file.read())

    def from_str(self, string):
        """ Create biopolymer form its string representation

        Args:
            string (:obj:`str`): string representation of the biopolymer

        Returns:
            :obj:`BpForm`: biopolymer form
        """
        class ParseTreeTransformer(lark.Transformer):
            def __init__(self, bp_form):
                super(ParseTreeTransformer, self).__init__()
                self.bp_form = bp_form
                bp_form.seq.clear()
                bp_form.crosslinks.clear()
                bp_form.circular = False

            @lark.v_args(inline=True)
            def start(self, *args):
                for arg_type, arg_val in args[0::2]:
                    if arg_type in ['seq', 'circular']:
                        setattr(self.bp_form, arg_type, arg_val)
                    else:
                        set_val = getattr(self.bp_form, arg_type)
                        for other in set_val:
                            if arg_val.is_equal(other):
                                str_arg_val = str(arg_val)
                                raise ValueError(f'`{arg_type}` cannot contain {str_arg_val} multiple times')
                        set_val.add(arg_val)
                return self.bp_form

            @lark.v_args(inline=True)
            def seq(self, *seq):
                monomer_seq = []
                for monomer in seq:
                    if isinstance(monomer, Monomer):
                        monomer_seq.append(monomer)
                return ('seq', monomer_seq)

            @lark.v_args(inline=True)
            def alphabet_monomer(self, code):
                code = code.value
                if code[0] == '{' and code[-1] == '}':
                    code = code[1:-1]
                monomer = self.bp_form.alphabet.monomers.get(code, None)
                if monomer is None:
                    raise ValueError('"{}" not in alphabet'.format(code))
                return monomer

            @lark.v_args(inline=True)
            def inline_monomer(self, *args):
                kwargs = {
                    'synonyms': SynonymSet(),
                    'identifiers': IdentifierSet(),
                    'monomers_position': set(),
                    'base_monomers': set(),
                    'backbone_bond_atoms': AtomList(),
                    'backbone_displaced_atoms': AtomList(),
                    'r_bond_atoms': AtomList(),
                    'r_displaced_atoms': AtomList(),
                    'l_bond_atoms': AtomList(),
                    'l_displaced_atoms': AtomList(),
                }
                for arg in args[1:-1:2]:
                    arg_name, arg_val = arg
                    if arg_name in ['id', 'name', 'structure', 'delta_mass', 'delta_charge', 'position', 'comments']:
                        if arg_name in kwargs:
                            raise ValueError('{} attribute cannot be repeated'.format(arg_name))
                        if arg_name == 'position':
                            kwargs['start_position'], kwargs['end_position'], kwargs['monomers_position'] = arg_val
                        else:
                            kwargs[arg_name] = arg_val
                    elif arg_name in ['synonyms', 'identifiers', 'monomers_position', 'base_monomers']:
                        kwargs[arg_name].add(arg_val)
                    elif arg_name in ['backbone_bond_atoms', 'backbone_displaced_atoms',
                                      'r_bond_atoms', 'r_displaced_atoms',
                                      'l_bond_atoms', 'l_displaced_atoms']:
                        kwargs[arg_name].append(arg_val)
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
                return ('identifiers', Identifier(args[-1].value[1:-1], args[-3].value[1:-1]))

            @lark.v_args(inline=True)
            def structure(self, *args):
                return ('structure', args[-2].value)

            @lark.v_args(inline=True)
            def backbone_bond_atom(self, *args):
                return ('backbone_bond_atoms', args[-1])

            @lark.v_args(inline=True)
            def backbone_displaced_atom(self, *args):
                return ('backbone_displaced_atoms', args[-1])

            @lark.v_args(inline=True)
            def r_bond_atom(self, *args):
                return ('r_bond_atoms', args[-1])

            @lark.v_args(inline=True)
            def r_displaced_atom(self, *args):
                return ('r_displaced_atoms', args[-1])

            @lark.v_args(inline=True)
            def l_bond_atom(self, *args):
                return ('l_bond_atoms', args[-1])

            @lark.v_args(inline=True)
            def l_displaced_atom(self, *args):
                return ('l_displaced_atoms', args[-1])

            @lark.v_args(inline=True)
            def atom(self, *args):
                atom = Atom(Monomer, args[0].value, position=int(float(args[1].value)))
                if len(args) >= 3:
                    atom.charge = int(float(args[2].value))
                return atom

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
                monomers_position = []

                before_dash = True
                for arg in args:
                    if isinstance(arg, lark.lexer.Token):
                        if arg.type == 'DASH':
                            before_dash = False
                        if arg.type == 'INT':
                            if before_dash:
                                start_position = int(float(arg.value))
                            else:
                                end_position = int(float(arg.value))
                    elif isinstance(arg, list):
                        monomers_position = arg

                return ('position', (start_position, end_position, monomers_position))

            @lark.v_args(inline=True)
            def monomers_position(self, *args):
                monomers = []
                for arg in args:
                    if arg.type in ['CHAR', 'CHARS']:
                        monomer = self.bp_form.alphabet.monomers.get(arg.value, None)
                        if monomer is None:
                            raise ValueError('"{}" not in alphabet'.format(arg.value))
                        monomers.append(monomer)
                return monomers

            @lark.v_args(inline=True)
            def base_monomer(self, *args):
                code = args[-2].value
                if code[0] == '"' and code[-1] == '"':
                    code = code[1:-1]
                monomer = self.bp_form.alphabet.monomers.get(code, None)
                if monomer is None:
                    raise ValueError('"{}" not in alphabet'.format(code))
                return ('base_monomers', monomer)

            @lark.v_args(inline=True)
            def comments(self, *args):
                return ('comments', args[-1].value[1:-1])

            @lark.v_args(inline=True)
            def crosslink(self, *args):
                bond = Bond()

                if len(args) <= 4:
                    args = []
                elif isinstance(args[3], tuple):
                    args = args[3:-1:2]
                else:
                    args = args[4:-1:2]

                for atom_type, atom in args:
                    atom_type_list = getattr(bond, atom_type)
                    atom_type_list.append(atom)

                return ('crosslinks', bond)

            @lark.v_args(inline=True)
            def crosslink_atom(self, *args):
                atom_type = args[0]
                monomer = int(float(args[2].value))
                element = args[3].value
                position = int(float(args[4].value))
                if len(args) >= 6:
                    charge = int(float(args[5].value))
                else:
                    charge = 0
                return atom_type, Atom(Monomer, monomer=monomer, element=element, position=position, charge=charge)

            @lark.v_args(inline=True)
            def crosslink_atom_type(self, *args):
                return args[0].value.replace('-', '_') + 's'

            @lark.v_args(inline=True)
            def circular(self, *args):
                return ('circular', True)

        tree = self._parser.parse(string)
        parse_tree_transformer = ParseTreeTransformer(self)
        return parse_tree_transformer.transform(tree)

    def get_canonical_seq(self, monomer_codes=None):
        """ Get IUPAC/IUBMB representation of a polymer with bases represented by the character codes
        of their parent monomers (e.g. methyl-2-adenosine is represented by 'A')

        Args:
            monomer_codes (:obj:`dict`, optional): dictionary that maps monomers to their codes

        Returns:
            :obj:`str`: IUPAC/IUBMB representation of a polymer
        """

        if monomer_codes is None:
            monomer_codes = {monomer: code for code, monomer in self.alphabet.monomers.items()}

        seq = ''
        for monomer in self.seq:
            seq += monomer.get_canonical_code(default_code=self.DEFAULT_FASTA_CODE, monomer_codes=monomer_codes)

        return seq

    def get_image(self, monomer_color=0x000000, backbone_color=0xff0000,
                  left_right_bond_color=0x00ff00, crosslink_bond_color=0x0000ff,
                  include_all_hydrogens=True, show_atom_nums=False,
                  width=200, height=200, image_format='svg', include_xml_header=True):
        """ Get image

        Args:
            monomer_color (:obj:`int`, optional): color to paint atoms involved in monomeric forms
            backbone_color (:obj:`int`, optional): color to paint atoms involved in backbones
            left_right_bond_color (:obj:`int`, optional): color to paint atoms involved in bond with monomeric form to left
            crosslink_bond_color (:obj:`int`, optional): color to paint atoms involved in crosslinks
            include_all_hydrogens (:obj:`bool`, optional): if :obj:`True`, show all hydrogens
            show_atom_nums (:obj:`bool`, optional): if :obj:`True`, show the numbers of the atoms
            width (:obj:`int`, optional): width in pixels
            height (:obj:`int`, optional): height in pixels
            image_format (:obj:`str`, optional): format of generated image {emf, eps, jpeg, msbmp, pdf, png, or svg}
            include_xml_header (:obj:`bool`, optional): if :obj:`True`, include XML header at the beginning of the SVG

        Returns:
            :obj:`object`: image
        """
        mol, atom_map = self.get_structure(include_all_hydrogens=include_all_hydrogens)
        el_table = openbabel.OBElementTable()

        atom_labels = []
        atom_sets = {}
        atom_sets[backbone_color] = {'positions': [], 'elements': [], 'color': backbone_color}
        bond_sets = {}
        bond_sets[left_right_bond_color] = {'positions': [], 'elements': [], 'color': left_right_bond_color}
        bond_sets[crosslink_bond_color] = {'positions': [], 'elements': [], 'color': crosslink_bond_color}

        # label monomers
        for i_monomer in atom_map.keys():
            i_atoms = sorted(atom_map[i_monomer]['monomer'].keys())
            if i_atoms:
                i_atom = atom_map[i_monomer]['monomer'][i_atoms[0]]
                atom = mol.GetAtom(i_atom)
                atom_labels.append({'position': i_atom,
                                    'element': el_table.GetSymbol(atom.GetAtomicNum()),
                                    'label': str(i_monomer),
                                    'color': monomer_color})

        # color backbones
        for i_monomer in atom_map.keys():
            for i_atom in atom_map[i_monomer]['backbone'].values():
                atom = mol.GetAtom(i_atom)
                atom_sets[backbone_color]['positions'].append(i_atom)
                atom_sets[backbone_color]['elements'].append(el_table.GetSymbol(atom.GetAtomicNum()))

        # color left/right bonds
        for i_monomer, (l_monomer, r_monomer) in enumerate(zip(self.seq[0:-1], self.seq[1:])):
            l_bond_atoms = l_monomer.r_bond_atoms or self.bond.r_bond_atoms
            r_bond_atoms = r_monomer.l_bond_atoms or self.bond.l_bond_atoms
            self._add_bonds_to_set(mol,
                                   l_bond_atoms,
                                   r_bond_atoms,
                                   atom_map, bond_sets[left_right_bond_color],
                                   el_table,
                                   i_monomer_1=i_monomer + 1,
                                   i_monomer_2=i_monomer + 2)

        # color crosslink bonds
        for crosslink in self.crosslinks:
            self._add_bonds_to_set(mol, crosslink.l_bond_atoms,
                                   crosslink.r_bond_atoms,
                                   atom_map, bond_sets[crosslink_bond_color],
                                   el_table)

        # remove unused atoms and bond sets
        if not atom_sets[backbone_color]['positions']:
            atom_sets.pop(backbone_color)
        if not bond_sets[left_right_bond_color]['positions']:
            bond_sets.pop(left_right_bond_color)
        if not bond_sets[crosslink_bond_color]['positions']:
            bond_sets.pop(crosslink_bond_color)

        cml = OpenBabelUtils.export(mol, 'cml')

        # draw molecule
        return draw_molecule(cml, 'cml', image_format=image_format,
                             atom_labels=atom_labels, atom_sets=atom_sets.values(), bond_sets=bond_sets.values(),
                             show_atom_nums=show_atom_nums,
                             width=width, height=height, include_xml_header=include_xml_header)

    def _add_bonds_to_set(self, mol, r_bond_atoms, l_bond_atoms, atom_map, bond_set, el_table,
                          i_monomer_1=None, i_monomer_2=None):
        for md_atom_1, md_atom_2 in zip(r_bond_atoms, l_bond_atoms):
            if i_monomer_1 is None:
                i_monomer_1 = md_atom_1.monomer
            else:
                i_monomer_1 = i_monomer_1

            if i_monomer_2 is None:
                i_monomer_2 = md_atom_2.monomer
            else:
                i_monomer_2 = i_monomer_2

            type_1 = md_atom_1.molecule.__name__.lower()
            type_2 = md_atom_2.molecule.__name__.lower()
            i_atom_1 = atom_map[i_monomer_1][type_1].get(md_atom_1.position, None)
            i_atom_2 = atom_map[i_monomer_2][type_2].get(md_atom_2.position, None)

            atom_1 = mol.GetAtom(i_atom_1)
            atom_2 = mol.GetAtom(i_atom_2)

            bond_set['positions'].append(
                [i_atom_1, i_atom_2])
            bond_set['elements'].append(
                [el_table.GetSymbol(atom_1.GetAtomicNum()),
                 el_table.GetSymbol(atom_2.GetAtomicNum())])


class BpFormFeature(object):
    """ A region (start and end positions) of a BpForm

    Attributes:
        form (:obj:`BpForm`): biopolymer form
        start_position (:obj:`int`): start position (1-base)
        end_position (:obj:`int`): end position (1-based)
    """

    def __init__(self, form, start_position, end_position):
        """
        Args:
            form (:obj:`BpForm`): biopolymer form
            start_position (:obj:`int`): start position (1-base)
            end_position (:obj:`int`): end position (1-based)
        """
        self.form = form
        self.start_position = start_position
        self.end_position = end_position

    @property
    def form(self):
        """ Get the biopolymer form

        Returns:
            :obj:`BpForm`: biopolymer form
        """
        return self._form

    @form.setter
    def form(self, value):
        """ Set the biopolymer form

        Args:
            value (:obj:`BpForm`): biopolymer form

        Raises:
            :obj:`ValueError`: form is not an instance of :obj:`BpForm`
        """
        if value is not None and not isinstance(value, BpForm):
            raise ValueError('`form` must be an instance of `BpForm` or None')

        if hasattr(self, '_form'):
            if value != self._form:
                if self._form is not None:
                    old_form = self._form
                    self._form = None
                    if self in old_form.features:
                        old_form.features.remove(self)
                if value is not None:
                    self._form = value
                    if self not in value.features:
                        value.features.add(self)
        else:
            self._form = value
            if value is not None:
                if self not in value.features:
                    value.features.add(self)

    @property
    def start_position(self):
        """ Get the start position

        Returns:
            :obj:`int`: start position
        """
        return self._start_position

    @start_position.setter
    def start_position(self, value):
        """ Set the start position

        Args:
            value (:obj:`int`): start position

        Raises:
            :obj:`ValueError`: start position is not a non-negative integer
        """
        if int(value) != value or value < 0:
            raise ValueError('`start_position` must be a non-negative integer')
        self._start_position = int(value)

    @property
    def end_position(self):
        """ Get the end position

        Returns:
            :obj:`int`: end position
        """
        return self._end_position

    @end_position.setter
    def end_position(self, value):
        """ Set the end position

        Args:
            value (:obj:`int`): end position

        Raises:
            :obj:`ValueError`: end position is not a non-negative integer
        """
        if int(value) != value or value < 0:
            raise ValueError('`end_position` must be a non-negative integer')
        self._end_position = int(value)


class BpFormFeatureSet(set):
    """ Set of features

    Attributes:
        form (:obj:`BpForm`): form
    """

    def __init__(self, form):
        """
        Args:
            form (:obj:`BpForm`): form
        """
        super(BpFormFeatureSet, self).__init__()
        self.form = form

    @property
    def form(self):
        """ Get the biopolymer form

        Returns:
            :obj:`BpForm`: biopolymer form
        """
        return self._form

    @form.setter
    def form(self, value):
        """ Set the biopolymer form

        Args:
            value (:obj:`BpForm`): biopolymer form

        Raises:
            :obj:`ValueError`: form is not an instance of :obj:`BpForm`
        """
        if not isinstance(value, BpForm):
            raise ValueError('`form` must be an instance of `BpForm`')

        if hasattr(self, '_form'):
            raise ValueError('`form` cannot be set outside constructor')

        self._form = value

    def add(self, feature):
        """ Add a feature

        Args:
            feature (:obj:`BpFormFeature`): feature

        Raises:
            :obj:`ValueError`: if the `feature` is not an instance of `BpFormFeature`
        """
        if not isinstance(feature, BpFormFeature):
            raise ValueError('`feature` must be an instance of `BpFormFeature`')
        if feature not in self:
            super(BpFormFeatureSet, self).add(feature)
            if feature.form != self.form:
                feature.form = self.form

    def remove(self, feature):
        """ Remove a feature

        Args:
            feature (:obj:`BpFormFeature`): feature
        """
        super(BpFormFeatureSet, self).remove(feature)
        if feature.form != None:
            feature.form = None

    def update(self, features):
        """ Add a set of features

        Args:
            features (iterable of :obj:`BpFormFeature`): features
        """
        for feature in features:
            self.add(feature)

    def symmetric_difference_update(self, other):
        """ Remove common elements with other and add elements from other not in self

        Args:
            other (:obj:`BpFormFeatureSet`): other set of features
        """
        for o in other:
            if o in self:
                self.remove(o)
            else:
                self.add(o)


class BpFormsWarning(UserWarning):
    """ BpForms warning """
    pass


def get_hydrogen_atom(parent_atom, bonding_hydrogens, i_monomer):
    """ Get a hydrogen atom attached to a parent atom

    Args:
        parent_atom (:obj:`openbabel.OBAtom`): parent atom
        bonding_hydrogens (:obj:`list`): hydrogens that have already been gotten
        i_monomer (:obj:`int`): index of parent monomer in sequence

    Returns:
        :obj:`openbabel.OBAtom`: hydrogen atom
    """
    for other_atom in openbabel.OBAtomAtomIter(parent_atom):
        if other_atom.GetAtomicNum() == 1:
            tmp = (i_monomer, other_atom.GetIdx())
            if tmp not in bonding_hydrogens:  # hydrogen
                bonding_hydrogens.append(tmp)
                return other_atom
    return None
