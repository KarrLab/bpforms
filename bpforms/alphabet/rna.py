""" Alphabet and BpForm to represent modified RNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import (Alphabet, AlphabetBuilder, Monomer, MonomerSequence, Backbone,
                          Bond, Atom, BpForm, Identifier, IdentifierSet, SynonymSet,
                          BpFormsWarning)
from wc_utils.util.chem import EmpiricalFormula, get_major_micro_species
import bs4
import csv
import io
import openbabel
import os.path
import pkg_resources
import requests
import requests_cache
import warnings

filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'rna.yml'))
rna_alphabet = Alphabet().from_yaml(filename)
# :obj:`Alphabet`: Alphabet for RNA nucleosides

canonical_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'rna.canonical.yml'))
canonical_rna_alphabet = Alphabet().from_yaml(canonical_filename)
# :obj:`Alphabet`: Alphabet for canonical RNA nucleosides


class RnaAlphabetBuilder(AlphabetBuilder):
    """ Build RNA alphabet from MODOMICS """

    INDEX_ENDPOINT = 'http://modomics.genesilico.pl/modifications/'
    INDEX_ASCII_ENDPOINT = 'http://modomics.genesilico.pl/modifications/?base=all&type=all&display_ascii=Display+as+ASCII'
    ENTRY_ENDPOINT = 'http://modomics.genesilico.pl/modifications/{}/'
    MAX_RETRIES = 5

    def run(self, ph=None, major_tautomer=False, path=filename):
        """ Build alphabet and, optionally, save to YAML file

        Args:
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            path (:obj:`str`, optional): path to save alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        return super(RnaAlphabetBuilder, self).run(ph=ph, major_tautomer=major_tautomer, path=path)

    def build(self, ph=None, major_tautomer=False):
        """ Build alphabet

        Args:
            ph (:obj:`float`, optional): pH at which to calculate the major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer

        Returns:
            :obj:`Alphabet`: alphabet
        """
        # initialize alphabet
        alphabet = Alphabet()

        # create canonical monomeric forms
        alphabet.from_yaml(canonical_filename)
        alphabet.id = 'rna'
        alphabet.name = 'MODOMICS RNA nucleosides'
        alphabet.description = ('The four canonical RNA nucleosides, plus the non-canonical RNA nucleosides in '
                                '<a href="http://modomics.genesilico.pl/modifications">MODOMICS</a>')

        # create requests session
        cache_name = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'rna'))
        session = requests_cache.core.CachedSession(cache_name, backend='sqlite', expire_after=None)
        session.mount('http://', requests.adapters.HTTPAdapter(max_retries=self.MAX_RETRIES))

        # get originating monomeric forms
        ascii_response = session.get(self.INDEX_ASCII_ENDPOINT)
        ascii_response.raise_for_status()
        stream = io.StringIO(ascii_response.text)
        stream.readline()
        reader = csv.reader(stream, delimiter='\t', quoting=csv.QUOTE_NONE)
        base_monomer_short_names = {}
        for monomer in reader:
            for i in range(1, 4):
                short_name = monomer[i]
                if short_name:
                    break
            if short_name in ['A', 'C', 'G', 'U']:
                continue
            base_monomer_short_names[short_name] = monomer[i + 2]

        # get index of nucleosides
        response = session.get(self.INDEX_ENDPOINT)
        response.raise_for_status()

        # get individual nucleosides and create monomeric forms
        doc = bs4.BeautifulSoup(response.text, 'html.parser')
        table = doc.find('table', {'class': 'datagrid'})
        tbody = table.find('tbody')
        mods = tbody.find_all('tr')
        monomer_short_names = {}
        invalid_nucleosides = []
        for i_mod, mod in enumerate(mods):
            if i_mod >= self._max_monomers:
                break

            cells = mod.find_all('td')
            new_nomenclature = cells[0].text
            name = cells[1].text
            short_name = cells[2].text
            abbrev = cells[3].text

            id = short_name

            chars = new_nomenclature
            if not chars:
                chars = id

            synonyms = SynonymSet()
            if abbrev:
                synonyms.add(abbrev)

            identifiers = IdentifierSet()
            identifiers.add(Identifier('modomics.short_name', short_name))
            if new_nomenclature:
                identifiers.add(Identifier('modomics.new_nomenclature', new_nomenclature))

            structure, more_identifiers = self.get_nucleoside_details(id, session)
            identifiers.update(more_identifiers)

            monomer = Monomer(
                id=chars,
                name=name,
                synonyms=synonyms,
                identifiers=identifiers,
                structure=structure,
            )

            if not monomer.structure:
                continue
            if '*' in monomer.export('smiles'):
                continue
            if ' (base)' in new_nomenclature:
                continue
            if ' (cap)' in name or ' cap' in name:
                continue
            if '-CoA)' in name:
                continue
            if new_nomenclature and new_nomenclature[-1] == 'N':
                continue

            conv = openbabel.OBConversion()
            assert conv.SetOutFormat('smi'), 'Unable to set format to SMILES'
            assert conv.SetInFormat('smi')
            conv.SetOptions('c', conv.OUTOPTIONS)

            smiles = conv.WriteString(monomer.structure).partition('\t')[0]
            if ph is not None:
                smiles = get_major_micro_species(smiles, 'smiles', 'smiles', ph, major_tautomer=major_tautomer)
            smiles_mol = openbabel.OBMol()
            conv.ReadString(smiles_mol, smiles)

            smiles = conv.WriteString(smiles_mol).partition('\t')[0]
            smiles_mol2 = openbabel.OBMol()
            conv.ReadString(smiles_mol2, smiles)

            monomer.structure = smiles_mol2

            if not self.is_valid_nucleoside(monomer):
                invalid_nucleosides.append(chars)
                continue

            alphabet.monomers[chars] = monomer
            monomer_short_names[short_name] = monomer

        if invalid_nucleosides:
            warnings.warn('The following compounds were ignored because they do not appear to be nucleosides:\n- {}'.format(
                '\n- '.join(invalid_nucleosides)), BpFormsWarning)

        for short_name, monomer in monomer_short_names.items():
            base_monomer = monomer_short_names.get(base_monomer_short_names.get(short_name, None), None)
            if base_monomer:
                monomer.base_monomers.add(base_monomer)

        # get major microspecies for each monomeric form
        self.get_major_micro_species(alphabet, ph=ph, major_tautomer=major_tautomer)

        # return alphabet
        return alphabet

    def get_nucleoside_details(self, id, session):
        """ Get the structure of a nucleoside in the MODOMICS database

        Args:
            id (:obj:`str`): id of nucleoside in MODOMICS database

        Returns:
            :obj:`openbabel.OBMol`: structure
            :obj:`IdentifierSet`: identifiers
        """
        response = session.get(self.ENTRY_ENDPOINT.format(id))
        response.raise_for_status()

        doc = bs4.BeautifulSoup(response.text, 'html.parser')

        table = doc.find(id='modification_details')
        tbody = table.find('tbody')
        rows = tbody.find_all('tr')
        mol = None
        identifiers = IdentifierSet()
        for row in rows:
            cells = row.find_all('td')
            if cells[0].text.startswith('SMILES'):
                smiles = cells[1].text
                if not smiles:
                    continue

                mol = openbabel.OBMol()
                conv = openbabel.OBConversion()
                assert conv.SetInFormat('smi')
                if not conv.ReadString(mol, smiles):
                    mol = None
                    continue

            elif cells[0].text.startswith('PubChem'):
                link = cells[1].find('a')
                identifiers.add(Identifier('pubchem.compound', link.text))

        return mol, identifiers

    def is_valid_nucleoside(self, monomer):
        """ Determine if nucleoside should be included in alphabet

        Args:
            monomer (:obj:`Monomer`): monomeric form

        Returns:
            :obj:`bool`: :obj:`True` if the monomeric form is a valid nucleoside
        """
        formula = monomer.get_formula()
        if formula.C < 9 or formula.O < 4 or formula.N < 2:
            return False

        atom_bs = []
        atom_rs = []
        for i_atom in range(1, monomer.structure.NumAtoms() + 1):
            atom = monomer.structure.GetAtom(i_atom)
            if self.is_backbone_atom(atom):
                atom_bs.append(atom)
            elif self.is_right_bond_atom(atom):
                atom_rs.append(atom)

        termini = []
        for atom_b in atom_bs:
            for atom_r in atom_rs:
                if self.is_terminus(atom_b, atom_r):
                    termini.append((atom_b, atom_r))

        if termini:
            atom_b_idx = termini[0][0].GetIdx()
            atom_r_idx = termini[0][1].GetIdx()
        else:
            if atom_bs and not atom_rs:
                atom_b_idx = atom_bs[0].GetIdx()
                atom_r_idx = None
            elif atom_rs and not atom_bs:
                atom_b_idx = None
                atom_r_idx = atom_rs[0].GetIdx()
            else:
                atom_b_idx = None
                atom_r_idx = None

        if atom_b_idx:
            monomer.backbone_bond_atoms = [Atom(Monomer, 'O', atom_b_idx)]
            monomer.backbone_displaced_atoms = [Atom(Monomer, 'H', atom_b_idx)]
        if atom_r_idx:
            monomer.right_bond_atoms = [Atom(Monomer, 'O', atom_r_idx)]
            monomer.right_displaced_atoms = [Atom(Monomer, 'H', atom_r_idx)]

        return True

    def is_terminus(self, b_atom, r_atom):
        """ Determine if a pair of atoms is a valid pair of linkage sites

        Args:
            b_atom (:obj:`openbabel.OBAtom`): potential backbone atom
            r_atom (:obj:`openbabel.OBAtom`): potential right bond atom

        Returns:
            :obj:`bool`: :obj:`True`, if the atoms are a valid pair of linkage sites
        """
        r_atom_2 = self.is_backbone_atom(b_atom)
        if not r_atom_2:
            return False

        if not self.is_right_bond_atom(r_atom):
            return False

        return r_atom_2.GetIdx() == r_atom.GetIdx()

    def is_backbone_atom(self, b_atom):
        """ Determine if an atom is a valid backbone linkage site

        Args:
            b_atom (:obj:`openbabel.OBAtom`): potential backbone atom

        Returns:
            :obj:`bool`: :obj:`True`, if the atom is a valid backbone linkage site
        """
        if b_atom.GetAtomicNum() != 8:
            return False
        if b_atom.GetFormalCharge() != 0:
            return False

        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(b_atom)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(b_atom)])
        other_atoms += [1] * (2 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        if other_atoms != [1, 6]:
            return False

        # get first C
        for other_atom in openbabel.OBAtomAtomIter(b_atom):
            if other_atom.GetAtomicNum() == 6:
                c_1 = other_atom
        if c_1.GetFormalCharge() != 0:
            return False
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(c_1)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(c_1)])
        other_atoms += [1] * (4 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        if other_atoms != [1, 1, 6, 8]:
            return False

        # get second C
        for other_atom in openbabel.OBAtomAtomIter(c_1):
            if other_atom.GetAtomicNum() == 6:
                c_2 = other_atom
        if c_2.GetFormalCharge() != 0:
            return False
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(c_2)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(c_2)])
        other_atoms += [1] * (4 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        if other_atoms != [1, 6, 6, 8]:
            return False

        # get third C
        for other_atom in openbabel.OBAtomAtomIter(c_2):
            if other_atom.GetAtomicNum() == 6 and other_atom.GetIdx() != c_1.GetIdx():
                c_3 = other_atom
        if c_3.GetFormalCharge() != 0:
            return False
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(c_3)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(c_3)])
        other_atoms += [1] * (4 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        if other_atoms != [1, 6, 6, 8]:
            return False

        # get second O
        for other_atom in openbabel.OBAtomAtomIter(c_3):
            if other_atom.GetAtomicNum() == 8:
                r_atom_2 = other_atom
        if r_atom_2.GetFormalCharge() != 0:
            return False

        return r_atom_2

    def is_right_bond_atom(self, r_atom):
        """ Determine if an atom is a valid right bond linkage site

        Args:
            b_atom (:obj:`openbabel.OBAtom`): potential right bond atom

        Returns:
            :obj:`bool`: :obj:`True`, if the atom is a valid right bond linkage site
        """
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(r_atom)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(r_atom)])
        other_atoms += [1] * (2 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        if other_atoms != [1, 6]:
            return False

        return True


class RnaForm(BpForm):
    """ RNA form """

    DEFAULT_FASTA_CODE = 'N'

    def __init__(self, seq=None, circular=False):
        """
        Args:
            seq (:obj:`MonomerSequence`, optional): sequence of monomeric forms of the DNA
            circular (:obj:`bool`, optional): if :obj:`True`, indicates that the biopolymer is circular
        """
        super(RnaForm, self).__init__(
            seq=seq, alphabet=rna_alphabet,
            backbone=Backbone(
                structure='OP([O-])([O-])=O',
                monomer_bond_atoms=[Atom(Backbone, element='P', position=2)],
                monomer_displaced_atoms=[Atom(Backbone, element='O', position=1), Atom(Backbone, element='H', position=1)]),
            bond=Bond(
                right_bond_atoms=[Atom(Monomer, element='O', position=None)],
                left_bond_atoms=[Atom(Backbone, element='P', position=2)],
                right_displaced_atoms=[Atom(Monomer, element='H', position=None)],
                left_displaced_atoms=[Atom(Backbone, element='O', position=3, charge=-1)]),
            circular=circular)


class CanonicalRnaForm(BpForm):
    """ Canonical RNA form """

    DEFAULT_FASTA_CODE = 'N'

    def __init__(self, seq=None, circular=False):
        """
        Args:
            seq (:obj:`MonomerSequence`, optional): sequence of monomeric forms of the DNA
            circular (:obj:`bool`, optional): if :obj:`True`, indicates that the biopolymer is circular
        """
        super(CanonicalRnaForm, self).__init__(
            seq=seq, alphabet=canonical_rna_alphabet,
            backbone=Backbone(
                structure='OP([O-])([O-])=O',
                monomer_bond_atoms=[Atom(Backbone, element='P', position=2)],
                monomer_displaced_atoms=[Atom(Backbone, element='O', position=1), Atom(Backbone, element='H', position=1)]),
            bond=Bond(
                right_bond_atoms=[Atom(Monomer, element='O', position=None)],
                left_bond_atoms=[Atom(Backbone, element='P', position=2)],
                right_displaced_atoms=[Atom(Monomer, element='H', position=None)],
                left_displaced_atoms=[Atom(Backbone, element='O', position=3, charge=-1)]),
            circular=circular)
