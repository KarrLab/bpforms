""" Alphabet and BpForm to represent modified RNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import (Alphabet, AlphabetBuilder, Monomer, MonomerSequence, Backbone,
                          Bond, Atom, BpForm, Identifier, IdentifierSet, SynonymSet,
                          BpFormsWarning)
from wc_utils.util.chem import EmpiricalFormula
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
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomer
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            path (:obj:`str`, optional): path to save alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        return super(RnaAlphabetBuilder, self).run(ph=ph, major_tautomer=major_tautomer, path=path)

    def build(self, ph=None, major_tautomer=False):
        """ Build alphabet

        Args:
            ph (:obj:`float`, optional): pH at which to calculate the major protonation state of each monomer
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer

        Returns:
            :obj:`Alphabet`: alphabet
        """
        # initialize alphabet
        alphabet = Alphabet()

        # create canonical monomers
        alphabet.from_yaml(canonical_filename)
        alphabet.id = 'rna'
        alphabet.name = 'MODOMICS RNA nucleosides'
        alphabet.description = ('The four canonical RNA nucleosides, plus the non-canonical RNA nucleosides in '
                                '<a href="http://modomics.genesilico.pl/modifications">MODOMICS</a>')

        # create requests session
        cache_name = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'rna'))
        session = requests_cache.core.CachedSession(cache_name, backend='sqlite', expire_after=None)
        session.mount('http://', requests.adapters.HTTPAdapter(max_retries=self.MAX_RETRIES))

        # get originating monomers
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

        # get individual nucleosides and create monomers
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

            if chars in alphabet.monomers:
                warnings.warn('Ignoring canonical monomer {}'.format(chars), BpFormsWarning)
                monomer_short_names[short_name] = alphabet.monomers[chars]
                continue

            monomer = Monomer(
                id=chars,
                name=name,
                synonyms=synonyms,
                identifiers=identifiers,
                structure=structure,
            )

            if not self.is_valid_nucleoside(new_nomenclature, short_name, name, monomer):
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

        # get major microspecies for each monomer
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

    def is_valid_nucleoside(self, new_nomenclature, short_name, name, monomer):
        """ Determine if nucleoside should be included in alphabet

        Args:
            new_nomenclature (:obj:`str`): new nomenclature
            short_name (:obj:`str`): short name
            name (:obj:`str`): name
            monomer (:obj:`Monomer`): monomer

        Returns:
            :obj:`bool`: :obj:`True` if monomer should be included in alphabet
        """
        if ' (base)' in new_nomenclature:
            return False
        if ' (cap)' in name or ' cap' in name:
            return False
        if '-CoA)' in name:
            return False
        if '(p' in short_name:
            return False
        if not monomer.structure:
            return False
        formula = monomer.get_formula()
        if formula.C < 9 or formula.O < 4 or formula.N < 2:
            return False
        return True


class RnaForm(BpForm):
    """ RNA form """

    DEFAULT_FASTA_CODE = 'N'

    def __init__(self, monomer_seq=None, circular=False):
        """
        Args:
            monomer_seq (:obj:`MonomerSequence`, optional): monomers of the DNA form
            circular (:obj:`bool`, optional): if :obj:`True`, indicates that the biopolymer is circular
        """
        super(RnaForm, self).__init__(
            monomer_seq=monomer_seq, alphabet=rna_alphabet,
            backbone=Backbone(
                structure='OP([O-])([O-])=O',
                backbone_bond_atoms=[Atom(Backbone, element='P', position=2)],
                backbone_displaced_atoms=[Atom(Backbone, element='O', position=1), Atom(Backbone, element='H', position=1)]),
            bond=Bond(
                left_bond_atoms=[Atom(Monomer, element='O', position=None)],
                right_bond_atoms=[Atom(Backbone, element='P', position=2)],
                left_displaced_atoms=[Atom(Monomer, element='H', position=None)],
                right_displaced_atoms=[Atom(Backbone, element='O', position=3, charge=-1)]),
            circular=circular)


class CanonicalRnaForm(BpForm):
    """ Canonical RNA form """

    DEFAULT_FASTA_CODE = 'N'

    def __init__(self, monomer_seq=None, circular=False):
        """
        Args:
            monomer_seq (:obj:`MonomerSequence`, optional): monomers of the DNA form
            circular (:obj:`bool`, optional): if :obj:`True`, indicates that the biopolymer is circular
        """
        super(CanonicalRnaForm, self).__init__(
            monomer_seq=monomer_seq, alphabet=canonical_rna_alphabet,
            backbone=Backbone(
                structure='OP([O-])([O-])=O',
                backbone_bond_atoms=[Atom(Backbone, element='P', position=2)],
                backbone_displaced_atoms=[Atom(Backbone, element='O', position=1), Atom(Backbone, element='H', position=1)]),
            bond=Bond(
                left_bond_atoms=[Atom(Monomer, element='O', position=None)],
                right_bond_atoms=[Atom(Backbone, element='P', position=2)],
                left_displaced_atoms=[Atom(Monomer, element='H', position=None)],
                right_displaced_atoms=[Atom(Backbone, element='O', position=3, charge=-1)]),
            circular=circular)
