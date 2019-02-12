""" Alphabet and BpForm to represent modified RNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Alphabet, AlphabetBuilder, Monomer, MonomerSequence, BpForm, Identifier, IdentifierSet, SynonymSet
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
# :obj:`Alphabet`: Alphabet for RNA nucleotides

canonical_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'rna.canonical.yml'))
canonical_rna_alphabet = Alphabet().from_yaml(canonical_filename)
# :obj:`Alphabet`: Alphabet for canonical RNA nucleotides


class RnaAlphabetBuilder(AlphabetBuilder):
    """ Build RNA alphabet from MODOMICS """

    INDEX_ENDPOINT = 'http://modomics.genesilico.pl/modifications/'
    INDEX_ASCII_ENDPOINT = 'http://modomics.genesilico.pl/modifications/?base=all&type=all&display_ascii=Display+as+ASCII'
    ENTRY_ENDPOINT = 'http://modomics.genesilico.pl/modifications/{}/'
    MAX_RETRIES = 5

    def run(self, path=filename):
        """ Build alphabet and, optionally, save to YAML file

        Args:
            path (:obj:`str`, optional): path to save alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        return super(RnaAlphabetBuilder, self).run(path)

    def build(self):
        """ Build alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        # initialize alphabet
        alphabet = Alphabet()

        # create canonical monomers
        alphabet.from_yaml(canonical_filename)
        alphabet.id = 'rna'
        alphabet.name = 'MODOMICS RNA nucleotides'
        alphabet.description = ('The four canonical RNA nucleotides, plus the modified RNA nucleotides in '
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
        orig_monomers = {}
        for monomer in reader:
            for i in range(1, 4):
                short_name = monomer[i]
                if short_name:
                    break
            orig_monomers[short_name] = monomer[i + 2]

        # get index of modifications
        response = session.get(self.INDEX_ENDPOINT)
        response.raise_for_status()

        # get individual modifications and create monomers
        doc = bs4.BeautifulSoup(response.text, 'html.parser')
        table = doc.find('table', {'class': 'datagrid'})
        tbody = table.find('tbody')
        mods = tbody.find_all('tr')
        for i_mod, mod in enumerate(mods):
            if i_mod >= self._max_monomers:
                break

            cells = mod.find_all('td')
            new_nomenclature = cells[0].text
            name = cells[1].text
            short_name = cells[2].text
            abbrev = cells[3].text

            if ' (base)' in new_nomenclature:
                continue

            id = short_name

            chars = new_nomenclature
            if not chars:
                chars = id

            synonyms = SynonymSet()
            if abbrev:
                synonyms.add(abbrev)

            identifiers = IdentifierSet()
            if short_name:
                identifiers.add(Identifier('modomics.short_name', short_name))
            if new_nomenclature:
                identifiers.add(Identifier('modomics.new_nomenclature', new_nomenclature))

            structure = self.get_modification_structure(id, session)

            if chars in alphabet.monomers:
                warnings.warn('Ignoring canonical monomer {}'.format(chars), UserWarning)
                continue

            alphabet.monomers[chars] = Monomer(
                id=chars,
                name=name,
                synonyms=synonyms,
                identifiers=identifiers,
                structure=structure,
                comments="Modification of {}.".format(orig_monomers[short_name])
            )

        # return alphabet
        return alphabet

    def get_modification_structure(self, id, session):
        """ Get the structure of a modified NMP in the MODOMICS database

        Args:
            id (:obj:`str`): id of modification in MODOMICS database

        Returns:
            :obj:`openbabel.OBMol`: structure
        """
        response = session.get(self.ENTRY_ENDPOINT.format(id))
        response.raise_for_status()

        doc = bs4.BeautifulSoup(response.text, 'html.parser')

        table = doc.find(id='modification_details')
        tbody = table.find('tbody')
        tr = tbody.find_all('tr')[-1]
        tds = tr.find_all('td')
        assert tds[0].text.startswith('SMILES'), "Wrong cell retrieved to parse structure"

        td = tds[1]
        smiles = td.text
        if not smiles:
            return None

        mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        assert conv.SetInFormat('smi')
        if not conv.ReadString(mol, smiles):
            return None

        return mol


class RnaForm(BpForm):
    """ RNA form """

    def __init__(self, monomer_seq=None):
        """
        Args:
            monomer_seq (:obj:`MonomerSequence`, optional): monomers of the DNA form
        """
        super(RnaForm, self).__init__(monomer_seq=monomer_seq, alphabet=rna_alphabet,
                                      bond_formula=EmpiricalFormula('OH') * -1, bond_charge=1)


class CanonicalRnaForm(BpForm):
    """ Canonical RNA form """

    def __init__(self, monomer_seq=None):
        """
        Args:
            monomer_seq (:obj:`MonomerSequence`, optional): monomers of the DNA form
        """
        super(CanonicalRnaForm, self).__init__(monomer_seq=monomer_seq, alphabet=canonical_rna_alphabet,
                                               bond_formula=EmpiricalFormula('OH') * -1, bond_charge=1)
