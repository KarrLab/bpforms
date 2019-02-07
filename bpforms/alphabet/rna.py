""" Alphabet and BpForm to represent modified RNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Alphabet, AlphabetBuilder, Base, BpForm, Identifier, IdentifierSet, SynonymSet
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

    INDEX_ENDPOINT = 'http://modomics.genesilico.pl/modifications/?base=all&type=all&display_ascii=Display+as+ASCII'
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

        # create canonical bases
        alphabet.from_yaml(canonical_filename)

        # create requests session
        cache_name = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'rna'))
        session = requests_cache.core.CachedSession(cache_name, backend='sqlite', expire_after=None)
        session.mount('http://', requests.adapters.HTTPAdapter(max_retries=self.MAX_RETRIES))

        # get index of modifications
        response = session.get(self.INDEX_ENDPOINT)
        response.raise_for_status()

        # get individual modifications and create bases
        reader = csv.DictReader(io.StringIO(response.text), delimiter='\t', quoting=csv.QUOTE_NONE)
        for mod in reader:
            if ' (base)' in mod['new_nomenclature']:
                continue

            id = mod['short_name']

            chars = mod['new_nomenclature']
            if not chars:
                chars = id

            synonyms = SynonymSet()
            if mod['rnamods_abbrev']:
                synonyms.add(mod['rnamods_abbrev'])

            identifiers = IdentifierSet()
            if mod['short_name']:
                identifiers.add(Identifier('modomics.short_name', mod['short_name']))
            if mod['new_nomenclature']:
                identifiers.add(Identifier('modomics.new_nomenclature', mod['new_nomenclature']))

            structure = self.get_modification_structure(id, session)

            if chars in alphabet.bases:
                warnings.warn('Ignoring canonical base {}'.format(chars), UserWarning)
                continue

            alphabet.bases[chars] = Base(
                id=chars,
                name=mod['name'],
                synonyms=synonyms,
                identifiers=identifiers,
                structure=structure,
                comments="Modification of {}.".format(mod['originating_base'])
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
        if table is None:
            return None

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

    def __init__(self, base_seq=None):
        """
        Args:
            base_seq (:obj:`BaseSequence`, optional): bases of the DNA form
        """
        super(RnaForm, self).__init__(base_seq=base_seq, alphabet=rna_alphabet,
                                      bond_formula=EmpiricalFormula('H') * -1, bond_charge=1)


class CanonicalRnaForm(BpForm):
    """ Canonical RNA form """

    def __init__(self, base_seq=None):
        """
        Args:
            base_seq (:obj:`BaseSequence`, optional): bases of the DNA form
        """
        super(CanonicalRnaForm, self).__init__(base_seq=base_seq, alphabet=canonical_rna_alphabet,
                                               bond_formula=EmpiricalFormula('H') * -1, bond_charge=1)
