""" Alphabet and BpForm to represent modified proteins

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Alphabet, AlphabetBuilder, Base, BpForm, Identifier, IdentifierSet, SynonymSet
from wc_utils.util.chem import EmpiricalFormula
from ftplib import FTP
from bs4 import BeautifulSoup
import requests
import requests_cache
import tempfile
import zipfile
import shutil
import os.path
import glob
import pkg_resources
import openbabel
import re
import warnings

filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.yml'))
protein_alphabet = Alphabet().from_yaml(filename)
# :obj:`Alphabet`: Alphabet for protein amino acids

canonical_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.canonical.yml'))
canonical_protein_alphabet = Alphabet().from_yaml(canonical_filename)
# :obj:`Alphabet`: Alphabet for canonical protein amino acids


class ProteinAlphabetBuilder(AlphabetBuilder):
    """ Build protein alphabet from RESID """

    MAX_RETRIES = 5

    def run(self, path=filename):
        """ Build alphabet and, optionally, save to YAML file

        Args:
            path (:obj:`str`, optional): path to save alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        return super(ProteinAlphabetBuilder, self).run(path)

    def build(self):
        """ Build alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        # initialize alphabet
        alphabet = Alphabet()

        # load canonical bases
        alphabet.from_yaml(canonical_filename)
        alphabet.id = 'protein'
        alphabet.name = 'Protein'
        alphabet.description = ('The 20 canonical bases, plus the modified bases in '
                                '<a href="https://pir.georgetown.edu/resid">RESID</a>')

        # get amino acid names from canonical list
        aa_names = []
        for base in alphabet.bases.keys():
            aa_names.append(alphabet.bases[base].name)

        # create tmp directory
        tmp_folder = tempfile.mkdtemp()

        # retrieve files from ftp.ebi.ac.uk and save them to tmpdir
        ftp = FTP('ftp.ebi.ac.uk')
        ftp.login()
        ftp.cwd('pub/databases/RESID/')

        local_filename = os.path.join(tmp_folder, 'models.zip')
        with open(local_filename, 'wb') as file:
            ftp.retrbinary('RETR %s' % 'models.zip', file.write)

        # quit ftp
        ftp.quit()

        # extract pdbs from models.zip into tmp folder
        with zipfile.ZipFile(os.path.join(tmp_folder, 'models.zip'), 'r') as z:
            z.extractall(tmp_folder)

        # create session to get metadata
        cache_name = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein'))
        session = requests_cache.core.CachedSession(cache_name, backend='sqlite', expire_after=None)
        session.mount('http://', requests.adapters.HTTPAdapter(max_retries=self.MAX_RETRIES))

        # extract name of the molecule from pdb file
        for file in glob.iglob(tmp_folder+'/*PDB'):
            with open(file, 'r') as f:
                names = []
                for line in f:
                    if re.match(r"^COMPND    ", line):
                        part1 = str(line[10:].strip())
                        names.append(part1)

                    # check if name is on two lines (when too long)
                    if re.match(r"^COMPND   1", line):
                        part2 = str(line[10:].strip())
                        names.append(part2)
                name = ''.join(names)
                id = re.split("[/.]", file)[3]
                structure = self.get_modification_structure(file)

                chebi_syn = self.get_modification_identifiers_synonyms(id, session)

                synonyms = SynonymSet()
                identifiers = IdentifierSet()

                if chebi_syn is not None:
                    if chebi_syn[2][1:-1]:
                        synonym = chebi_syn[2]
                        for elm in synonym:                     
                            synonyms.add(elm)

                    if chebi_syn[0] != 0:
                        chebi = chebi_syn[0]
                        identifiers.add(Identifier('chebi', 'CHEBI:'+str(chebi)))

                    if chebi_syn[1] != 0:
                        identifiers.add(Identifier('PDBHET', chebi_syn[1]))


            if name in aa_names:
                warnings.warn('Ignoring canonical base {}'.format(name), UserWarning)
                continue

            alphabet.bases[id] = Base(
                id=id,
                name=name,
                synonyms=synonyms,
                identifiers=identifiers,
                structure=structure,
                # comments="Modification of {}.".format(mod['originating_base'])
            )

        # remove tmp folder
        shutil.rmtree(tmp_folder)


        return alphabet

    def get_modification_structure(self, input_pdb):
        """ Get the structure of a modified AA from pdb structure

        Args:
            input_pdb (:obj:`str`): pdb structure

        Returns:
            :obj:`openbabel.OBMol`: structure
        """

        # get InChI from pdb structure
        mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        handler = openbabel.OBMessageHandler()
        handler.SetOutputLevel(2)
        conv.ReadFile(mol, input_pdb)
        conv.SetOutFormat('inchi')

        return mol

    def get_modification_identifiers_synonyms(self, id, session):
        """ Get the chebi ID and synonyms of a modified AA from pdb structure

        Args:
            input_pdb (:obj:`str`): id of RESID entry

        Returns:
            :obj: `list of :obj:`str`: list of synonyms
            :obj: `int`: ChEBI id 
            :obj: `str`: PDB HETATM name
        """

        page = session.get('https://pir.georgetown.edu/cgi-bin/resid?id='+id)
        soup = BeautifulSoup(page.text, features="lxml")

        element = soup.select('p.annot')
        names = element[0].get_text()
        cross_references = element[1].get_text()

        # get synonyms
        if 'Alternate names' in names:
            l = re.split("[:;]", names.strip())[1:]
            synonyms = list(map(lambda x: x.strip(), l))

        # get ChEBI id and HETATM name
        i_chebi = 0
        pdbhet = 0
        if 'Cross-references' in cross_references:
            l = re.split("[:;]", cross_references.strip())[1:]
            l2 = list(map(lambda x: x.strip(), l))
            if 'ChEBI' in l2:
                i_chebi = int(l2[l2.index('ChEBI')+1])
            else:
                i_chebi = 0

            if 'PDBHET' in l2:
                pdbhet = l2[l2.index('PDBHET')+1]
            else:
                pdbhet = 0

            return i_chebi, pdbhet, synonyms


class ProteinForm(BpForm):
    """ Protein form """

    def __init__(self, base_seq=None):
        """
        Args:
            base_seq (:obj:`BaseSequence`, optional): bases of the protein form
        """
        super(ProteinForm, self).__init__(base_seq=base_seq, alphabet=protein_alphabet,
                                          bond_charge=0)


class CanonicalProteinForm(BpForm):
    """ Canonical protein form """

    def __init__(self, base_seq=None):
        """
        Args:
            base_seq (:obj:`BaseSequence`, optional): bases of the protein form
        """
        super(CanonicalProteinForm, self).__init__(base_seq=base_seq, alphabet=canonical_protein_alphabet,
                                                   bond_formula=EmpiricalFormula('H2O') * -1, bond_charge=0)
