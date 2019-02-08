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
import tempfile
import zipfile
import shutil
import os.path
import glob
import pkg_resources
import openbabel
import re

filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.yml'))
protein_alphabet = Alphabet().from_yaml(filename)
# :obj:`Alphabet`: Alphabet for protein amino acids

canonical_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.canonical.yml'))
canonical_protein_alphabet = Alphabet().from_yaml(canonical_filename)
# :obj:`Alphabet`: Alphabet for canonical protein amino acids

class ProteinAlphabetBuilder(AlphabetBuilder):
    """ Build protein alphabet from RESID """

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

        # create tmp directory
        tmp_folder = tempfile.mkdtemp()

        # retrieve files from ftp.ebi.ac.uk and save them to tmpdir
        ftp = FTP('ftp.ebi.ac.uk')
        ftp.login()
        ftp.cwd('pub/databases/RESID/')
        filenames = ftp.nlst()

        for filename in filenames:
            local_filename = os.path.join(tmp_folder, filename)
            with open(local_filename, 'wb') as file:
               ftp.retrbinary('RETR %s' %filename, file.write)

        # quit ftp
        ftp.quit()
        
        # extract pdbs from models.zip into tmp folder
        with zipfile.ZipFile(os.path.join(tmp_folder,'models.zip'),'r') as z:
            z.extractall(tmp_folder)

        for file in glob.iglob(tmp_folder+'/*PDB'):
            with open(file, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if re.match(r"^COMPND",line) is not None :
                        name = str.split(line)[1]

                id = re.split("[/.]",file)[3]
                # structure = self.get_modification_structure(file, id)

            # problem for now check later
            # if name in alphabet:
                # warnings.warn('Ignoring canonical base {}'.format(name), UserWarning)
                # continue

        # remove tmp folder
        shutil.rmtree(tmp_folder)

        alphabet.bases.A = Base(
            id='ALA',
            name="alanine",
            synonyms=SynonymSet([
                'L-alpha-alanine',
            ]),
            identifiers=IdentifierSet([
                Identifier('pubchem.compound', '7311724'),
                Identifier('metacyc.compound', 'L-ALPHA-ALANINE'),
                Identifier('chebi', 'CHEBI:57972'),
            ]),
            structure='InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1',
        )

        return alphabet

    def get_modification_structure(self, input_pdb, id):
        """ Get the structure of a modified AA from pdb structure

        Args:
            input_pdb (:obj:`str`): pdb structure

        Returns:
            :obj:`openbabel.OBMol`: structure
            :obj: `list of :obj:`str`: list of synonyms
            :obj:`int`: ChEBI id 
        """

        # get InChI from pdb structure
        mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        handler = openbabel.OBMessageHandler()
        handler.SetOutputLevel(0)
        conv.ReadFile(mol, input_pdb)

        # get additional info from page of each modified amino acid
        page = requests.get('https://pir.georgetown.edu/cgi-bin/resid?id='+id)
        soup = BeautifulSoup(page.text, features="lxml")

        element = soup.select('p.annot')
        names = element[0].get_text()
        cross_references = element[1].get_text()

        # get synonyms
        if 'Alternate names' in names:
            l = re.split("[:;]",names.strip())[1:]
            synonyms = list(map(lambda x:x.strip(),l))

        # get ChEBI id
        if 'Cross-references' in cross_references:
            l = re.split("[:;]",cross_references.strip())[1:]
            l2 = list(map(lambda x:x.strip(),l))
            if 'ChEBI' in l2:
                i_chebi = l2[l2.index('ChEBI')+1]
                return mol, synonyms, i_chebi


class ProteinForm(BpForm):
    """ Protein form """

    def __init__(self, base_seq=None):
        """
        Args:
            base_seq (:obj:`BaseSequence`, optional): bases of the DNA form
        """
        super(ProteinForm, self).__init__(base_seq=base_seq, alphabet=protein_alphabet,
                                          bond_formula=EmpiricalFormula('H2O') * -1, bond_charge=0)


class CanonicalProteinForm(BpForm):
    """ Canonical protein form """

    def __init__(self, base_seq=None):
        """
        Args:
            base_seq (:obj:`BaseSequence`, optional): bases of the DNA form
        """
        super(CanonicalProteinForm, self).__init__(base_seq=base_seq, alphabet=canonical_protein_alphabet,
                                                   bond_formula=EmpiricalFormula('H2O') * -1, bond_charge=0)

