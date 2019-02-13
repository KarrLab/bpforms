""" Alphabet and BpForm to represent modified proteins

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Alphabet, AlphabetBuilder, Monomer, MonomerSequence, BpForm, Identifier, IdentifierSet, SynonymSet
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

        # load canonical monomers
        alphabet.from_yaml(canonical_filename)
        alphabet.id = 'protein'
        alphabet.name = 'RESID protein amino acids'
        alphabet.description = ('The 20 canonical amino acids, plus the non-canonical amino acids in '
                                '<a href="https://pir.georgetown.edu/resid">RESID</a>')

        # get amino acid names from canonical list
        canonical_aas = {}
        for monomer in alphabet.monomers.keys():
            canonical_aas[alphabet.monomers[monomer].name] = alphabet.monomers[monomer]

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
                structure = self.get_monomer_structure(file)

                chebi_syn = self.get_monomer_identifiers_synonyms(id, session)

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

                # removing modified monomers where metal present in structure because:
                # - inchi structure generated separates each non covalently bound parts of the monomer
                # - for many cases theses structures consist of a group of modified monomers coordinating
                # a metal, and not a single PTM monomer per se
                inchi_formula = self.get_monomer_inchi_formula(structure)
                if '.' in inchi_formula:
                    warnings.warn('Ignoring metal coordinated monomer {}'.format(name), UserWarning)
                    continue                

                if name in canonical_aas:
                    canonical_aas[name].structure = structure
                    warnings.warn('Ignoring canonical monomer {}'.format(name), UserWarning)
                    continue

                monomer = Monomer(
                    id=id,
                    name=name,
                    synonyms=synonyms,
                    identifiers=identifiers,
                    structure=structure,
                )
                alphabet.monomers[id] = monomer

        # remove tmp folder
        shutil.rmtree(tmp_folder)

        return alphabet

    def get_monomer_structure(self, pdb_filename):
        """ Get the structure of an amino acid from a PDB file

        Args:
            pdb_filename (:obj:`str`): path to PDB file with structure

        Returns:
            :obj:`openbabel.OBMol`: structure
        """

        # get InChI from pdb structure
        pdb_mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        assert conv.SetInFormat('pdb'), 'Unable to set format to PDB'
        conv.ReadFile(pdb_mol, pdb_filename)
        assert conv.SetOutFormat('inchi'), 'Unable to set format to InChI'
        inchi = conv.WriteString(pdb_mol)

        # create molecule from InChI -- necessary to sanitize molecule from PDB
        inchi_mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        assert conv.SetInFormat('inchi'), 'Unable to set format to InChI'
        conv.ReadString(inchi_mol, inchi)

        return inchi_mol

    def get_monomer_identifiers_synonyms(self, id, session):
        """ Get the CHEBI ID and synonyms of an amino acid from its RESID webpage

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

    def get_monomer_inchi_formula(self, mol):
        """ Get the inchi formula of a modified amino acid

        Args:
            mol (:obj:`openbabel.OBMol`): inchi structure

        Returns:
            :obj:`int`: inchi chemical formula
        """

        conv = openbabel.OBConversion()
        conv.SetOutFormat('inchi')
        inchi = conv.WriteString(mol).strip()
        inchi_formula = re.split('[/]', inchi)[1]

        return inchi_formula

class ProteinForm(BpForm):
    """ Protein form """

    def __init__(self, monomer_seq=None):
        """
        Args:
            monomer_seq (:obj:`MonomerSequence`, optional): monomers of the protein form
        """
        super(ProteinForm, self).__init__(monomer_seq=monomer_seq, alphabet=protein_alphabet,
                                          backbone_formula=EmpiricalFormula('O'), backbone_charge=0,
                                          bond_formula=EmpiricalFormula('H2O') * -1, bond_charge=0)


class CanonicalProteinForm(BpForm):
    """ Canonical protein form """

    def __init__(self, monomer_seq=None):
        """
        Args:
            monomer_seq (:obj:`MonomerSequence`, optional): monomers of the protein form
        """
        super(CanonicalProteinForm, self).__init__(monomer_seq=monomer_seq, alphabet=canonical_protein_alphabet,
                                                   backbone_formula=EmpiricalFormula('O'), backbone_charge=0,
                                                   bond_formula=EmpiricalFormula('H2O') * -1, bond_charge=0)
