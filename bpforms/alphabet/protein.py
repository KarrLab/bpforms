""" Alphabet and BpForm to represent modified proteins

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Alphabet, AlphabetBuilder, Monomer, MonomerSequence, BpForm, Identifier, IdentifierSet, SynonymSet
from bs4 import BeautifulSoup
from ftplib import FTP
from wc_utils.util.chem import EmpiricalFormula
import glob
import openbabel
import os.path
import pkg_resources
import re
import requests
import requests_cache
import shutil
import tempfile
import warnings
import zipfile

filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.yml'))
protein_alphabet = Alphabet().from_yaml(filename)
# :obj:`Alphabet`: Alphabet for protein amino acids

canonical_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.canonical.yml'))
canonical_protein_alphabet = Alphabet().from_yaml(canonical_filename)
# :obj:`Alphabet`: Alphabet for canonical protein amino acids


class ProteinAlphabetBuilder(AlphabetBuilder):
    """ Build protein alphabet from RESID """

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
        return super(ProteinAlphabetBuilder, self).run(ph=ph, major_tautomer=major_tautomer, path=path)

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
        for monomer in alphabet.monomers.values():
            canonical_aas[monomer.name] = monomer

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
        base_monomers = {}
        monomer_ids = {}
        curated_entries = []
        for file in glob.iglob(tmp_folder+'/*PDB'):
            number_hn = 0
            number_co = 0
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

                    # count number of HN and O from peptide bond atom names nomenclature
                    if str.split(line)[0] == 'ATOM' and str.split(line)[2] == 'HN':
                        number_hn += 1
                    if str.split(line)[0] == 'ATOM' and (str.split(line)[2] == 'O'):
                        number_co += 1
                name = ''.join(names)
                id = re.split("[/.]", file)[3]
                structure = self.get_monomer_structure(name, file)
            if not structure:
                continue

            # get the ids of modified monomer with more than one peptide bond
            if number_hn > 1 and number_co > 1:
                curated_entries.append(id)

            code, synonyms, identifiers, base_monomer_ids, comments = self.get_monomer_details(id, session)

            if code is None or code in alphabet.monomers:
                code = id

            if name in canonical_aas:
                monomer_ids[id] = canonical_aas[name]
                canonical_aas[name].structure = structure
                warnings.warn('Ignoring canonical monomer {}'.format(name), UserWarning)
                continue

            monomer = Monomer(
                id=id,
                name=name,
                synonyms=synonyms,
                identifiers=identifiers,
                structure=structure,
                comments=comments,
            )
            alphabet.monomers[code] = monomer


            monomer_ids[id] = monomer
            base_monomers[monomer] = base_monomer_ids

        for monomer, base_monomer_ids in base_monomers.items():
            for base_monomer_id in base_monomer_ids:
                base_monomer = monomer_ids.get(base_monomer_id, None)
                monomer.base_monomers.add(base_monomer)

        # remove tmp folder
        shutil.rmtree(tmp_folder)

        return alphabet

    def get_monomer_structure(self, name, pdb_filename):
        """ Get the structure of an amino acid from a PDB file

        Args:
            name (:obj:`str`): monomer name
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

        # removing modified monomers where metal present in structure because:
        # - inchi structure generated separates each non covalently bound parts of the monomer
        # - for many cases theses structures consist of a group of modified monomers coordinating
        # a metal, and not a single PTM monomer per se
        formula = inchi.split('/')[1]
        if '.' in formula:
            warnings.warn('Ignoring metal coordinated monomer {}'.format(name), UserWarning)
            return None

        # create molecule from InChI -- necessary to sanitize molecule from PDB
        inchi_mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        assert conv.SetInFormat('inchi'), 'Unable to set format to InChI'
        conv.ReadString(inchi_mol, inchi)

        return inchi_mol

    def get_monomer_details(self, id, session):
        """ Get the CHEBI ID and synonyms of an amino acid from its RESID webpage

        Args:
            input_pdb (:obj:`str`): id of RESID entry

        Returns:
            :obj:`str`: code            
            :obj:`SynonymSet`: set of synonyms
            :obj:`IdentifierSet`: set of identifiers
            :obj:`set` of :obj:`str`: ids of base monomers
            :obj:`str`: comments
        """

        page = session.get('https://pir.georgetown.edu/cgi-bin/resid?id='+id)
        soup = BeautifulSoup(page.text, features="lxml")

        paragraphs = soup.select('p.annot')
        code = None
        synonyms = SynonymSet()
        identifiers = IdentifierSet()
        base_monomer_ids = set()
        comments = []
        for paragraph in paragraphs:
            text = paragraph.get_text()

            # code
            if 'Sequence code: ' in text and code is None:
                _, _, code = text.partition('Sequence code: ')
                code, _, _ = code.partition('#')
                code, _, _ = code.partition('\n')
                code = code.strip()

            # get synonyms
            if 'Alternate names:' in text:
                l = re.split("[:;]", text.strip())[1:]
                synonyms.update(map(lambda x: x.strip(), l))

            if 'Systematic name:' in text:
                _, _, systematic_name = text.partition('Systematic name:')
                systematic_name, _, _ = systematic_name.partition('\n')
                synonyms.add(systematic_name.strip())

            # ChEBI id and HETATM name
            if 'Cross-references: ' in text:
                l = re.split("[:;]", text.strip())[1:]
                l2 = list(map(lambda x: x.strip(), l))
                if 'CAS' in l2:
                    identifiers.add(Identifier('cas', l2[l2.index('CAS')+1]))

                if 'ChEBI' in l2:
                    identifiers.add(Identifier('chebi', 'CHEBI:' + l2[l2.index('ChEBI')+1]))

                if 'PDBHET' in l2:
                    identifiers.add(Identifier('pdb.ligand', l2[l2.index('PDBHET')+1]))

                if 'PSI-MOD' in l2:
                    identifiers.add(Identifier('mod', 'MOD:' + l2[l2.index('PSI-MOD')+1]))

                if 'GO' in l2:
                    identifiers.add(Identifier('go', 'GO:' + l2[l2.index('GO')+1]))

            # base amino acid
            if 'Based on ' in text:
                base_monomer_ids.update(text.partition('Based on ')[2].split('+'))

            # comments
            if 'Comment: ' in text:
                comments.append(text.partition('Comment: ')[2].strip())

            if 'Generating Enzyme:' in text:
                _, _, comment = text.partition('Generating Enzyme:')
                comment = comment.strip()
                comment = 'Generating Enzyme: ' + comment
                if comment[-1] != '.':
                    comment += '.'
                comments.append(comment)

        if comments:
            comments = ' '.join(comments)
        else:
            comments = None
        return code, synonyms, identifiers, base_monomer_ids, comments


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
