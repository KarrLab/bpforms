""" Alphabet and BpForm to represent modified proteins

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import (Alphabet, AlphabetBuilder, Monomer, MonomerSequence, Backbone,
                          Bond, Atom, BpForm, Identifier, IdentifierSet, SynonymSet)
from bs4 import BeautifulSoup
from ftplib import FTP
from wc_utils.util.chem import EmpiricalFormula, OpenBabelUtils, get_major_micro_species
import glob
import openbabel
import os.path
import pkg_resources
import re
import requests
import requests_cache
import warnings
import zipfile

filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.yml'))
protein_alphabet = Alphabet().from_yaml(filename)
# :obj:`Alphabet`: Alphabet for protein residues

canonical_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.canonical.yml'))
canonical_protein_alphabet = Alphabet().from_yaml(canonical_filename)
# :obj:`Alphabet`: Alphabet for canonical protein residues

# multiple_peptide_bonds_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.multiple-peptide-bonds.yml'))
# multiple_peptide_bonds_protein_alphabet = Alphabet().from_yaml(multiple_peptide_bonds_filename)
# :obj:`Alphabet`: Alphabet for protein residues that contain multiple peptide bonds


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
        # return super(ProteinAlphabetBuilder, self).run(ph=ph, major_tautomer=major_tautomer, path=path)

        alphabet = self.build()
        # if ph is not None:
        #    alphabet.get_major_micro_species(ph, major_tautomer=major_tautomer)
        if path:
            self.save(alphabet, path)
        return alphabet

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
        alphabet.name = 'RESID protein residues'
        alphabet.description = ('The 20 canonical protein residues, plus the non-canonical protein residues in '
                                '<a href="https://pir.georgetown.edu/resid">RESID</a>')

        # get amino acid names from canonical list
        canonical_aas = {}
        for monomer in alphabet.monomers.values():
            canonical_aas[monomer.name] = monomer

        # create directory
        pdb_dir = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.pdb'))
        if not os.path.isdir(pdb_dir):
            os.mkdir(pdb_dir)

        # download and unzip PDB files
        zip_filename = os.path.join(pdb_dir, 'models.zip')
        if not os.path.isfile(zip_filename):
            # retrieve files from ftp.ebi.ac.uk and save them to tmpdir
            ftp = FTP('ftp.ebi.ac.uk')
            ftp.login()
            ftp.cwd('pub/databases/RESID/')

            with open(zip_filename, 'wb') as file:
                ftp.retrbinary('RETR %s' % 'models.zip', file.write)

            # quit ftp
            ftp.quit()

            # extract pdbs from models.zip into tmp folder
            with zipfile.ZipFile(os.path.join(pdb_dir, 'models.zip'), 'r') as z:
                z.extractall(pdb_dir)

        # create session to get metadata
        cache_name = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein'))
        session = requests_cache.core.CachedSession(cache_name, backend='sqlite', expire_after=None)
        session.mount('http://', requests.adapters.HTTPAdapter(max_retries=self.MAX_RETRIES))

        # extract name of the molecule from pdb file
        base_monomers = {}
        monomer_ids = {}
        for file in glob.iglob(pdb_dir + '/*.PDB'):
            id = re.split("[/.]", file)[3]
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

            result = self.get_monomer_structure(name, file)
            if result is None:
                warnings.warn('Ignoring monomer {} that has no structure'.format(id), UserWarning)
                continue
            structure, index_n, index_c = result

            code, synonyms, identifiers, base_monomer_ids, comments = self.get_monomer_details(id, session)
            if code is None or code in alphabet.monomers:
                code = id

            if name in canonical_aas:
                canonical_aa = canonical_aas[name]
                monomer_ids[id] = canonical_aa
                canonical_aa.structure = structure
                canonical_aa.monomer_bond_atoms[0].position = index_c
                canonical_aa.monomer_displaced_atoms[0].position = index_c
                canonical_aa.left_bond_atoms[0].position = index_c
                canonical_aa.right_bond_atoms[0].position = index_n
                canonical_aa.right_displaced_atoms[0].position = index_n
                canonical_aa.right_displaced_atoms[1].position = index_n
                warnings.warn('Updated canonical monomer {}'.format(name), UserWarning)
                continue

            monomer = Monomer(
                id=id,
                name=name,
                synonyms=synonyms,
                identifiers=identifiers,
                structure=structure,
                comments=comments,
                monomer_bond_atoms=[Atom(Monomer, element='C', position=index_c)],
                monomer_displaced_atoms=[Atom(Monomer, element='H', position=index_c)],
                left_bond_atoms=[Atom(Monomer, element='C', position=index_c)],
                right_bond_atoms=[Atom(Monomer, element='N', position=index_n, charge=-1)],
                right_displaced_atoms=[Atom(Monomer, element='H', position=index_n, charge=1),
                                       Atom(Monomer, element='H', position=index_n)],
            )
            alphabet.monomers[code] = monomer

            monomer_ids[id] = monomer
            base_monomers[monomer] = base_monomer_ids

        for monomer, base_monomer_ids in base_monomers.items():
            for base_monomer_id in base_monomer_ids:
                base_monomer = monomer_ids.get(base_monomer_id, None)
                monomer.base_monomers.add(base_monomer)

        return alphabet

    def get_monomer_structure(self, name, pdb_filename):
        """ Get the structure of an amino acid from a PDB file
        where N from NH moeity and C from CO moeity are labeled as N15 and C13 isotopes

        Args:
            name (:obj:`str`): monomer name
            pdb_filename (:obj:`str`): path to PDB file with structure

        Returns:
            :obj:`openbabel.OBMol`: structure
        """
        pdb_mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        assert conv.SetInFormat('pdb'), 'Unable to set format to PDB'
        conv.ReadFile(pdb_mol, pdb_filename)

        # removing modified monomers where metal present in structure because:
        # - inchi structure generated separates each non covalently bound parts of the monomer
        # - for many cases theses structures consist of a group of modified monomers coordinating
        # a metal, and not a single PTM monomer per se
        assert conv.SetOutFormat('inchi'), 'Unable to set format to InChI'
        inchi = conv.WriteString(pdb_mol)
        formula = inchi.split('/')[1]
        if '.' in formula:
            warnings.warn('Ignoring metal coordinated monomer {}'.format(name), UserWarning)
            return None

        assert conv.SetOutFormat('smi'), 'Unable to set format to SMILES'
        conv.SetOptions('c', conv.OUTOPTIONS)
        smiles = conv.WriteString(pdb_mol).partition('\t')[0]
        smiles = get_major_micro_species(smiles, 'smiles', 7.4, major_tautomer=True)
        smiles_mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        assert conv.SetInFormat('smi'), 'Unable to set format to SMILES'
        conv.ReadString(smiles_mol, smiles)
        smiles = conv.WriteString(smiles_mol)
        conv.ReadString(smiles_mol, smiles)

        # count the total number of atoms in molecule and loop over each atom
        atomcount = pdb_mol.NumAtoms()
        res = pdb_mol.GetResidue(0)
        countN = 0
        countC = 0

        for i in range(1, atomcount + 1):

            # since N-HN and C-O of peptide bonds must be consecutive in pdb file, get atom_i and atom_i+1
            if pdb_mol.GetAtom(i + 1):
                atom1 = pdb_mol.GetAtom(i)
                atom2 = pdb_mol.GetAtom(i + 1)

                # exception for Proline based-residue where no HN is present
                if res.GetName() == 'Pro':
                    if res.GetNumAtoms() == 1:
                        for res in openbabel.OBResidueIter(pdb_mol):
                            for atom in openbabel.OBResidueAtomIter(res):
                                if atom1.GetType() == 'N3' and atom.GetIdx() == atom1.GetIdx():
                                    index_n = atom1.GetIdx()
                                    atom1.SetIsotope(15)
                                if atom1.GetType() == 'C2' and atom2.GetType() == 'O2' and atom.GetIdx() == atom1.GetIdx():
                                    index_c = atom1.GetIdx()
                                    atom1.SetIsotope(13)
                    else:
                        if atom1.GetType() == 'N3' and (res.GetAtomID(atom1)).strip() == 'N':
                            index_n = atom1.GetIdx()
                            atom1.SetIsotope(15)
                        if atom1.GetType() == 'C2' and atom2.GetType() == 'O2' and (res.GetAtomID(atom1)).strip() == 'C':
                            index_c = atom1.GetIdx()
                            atom1.SetIsotope(13)

                # need to check first if residue numbering is present in a correct way,
                # some residues do not have a proper residue number and name but
                # just a duplicate of the atom number
                if res.GetNumAtoms() == 1:
                    for res_i in openbabel.OBResidueIter(pdb_mol):
                        for atom in openbabel.OBResidueAtomIter(res_i):
                            if atom1.GetType() == 'N3' and atom2.GetType() == 'H' and atom.GetIdx() == atom1.GetIdx() and countN == 0:
                                if pdb_mol.GetAtom(i+2) and (pdb_mol.GetAtom(i+2)).GetType() == 'C3':
                                    countN = 1
                                    index_n = atom1.GetIdx()
                                    atom1.SetIsotope(15)
                            if atom1.GetType() == 'C2' and atom2.GetType() == 'O2' and atom.GetIdx() == atom1.GetIdx() and countC == 0:
                                countC = 1
                                index_c = atom1.GetIdx()
                                atom1.SetIsotope(13)

                if res.GetName() != 'Pro' and res.GetNumAtoms() != 1:
                    # check types according to openbabel atom types for N and HN
                    if atom1.GetType() == 'N3' and atom2.GetType() == 'H' and (res.GetAtomID(atom1)).strip() == 'N' and countN == 0:
                        if pdb_mol.GetAtom(i+2) and (pdb_mol.GetAtom(i+2)).GetType() == 'C3':
                            countN = 1
                            index_n = atom1.GetIdx()
                            atom1.SetIsotope(15)
                    # check types according to openbabel atom types for C and O
                    if atom1.GetType() == 'C2' \
                            and atom2.GetType() == 'O2' \
                            and (res.GetAtomID(atom1)).strip() == 'C' \
                            and (res.GetAtomID(atom2)).strip() == 'O' \
                            and countC == 0:
                        countC = 1
                        index_c = atom1.GetIdx()
                        atom1.SetIsotope(13)

        assert conv.SetOutFormat('inchi'), 'Unable to set format to InChI'
        inchi_isotopes = conv.WriteString(pdb_mol)

        # removing modified monomers where metal present in structure because:
        # - inchi structure generated separates each non covalently bound parts of the monomer
        # - for many cases theses structures consist of a group of modified monomers coordinating
        # a metal, and not a single PTM monomer per se
        formula = inchi_isotopes.split('/')[1]
        if '.' in formula:
            warnings.warn('Ignoring metal coordinated monomer {}'.format(name), UserWarning)
            return None

        # create molecule from SMILES -- necessary to sanitize molecule from PDB
        assert conv.SetOutFormat('smi'), 'Unable to set format to SMILES'
        conv.SetOptions('c', conv.OUTOPTIONS)
        smiles_isotopes = conv.WriteString(pdb_mol).partition('\t')[0]
        smiles_isotopes = get_major_micro_species(smiles_isotopes, 'smiles', 7.4, major_tautomer=True)
        smiles_mol_isotopes = openbabel.OBMol()
        assert conv.SetInFormat('smi'), 'Unable to set format to SMILES'
        conv.ReadString(smiles_mol_isotopes, smiles_isotopes)
        smiles_isotopes = conv.WriteString(smiles_mol_isotopes)
        conv.ReadString(smiles_mol_isotopes, smiles_isotopes)

        i_n_isotopes = None
        i_c_isotopes = None
        i_heavy = 0
        for i_atom in range(1, smiles_mol_isotopes.NumAtoms() + 1):
            atom = smiles_mol_isotopes.GetAtom(i_atom)
            if atom.GetAtomicNum() > 1:
                i_heavy += 1
            if atom.GetIsotope() == 15:
                i_n_isotopes = i_heavy
            if atom.GetIsotope() == 13:
                i_c_isotopes = i_heavy

        i_n = None
        i_c = None
        i_heavy = 0
        for i_atom in range(1, smiles_mol.NumAtoms() + 1):
            atom = smiles_mol.GetAtom(i_atom)
            if atom.GetAtomicNum() > 1:
                i_heavy += 1
                if i_heavy == i_n_isotopes:
                    i_n = i_atom
                    assert atom.GetAtomicNum() == 7
                if i_heavy == i_c_isotopes:
                    i_c = i_atom
                    assert atom.GetAtomicNum() == 6

        if None in (i_n, i_c):
            warnings.warn('Ignoring monomer {} without bonding possibility'.format(name), UserWarning)
            return None

        return (smiles_mol, i_n, i_c)

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

        page = session.get('https://proteininformationresource.org/cgi-bin/resid?id='+id)
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

    DEFAULT_FASTA_CODE = 'X'

    def __init__(self, monomer_seq=None, circular=False):
        """
        Args:
            monomer_seq (:obj:`MonomerSequence`, optional): monomers of the protein form
            circular (:obj:`bool`, optional): if :obj:`True`, indicates that the biopolymer is circular
        """
        super(ProteinForm, self).__init__(
            monomer_seq=monomer_seq, alphabet=protein_alphabet,
            backbone=Backbone(
                structure='[OH-]',
                backbone_bond_atoms=[Atom(Backbone, element='O', position=1)],
                backbone_displaced_atoms=[Atom(Backbone, element='H', position=1)]),
            bond=Bond(
                left_bond_atoms=[Atom(Monomer, element='C', position=None)],
                right_bond_atoms=[Atom(Monomer, element='N', position=None)],
                left_displaced_atoms=[Atom(Backbone, element='O', position=1)],
                right_displaced_atoms=[Atom(Monomer, element='H', position=None), Atom(Monomer, element='H', position=None)]),
            circular=circular)


class CanonicalProteinForm(BpForm):
    """ Canonical protein form """

    DEFAULT_FASTA_CODE = 'X'

    def __init__(self, monomer_seq=None, circular=False):
        """
        Args:
            monomer_seq (:obj:`MonomerSequence`, optional): monomers of the protein form
            circular (:obj:`bool`, optional): if :obj:`True`, indicates that the biopolymer is circular
        """
        super(CanonicalProteinForm, self).__init__(
            monomer_seq=monomer_seq, alphabet=canonical_protein_alphabet,
            backbone=Backbone(
                structure='[OH-]',
                backbone_bond_atoms=[Atom(Backbone, element='O', position=1)],
                backbone_displaced_atoms=[Atom(Backbone, element='H', position=1)]),
            bond=Bond(
                left_bond_atoms=[Atom(Monomer, element='C', position=None)],
                right_bond_atoms=[Atom(Monomer, element='N', position=None)],
                left_displaced_atoms=[Atom(Backbone, element='O', position=1)],
                right_displaced_atoms=[Atom(Monomer, element='H', position=None), Atom(Monomer, element='H', position=None)]),
            circular=circular)
