""" Alphabet and BpForm to represent modified proteins

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import (Alphabet, AlphabetBuilder, Monomer, MonomerSequence, Backbone,
                          Bond, Atom, BpForm, Identifier, IdentifierSet, SynonymSet,
                          BpFormsWarning)
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
        return super(ProteinAlphabetBuilder, self).run(ph=ph, major_tautomer=major_tautomer, path=path)

    def build(self, ph=None, major_tautomer=False):
        """ Build alphabet

        Args:
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomer
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer

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
        for file in glob.iglob(os.path.join(pdb_dir, '*.PDB')):
            id, _ = os.path.splitext(os.path.basename(file))
            names = []
            with open(file, 'r') as f:
                for line in f:
                    if re.match(r"^COMPND    ", line):
                        part1 = str(line[10:].strip())
                        names.append(part1)

                    # check if name is on two lines (when too long)
                    if re.match(r"^COMPND   1", line):
                        part2 = str(line[10:].strip())
                        names.append(part2)
            name = ''.join(names)

            result = self.get_monomer_structure(name, file, ph=ph, major_tautomer=major_tautomer)
            if result is None:
                warnings.warn('Ignoring monomer {} that has no structure'.format(id), BpFormsWarning)
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
                warnings.warn('Updated canonical monomer {}'.format(name), BpFormsWarning)
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
                if base_monomer == None:
                    warnings.warn('Base {} for {} is invalid'.format(base_monomer_id, monomer.id), BpFormsWarning)
                else:
                    monomer.base_monomers.add(base_monomer)

        return alphabet

    def get_monomer_structure(self, name, pdb_filename, ph=None, major_tautomer=False):
        """ Get the structure of an amino acid from a PDB file

        Args:
            name (:obj:`str`): monomer name
            pdb_filename (:obj:`str`): path to PDB file with structure
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomer
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer

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
            warnings.warn('Ignoring metal coordinated monomer {}'.format(name), BpFormsWarning)
            return None

        assert conv.SetOutFormat('smi')
        conv.SetOptions('c', conv.OUTOPTIONS)
        smiles = conv.WriteString(pdb_mol).partition('\t')[0]
        if ph is not None:
            smiles = get_major_micro_species(smiles, 'smiles', 'smiles', ph, major_tautomer=major_tautomer)
        smiles_mol = openbabel.OBMol()
        assert conv.SetInFormat('smi')
        conv.ReadString(smiles_mol, smiles)

        assert conv.SetOutFormat('smi'), 'Unable to set format to SMILES'
        conv.SetOptions('c', conv.OUTOPTIONS)
        smiles = conv.WriteString(smiles_mol).partition('\t')[0]
        assert conv.SetInFormat('smi'), 'Unable to set format to SMILES'
        conv.ReadString(smiles_mol, smiles)

        # find N and C termini
        atom_ns = []
        atom_cs = []
        for i_atom in range(1, smiles_mol.NumAtoms() + 1):
            atom = smiles_mol.GetAtom(i_atom)
            if self.is_n_terminus(atom):
                atom_ns.append(atom)
            elif self.is_c_terminus(atom):
                atom_cs.append(atom)

        termini = []
        for atom_n in atom_ns:
            for atom_c in atom_cs:
                if self.is_terminus(atom_n, atom_c):
                    termini.append((atom_n, atom_c))

        if not termini:
            warnings.warn('Ignoring monomer {} without N- and C-termini'.format(name), BpFormsWarning)
            return None
        if len(termini) > 1:
            warnings.warn('Ignoring monomer {} with multiple N- and C-termini'.format(name), BpFormsWarning)
            return None
        atom_n, atom_c = termini[0]

        return (smiles_mol, atom_n.GetIdx(), atom_c.GetIdx())

    def is_n_terminus(self, atom):
        """ Determine if an atom is an N-terminus

        Args:
            atom (:obj:`openbabel.OBAtom`): atom

        Returns:
            :obj:`bool`: :obj:`True` if the atom is an N-terminus
        """
        if atom is None:
            return False

        # check atom is nitrogen
        if atom.GetAtomicNum() != 7:
            return False

        # check atom has charge + 1
        if atom.GetFormalCharge() != 1:
            return False

        # check atom is single-bonded to at least 1 carbon
        has_single_c = False
        for bond in openbabel.OBAtomBondIter(atom):
            other_atom = bond.GetBeginAtom()
            if other_atom == atom:
                other_atom = bond.GetEndAtom()
            if bond.GetBondOrder() == 1 and other_atom.GetAtomicNum() == 6:
                has_single_c = True
        if not has_single_c:
            return False

        # check atom bonded to at least 2 hydrogens
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(atom)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(atom)])
        other_atoms += [1] * (4 - tot_bond_order)
        if other_atoms.count(1) < 2:
            return False

        # return True
        return True

    def is_c_terminus(self, atom):
        """ Determine if an atom is an C-terminus

        Args:
            atom (:obj:`openbabel.OBAtom`): atom

        Returns:
            :obj:`bool`: :obj:`True` if the atom is an C-terminus
        """
        if atom is None:
            return False

        # check atom is hydrogen
        if atom.GetAtomicNum() != 6:
            return False

        # check atom has charge 0
        if atom.GetFormalCharge() != 0:
            return False

        # check atom is bonded to 1 oxygen, 1 carbon, and 1 hydrogen atom
        other_atoms = sorted([other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(atom)])
        if set(other_atoms).difference(set([1, 6, 8])):
            return False
        if 6 not in other_atoms or 8 not in other_atoms:
            return False
        if len(other_atoms) > 2 and other_atoms[-2] != 1:
            return False

        # check bonds to carbon and hydrogen are single bonds and bond to oxygen is a double bond
        for bond in openbabel.OBAtomBondIter(atom):
            other_atom = bond.GetBeginAtom()
            if other_atom == atom:
                other_atom = bond.GetEndAtom()
            if other_atom.GetAtomicNum() != 8 and bond.GetBondOrder() != 1:
                return False
            if other_atom.GetAtomicNum() == 8 and bond.GetBondOrder() != 2:
                return False

        # return True
        return True

    def is_terminus(self, atom_n, atom_c):
        """ Determine if a pair of atoms are N- and C-termini

        Args:
            atom_n (:obj:`openbabel.OBAtom`): potential N-terminus
            atom_c (:obj:`openbabel.OBAtom`): potential C-terminus

        Returns:
            :obj:`bool`: :obj:`True`, if the atoms are N- and C-termini
        """
        atom_n_neighbors = set((atom.GetIdx(), atom.GetAtomicNum()) for atom in openbabel.OBAtomAtomIter(atom_n))
        atom_c_neighbors = set((atom.GetIdx(), atom.GetAtomicNum()) for atom in openbabel.OBAtomAtomIter(atom_c))
        connecting_atoms = atom_n_neighbors.intersection(atom_c_neighbors)

        for idx, atomic_num in connecting_atoms:
            if atomic_num == 6:
                return True

        return False

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
