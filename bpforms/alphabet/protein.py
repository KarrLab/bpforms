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
from xml.etree.ElementTree import ElementTree
import glob
import jnius
import openbabel
import os.path
import pkg_resources
import re
import requests
import requests_cache
import tarfile
import warnings
import zipfile

filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.yml'))
protein_alphabet = Alphabet().from_yaml(filename)
# :obj:`Alphabet`: Alphabet for protein residues

canonical_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.canonical.yml'))
canonical_protein_alphabet = Alphabet().from_yaml(canonical_filename)
# :obj:`Alphabet`: Alphabet for canonical protein residues


class ProteinAlphabetBuilder(AlphabetBuilder):
    """ Build protein alphabet from `RESID <https://proteininformationresource.org/resid/>`_ """

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
        return super(ProteinAlphabetBuilder, self).run(ph=ph, major_tautomer=major_tautomer, path=path)

    def build(self, ph=None, major_tautomer=False):
        """ Build alphabet

        Args:
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer

        Returns:
            :obj:`Alphabet`: alphabet
        """
        # initialize alphabet
        alphabet = Alphabet()

        # load canonical monomeric forms
        alphabet.from_yaml(canonical_filename)
        alphabet.id = 'protein'
        alphabet.name = 'protein residues'
        alphabet.description = ('The 20 canonical protein residues, plus the non-canonical protein residues in the '
                                '<a href="http://www.wwpdb.org/data/ccd">PDB Chemical Component Dictionary</a> '
                                'and <a href="https://pir.georgetown.edu/resid">RESID</a>')

        # build from sources
        self.build_from_resid(alphabet, ph=ph, major_tautomer=major_tautomer)
        self.build_from_pdb(alphabet, ph=ph, major_tautomer=major_tautomer)

        # return alphabet
        return alphabet

    def build_from_resid(self, alphabet, ph=None, major_tautomer=False):
        """ Build alphabet from RESID

        Args:
            alphabet (:obj:`Alphabet`): alphabet
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
        """
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

            result = self.get_resid_monomer_structure(name, file, ph=ph, major_tautomer=major_tautomer)
            if result is None:
                warnings.warn('Ignoring monomeric form {} that has no structure'.format(id), BpFormsWarning)
                continue
            structure, index_n, index_c = result

            code, synonyms, identifiers, base_monomer_ids, comments = self.get_resid_monomer_details(id, session)
            if code is None or code in alphabet.monomers:
                code = id

            if name in canonical_aas:
                canonical_aa = canonical_aas[name]
                monomer_ids[id] = canonical_aa

                if name not in ['L-arginine', 'L-asparagine', 'L-valine']:
                    canonical_aa.structure = structure
                    canonical_aa.backbone_bond_atoms[0].position = index_c
                    canonical_aa.backbone_displaced_atoms[0].position = index_c
                    canonical_aa.right_bond_atoms[0].position = index_c
                    canonical_aa.left_bond_atoms[0].position = index_n
                    canonical_aa.left_displaced_atoms[0].position = index_n
                    canonical_aa.left_displaced_atoms[1].position = index_n
                    warnings.warn('Updated canonical monomeric form {}'.format(name), BpFormsWarning)

                continue

            monomer = Monomer(
                id=id,
                name=name,
                synonyms=synonyms,
                identifiers=identifiers,
                structure=structure,
                comments=comments,
            )
            if index_c:
                monomer.backbone_bond_atoms.append(Atom(Monomer, element='C', position=index_c))
                monomer.backbone_displaced_atoms.append(Atom(Monomer, element='H', position=index_c))
                monomer.right_bond_atoms.append(Atom(Monomer, element='C', position=index_c))
            if index_n:
                monomer.left_bond_atoms.append(Atom(Monomer, element='N', position=index_n, charge=-1))
                monomer.left_displaced_atoms.append(Atom(Monomer, element='H', position=index_n, charge=1))
                monomer.left_displaced_atoms.append(Atom(Monomer, element='H', position=index_n))

            alphabet.monomers[code] = monomer

            monomer_ids[id] = monomer
            base_monomers[monomer] = base_monomer_ids

        for monomer, base_monomer_ids in base_monomers.items():
            for base_monomer_id in base_monomer_ids:
                base_monomer = monomer_ids.get(base_monomer_id, None)
                if base_monomer == None:
                    warnings.warn('Base {} for {} is invalid'.format(base_monomer_id, monomer.id),
                                  BpFormsWarning)  # pragma: no cover # all parent entries are retained
                else:
                    monomer.base_monomers.add(base_monomer)

    def get_resid_monomer_structure(self, name, pdb_filename, ph=None, major_tautomer=False):
        """ Get the structure of an amino acid from a PDB file

        Args:
            name (:obj:`str`): name of monomeric form
            pdb_filename (:obj:`str`): path to PDB file with structure
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer

        Returns:
            :obj:`openbabel.OBMol`: structure
            :obj:`int`: index of atom of N terminus
            :obj:`int`: index of atom of C terminus
        """
        pdb_mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        assert conv.SetInFormat('pdb'), 'Unable to set format to PDB'
        conv.ReadFile(pdb_mol, pdb_filename)

        # removing modified monomeric forms where metal present in structure because:
        # - inchi structure generated separates each non covalently bound parts of the monomeric form
        # - for many cases theses structures consist of a group of modified monomeric form coordinating
        # a metal, and not a single PTM monomeric form per se
        assert conv.SetOutFormat('inchi'), 'Unable to set format to InChI'
        inchi = conv.WriteString(pdb_mol)
        formula = inchi.split('/')[1]
        if '.' in formula:
            warnings.warn('Ignoring metal coordinated monomeric form {}'.format(name), BpFormsWarning)
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

        return self.get_termini(smiles_mol)

    def build_from_pdb(self, alphabet, ph=None, major_tautomer=False):
        """ Build alphabet from `PDB Chemical Component Dictionary <http://www.wwpdb.org/data/ccd>`_

        Args:
            alphabet (:obj:`Alphabet`): alphabet
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
        """
        filename = pkg_resources.resource_filename('bpforms',
                                                   os.path.join('alphabet', 'PDB', 'components-pub-xml.tar.gz'))
        if not os.path.isfile(filename):
            response = requests.get('http://ligand-expo.rcsb.org/dictionaries/components-pub-xml.tar.gz')
            response.raise_for_status()
            with open(filename, 'wb') as file:
                file.write(response.content)

        smiles_to_monomer = {}
        for monomer in alphabet.monomers.values():
            smiles_to_monomer[monomer.export('smiles')] = monomer

        conv = openbabel.OBConversion()
        assert conv.SetInFormat('smi'), 'Unable to set format to SMILES'
        assert conv.SetOutFormat('smi'), 'Unable to set format to SMILES'
        conv.SetOptions('c', conv.OUTOPTIONS)

        ns = '{http://pdbml.pdb.org/schema/pdbx-v40.xsd}'
        with tarfile.open(filename, 'r:gz') as tar_file:
            i_file = 0
            n_files = len(tar_file.getmembers())
            n_monomers = 0
            for file_info in tar_file:
                i_file += 1
                if i_file % 100 == 1:
                    print('Processing file {} of {}'.format(i_file, n_files))

                if os.path.splitext(file_info.name)[-1] != '.xml':
                    continue

                xml_file = tar_file.extractfile(file_info)
                xml_root = ElementTree().parse(xml_file)
                xml_group = xml_root.find(ns + 'chem_compCategory')
                if xml_group is None:
                    continue  # pragma: no cover # element is always present

                xml_comp = xml_group.find(ns + 'chem_comp')
                if xml_comp is None:
                    continue  # pragma: no cover # element is always present
                id = xml_comp.get('id')
                identifiers = IdentifierSet([Identifier('pdb-ccd', id)])

                xml_el = xml_comp.find(ns + 'pdbx_release_status')
                if xml_el is None or xml_el.text != 'REL':
                    continue

                xml_el = xml_comp.find(ns + 'type')
                if xml_el is None or xml_el.text not in ['L-peptide linking',
                                                         'L-peptide COOH carboxy terminus',
                                                         'L-peptide NH3 amino terminus']:
                    continue

                xml_el = xml_comp.find(ns + 'pdbx_ambiguous_flag')
                if xml_el is None or xml_el.text != 'N':
                    continue  # pragma: no cover # element is always present

                xml_el = xml_comp.find(ns + 'name')
                if xml_el is None:
                    name = None  # pragma: no cover # element is always present
                else:
                    name = xml_el.text.lower()

                synonyms = SynonymSet()
                xml_el = xml_comp.find(ns + 'one_letter_code')
                if xml_el is not None:
                    synonyms.add(xml_el.text)

                for xml_group in xml_root.findall(ns + 'pdbx_chem_comp_identifierCategory'):
                    for xml_subgroup in xml_group.findall(ns + 'pdbx_chem_comp_identifier'):
                        if xml_subgroup.get('type') == 'SYSTEMATIC NAME':
                            synonyms.add(xml_subgroup.find(ns + 'identifier').text)

                smiles = None
                for xml_group in xml_root.findall(ns + 'pdbx_chem_comp_descriptorCategory'):
                    for xml_subgroup in xml_group.findall(ns + 'pdbx_chem_comp_descriptor'):
                        if xml_subgroup.get('type') == 'SMILES_CANONICAL' and \
                                xml_subgroup.get('program') == 'OpenEye OEToolkits':
                            smiles = xml_subgroup.find(ns + 'descriptor').text
                if smiles is None:
                    continue  # pragma: no cover # element is always present

                if ph:
                    try:
                        smiles = get_major_micro_species(smiles, 'smiles', 'smiles', ph, major_tautomer=major_tautomer)
                    except jnius.JavaException:
                        continue

                mol = openbabel.OBMol()
                conv.ReadString(mol, smiles)
                smiles = conv.WriteString(mol).partition('\t')[0]

                mol = openbabel.OBMol()
                conv.ReadString(mol, smiles)
                smiles = conv.WriteString(mol).partition('\t')[0]

                monomer = smiles_to_monomer.get(smiles, None)
                n_monomers += 1
                if monomer is not None:
                    monomer.synonyms.add(name)
                    monomer.synonyms.update(synonyms)
                    monomer.identifiers.update(identifiers)
                else:
                    monomer = Monomer(id=id, name=name, synonyms=synonyms,
                                      identifiers=identifiers, structure=smiles)

                    _, i_n, i_c = self.get_termini(monomer.structure)
                    if i_c:
                        monomer.backbone_bond_atoms.append(Atom(Monomer, element='C', position=i_c))
                        monomer.backbone_displaced_atoms.append(Atom(Monomer, element='H', position=i_c))
                        monomer.right_bond_atoms.append(Atom(Monomer, element='C', position=i_c))
                    if i_n:
                        monomer.left_bond_atoms.append(Atom(Monomer, element='N', position=i_n, charge=-1))
                        monomer.left_displaced_atoms.append(Atom(Monomer, element='H', position=i_n, charge=1))
                        monomer.left_displaced_atoms.append(Atom(Monomer, element='H', position=i_n))

                    assert id not in alphabet.monomers
                    alphabet.monomers[id] = monomer

                if n_monomers == self._max_monomers:
                    break

    def get_termini(self, mol):
        """ Get indices of atoms of N and C termini

        Args:
            mol (:obj:`openbabel.OBMol`): molecule

        Returns:
            :obj:`openbabel.OBMol`: structure
            :obj:`int`: index of atom of N terminus
            :obj:`int`: index of atom of C terminus
        """
        # find N and C termini
        atom_ns = []
        atom_cs = []
        for i_atom in range(1, mol.NumAtoms() + 1):
            atom = mol.GetAtom(i_atom)
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
            if atom_cs and not atom_ns:
                termini.append((None, atom_cs[0]))
            if atom_ns and not atom_cs:
                termini.append((atom_ns[0], None))

        if termini:
            atom_n, atom_c = termini[0]

            if atom_n:
                idx_n = atom_n.GetIdx()
            else:
                idx_n = None

            if atom_c:
                idx_c = atom_c.GetIdx()
            else:
                idx_c = None
        else:
            idx_n = None
            idx_c = None

        return (mol, idx_n, idx_c)

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

    def get_resid_monomer_details(self, id, session):
        """ Get the CHEBI ID and synonyms of an amino acid from its RESID webpage

        Args:
            input_pdb (:obj:`str`): id of RESID entry

        Returns:
            :obj:`str`: code
            :obj:`SynonymSet`: set of synonyms
            :obj:`IdentifierSet`: set of identifiers
            :obj:`set` of :obj:`str`: ids of base monomeric forms
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

    def __init__(self, seq=None, circular=False):
        """
        Args:
            seq (:obj:`MonomerSequence`, optional): sequence of monomeric forms of the protein
            circular (:obj:`bool`, optional): if :obj:`True`, indicates that the biopolymer is circular
        """
        super(ProteinForm, self).__init__(
            seq=seq, alphabet=protein_alphabet,
            backbone=Backbone(
                structure='[OH-]',
                monomer_bond_atoms=[Atom(Backbone, element='O', position=1)],
                monomer_displaced_atoms=[Atom(Backbone, element='H', position=1)]),
            bond=Bond(
                right_bond_atoms=[Atom(Monomer, element='C', position=None)],
                left_bond_atoms=[Atom(Monomer, element='N', position=None)],
                right_displaced_atoms=[Atom(Backbone, element='O', position=1)],
                left_displaced_atoms=[Atom(Monomer, element='H', position=None), Atom(Monomer, element='H', position=None)]),
            circular=circular)


class CanonicalProteinForm(BpForm):
    """ Canonical protein form """

    DEFAULT_FASTA_CODE = 'X'

    def __init__(self, seq=None, circular=False):
        """
        Args:
            seq (:obj:`MonomerSequence`, optional): sequence of monomeric forms of the protein
            circular (:obj:`bool`, optional): if :obj:`True`, indicates that the biopolymer is circular
        """
        super(CanonicalProteinForm, self).__init__(
            seq=seq, alphabet=canonical_protein_alphabet,
            backbone=Backbone(
                structure='[OH-]',
                monomer_bond_atoms=[Atom(Backbone, element='O', position=1)],
                monomer_displaced_atoms=[Atom(Backbone, element='H', position=1)]),
            bond=Bond(
                right_bond_atoms=[Atom(Monomer, element='C', position=None)],
                left_bond_atoms=[Atom(Monomer, element='N', position=None)],
                right_displaced_atoms=[Atom(Backbone, element='O', position=1)],
                left_displaced_atoms=[Atom(Monomer, element='H', position=None), Atom(Monomer, element='H', position=None)]),
            circular=circular)
