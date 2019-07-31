""" Alphabet and BpForm to represent modified proteins

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import (Alphabet, AlphabetBuilder, Monomer, MonomerSequence,
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
# :obj:`Alphabet`: Alphabet for protein amino acids

canonical_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.canonical.yml'))
canonical_protein_alphabet = Alphabet().from_yaml(canonical_filename)
# :obj:`Alphabet`: Alphabet for canonical protein amino acids


class ProteinAlphabetBuilder(AlphabetBuilder):
    """ Build protein alphabet from the `PDB Chemical Component Dictionary <http://www.wwpdb.org/data/ccd>` and
    `RESID <https://proteininformationresource.org/resid/>`_ """

    MAX_RETRIES = 5

    def run(self, ph=None, major_tautomer=False, dearomatize=False, path=filename):
        """ Build alphabet and, optionally, save to YAML file

        Args:
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule
            path (:obj:`str`, optional): path to save alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        return super(ProteinAlphabetBuilder, self).run(ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize, path=path)

    def build(self, ph=None, major_tautomer=False, dearomatize=False):
        """ Build alphabet

        Args:
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule

        Returns:
            :obj:`Alphabet`: alphabet
        """
        # initialize alphabet
        alphabet = Alphabet()

        # load canonical monomeric forms
        alphabet.from_yaml(canonical_filename)
        alphabet.id = 'protein'
        alphabet.name = 'Protein amino acids'
        alphabet.description = ('The canonical protein amino acids, plus non-canonical amino acids based on '
                                '<a href="http://www.wwpdb.org/data/ccd">PDB Chemical Component Dictionary</a> '
                                'and <a href="https://pir.georgetown.edu/resid">RESID</a>')

        # build from sources
        self.build_from_resid(alphabet, ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)
        self.build_from_pdb(alphabet, ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)
        self.build_from_mod(alphabet, ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)

        # corrections
        alphabet.monomers.A.identifiers.add(Identifier('chebi', 'CHEBI:46217'))
        alphabet.monomers.C.identifiers.add(Identifier('chebi', 'CHEBI:29950'))
        alphabet.monomers.D.identifiers.add(Identifier('chebi', 'CHEBI:29958'))
        alphabet.monomers.E.identifiers.add(Identifier('chebi', 'CHEBI:29972'))
        alphabet.monomers.F.identifiers.add(Identifier('chebi', 'CHEBI:29997'))
        alphabet.monomers.G.identifiers.add(Identifier('chebi', 'CHEBI:29947'))
        alphabet.monomers.H.identifiers.add(Identifier('chebi', 'CHEBI:29979'))
        alphabet.monomers.I.identifiers.add(Identifier('chebi', 'CHEBI:30009'))
        alphabet.monomers.K.identifiers.add(Identifier('chebi', 'CHEBI:29967'))
        alphabet.monomers.L.identifiers.add(Identifier('chebi', 'CHEBI:30006'))
        alphabet.monomers.M.identifiers.add(Identifier('chebi', 'CHEBI:16044'))
        alphabet.monomers.N.identifiers.add(Identifier('chebi', 'CHEBI:50347'))
        alphabet.monomers.P.identifiers.add(Identifier('chebi', 'CHEBI:50342'))
        alphabet.monomers.Q.identifiers.add(Identifier('chebi', 'CHEBI:30011'))
        alphabet.monomers.R.identifiers.add(Identifier('chebi', 'CHEBI:29952'))
        alphabet.monomers.S.identifiers.add(Identifier('chebi', 'CHEBI:29999'))
        alphabet.monomers.T.identifiers.add(Identifier('chebi', 'CHEBI:30013'))
        alphabet.monomers.V.identifiers.add(Identifier('chebi', 'CHEBI:30015'))
        alphabet.monomers.W.identifiers.add(Identifier('chebi', 'CHEBI:29954'))
        alphabet.monomers.Y.identifiers.add(Identifier('chebi', 'CHEBI:46858'))
        alphabet.monomers.AA0431.identifiers.add(Identifier('mod', 'MOD:00160'))
        alphabet.monomers.AA0560.identifiers.remove(Identifier('pdb.ligand', 'NAG'))
        alphabet.monomers.AA0027.identifiers.remove(Identifier('cas', '17576'))
        alphabet.monomers.AA0230.identifiers.remove(Identifier('pdb.ligand', 'NO'))
        alphabet.monomers.AA0232.identifiers.remove(Identifier('pdb.ligand', 'OTD'))
        alphabet.monomers.D.base_monomers.clear()
        for monomer in alphabet.monomers.values():
            for identifier in list(monomer.identifiers):
                if identifier.ns == 'pdb.ligand' and identifier.id == 'ACE':
                    monomer.identifiers.remove(identifier)
        for monomer in alphabet.monomers.values():
            for base_monomer in list(monomer.base_monomers):
                if base_monomer.id == 'AA0025':
                    monomer.base_monomers.remove(base_monomer)
                    monomer.base_monomers.add(alphabet.monomers.C)

        # save report
        n_only = []
        c_only = []
        no_termini = []
        for code, monomer in alphabet.monomers.items():
            if monomer.l_bond_atoms and not monomer.r_bond_atoms:
                n_only.append((code, monomer.export('smiles', options=('c',))))
            elif not monomer.l_bond_atoms and monomer.r_bond_atoms:
                c_only.append((code, monomer.export('smiles', options=('c',))))
            elif not monomer.l_bond_atoms and not monomer.r_bond_atoms:
                no_termini.append((code, monomer.export('smiles', options=('c',))))
        n_only.sort()
        c_only.sort()
        no_termini.sort()

        filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'protein.report.txt'))
        with open(filename, 'w') as file:
            file.write(('{} monomers only have N termini:\n'
                        '  Code\tSMILES\n'
                        '  {}\n\n').format(
                len(n_only), '\n  '.join('\t'.join(txt) for txt in n_only)))

            file.write(('{} monomers only have C termini:\n'
                        '  Code\tSMILES\n'
                        '  {}\n\n').format(
                len(c_only), '\n  '.join('\t'.join(txt) for txt in c_only)))

            file.write(('{} monomers have no termini:\n'
                        '  Code\tSMILES\n'
                        '  {}\n\n').format(
                len(no_termini), '\n  '.join('\t'.join(txt) for txt in no_termini)))

        # transform from residue to amino acids

        # return alphabet
        return alphabet

    def build_from_resid(self, alphabet, ph=None, major_tautomer=False, dearomatize=False):
        """ Build alphabet from RESID

        Args:
            alphabet (:obj:`Alphabet`): alphabet
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule
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

            result = self.get_resid_monomer_structure(name, file, ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)
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
                    canonical_aa.l_bond_atoms = []
                    canonical_aa.l_displaced_atoms = []
                    canonical_aa.r_bond_atoms = []
                    canonical_aa.r_displaced_atoms = []
                    self.set_termini(structure, canonical_aa, index_n, index_c)

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

            self.set_termini(structure, monomer, index_n, index_c)

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

    def get_resid_monomer_structure(self, name, pdb_filename, ph=None, major_tautomer=False, dearomatize=False):
        """ Get the structure of an amino acid from a PDB file

        Args:
            name (:obj:`str`): name of monomeric form
            pdb_filename (:obj:`str`): path to PDB file with structure
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule

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
            smiles = get_major_micro_species(smiles, 'smiles', 'smiles', ph, major_tautomer=major_tautomer, dearomatize=dearomatize)
        smiles_mol = openbabel.OBMol()
        assert conv.SetInFormat('smi')
        conv.ReadString(smiles_mol, smiles)

        assert conv.SetOutFormat('smi'), 'Unable to set format to SMILES'
        conv.SetOptions('c', conv.OUTOPTIONS)
        smiles = conv.WriteString(smiles_mol).partition('\t')[0]
        assert conv.SetInFormat('smi'), 'Unable to set format to SMILES'
        conv.ReadString(smiles_mol, smiles)

        _, i_n, i_c = self.get_termini(smiles_mol, residue=False)
        if i_n is None or i_c is None:
            smiles_mol_2 = openbabel.OBMol()
            conv.ReadString(smiles_mol_2, smiles)
            _, i_n_2, i_c_2 = self.get_termini(smiles_mol_2)
            if ((i_n_2 is not None) + (i_c_2 is not None)) > ((i_n is not None) + (i_c is not None)):
                smiles = conv.WriteString(smiles_mol_2)

                smiles_mol = openbabel.OBMol()
                assert conv.SetInFormat('smi')
                conv.ReadString(smiles_mol, smiles)
                smiles = conv.WriteString(smiles_mol_2)

                smiles_mol = openbabel.OBMol()
                assert conv.SetInFormat('smi')
                conv.ReadString(smiles_mol, smiles)

                _, i_n, i_c = self.get_termini(smiles_mol, residue=False)

        return (smiles_mol, i_n, i_c)

    def build_from_pdb(self, alphabet, ph=None, major_tautomer=False, dearomatize=False):
        """ Build alphabet from `PDB Chemical Component Dictionary <http://www.wwpdb.org/data/ccd>`_

        Args:
            alphabet (:obj:`Alphabet`): alphabet
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule
        """
        dirname = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'PDB'))
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
        filename = os.path.join(dirname, 'components-pub-xml.tar.gz')
        if not os.path.isfile(filename):
            response = requests.get('http://ligand-expo.rcsb.org/dictionaries/components-pub-xml.tar.gz')
            response.raise_for_status()
            with open(filename, 'wb') as file:
                file.write(response.content)

        smiles_to_monomer = {}
        pdb_ligand_to_monomer = {}
        name_to_monomer = {}
        ambiguous_names = []
        for monomer in alphabet.monomers.values():
            smiles_to_monomer[monomer.export('smiles', options=('c',))] = monomer

            for identifier in monomer.identifiers:
                if identifier.ns == 'pdb.ligand':
                    pdb_ligand_to_monomer[identifier.id] = monomer
                    break

            if monomer.name is not None:
                if monomer.name in name_to_monomer:
                    ambiguous_names.append(monomer.name.lower())
                name_to_monomer[monomer.name.lower()] = monomer
            for synonym in monomer.synonyms:
                if synonym in name_to_monomer:
                    ambiguous_names.append(synonym.lower())
                name_to_monomer[synonym.lower()] = monomer

        for name in ambiguous_names:
            if name in name_to_monomer:
                name_to_monomer.pop(name)

        smiles_conv = openbabel.OBConversion()
        assert smiles_conv.SetInFormat('smi'), 'Unable to set format to SMILES'
        assert smiles_conv.SetOutFormat('smi'), 'Unable to set format to SMILES'
        smiles_conv.SetOptions('c', smiles_conv.OUTOPTIONS)

        inchi_conv = openbabel.OBConversion()
        assert inchi_conv.SetInFormat('smi'), 'Unable to set format to SMILES'
        assert inchi_conv.SetOutFormat('inchi'), 'Unable to set format to InChI'

        ns = '{http://pdbml.pdb.org/schema/pdbx-v40.xsd}'
        with tarfile.open(filename, 'r:gz') as tar_file:
            i_file = 0
            n_files = len(tar_file.getmembers())
            n_monomers = 0
            base_monomers = {}
            same_structures = []
            same_pdb_ids = []
            same_names = []
            potential_incorrect_merges = []
            new_monomers = []
            pdb_ccd_id_to_monomers = {}
            for file_info in tar_file:
                i_file += 1
                if i_file % 1000 == 1:
                    print('Processing file {} of {}'.format(i_file, n_files))

                if os.path.splitext(file_info.name)[-1] != '.xml':
                    continue

                xml_file = tar_file.extractfile(file_info)
                xml_root = ElementTree().parse(xml_file)
                xml_group = xml_root.find(ns + 'chem_compCategory')
                if xml_group is None:
                    continue  # pragma: no cover # element is always present

                # get id
                xml_comp = xml_group.find(ns + 'chem_comp')
                if xml_comp is None:
                    continue  # pragma: no cover # element is always present
                id = xml_comp.get('id')
                identifiers = IdentifierSet([Identifier('pdb-ccd', id)])

                # check that compound has been released, is an amino acid, and is no ambiguous
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

                # get name
                xml_el = xml_comp.find(ns + 'name')
                if xml_el is None:
                    name = None  # pragma: no cover # element is always present
                else:
                    name = xml_el.text.lower()

                # retrieve synonyms
                synonyms = SynonymSet()
                xml_el = xml_comp.find(ns + 'one_letter_code')
                if xml_el is not None:
                    synonyms.add(xml_el.text)

                for xml_group in xml_root.findall(ns + 'pdbx_chem_comp_identifierCategory'):
                    for xml_subgroup in xml_group.findall(ns + 'pdbx_chem_comp_identifier'):
                        if xml_subgroup.get('type') == 'SYSTEMATIC NAME':
                            synonyms.add(xml_subgroup.find(ns + 'identifier').text)

                xml_el = xml_comp.find(ns + 'mon_nstd_parent_comp_id')
                if xml_el is not None:
                    base_monomers[id] = xml_el.text

                # retrieve structure
                smiles = None
                for xml_group in xml_root.findall(ns + 'pdbx_chem_comp_descriptorCategory'):
                    for xml_subgroup in xml_group.findall(ns + 'pdbx_chem_comp_descriptor'):
                        if xml_subgroup.get('type') == 'SMILES_CANONICAL' and \
                                xml_subgroup.get('program') == 'OpenEye OEToolkits':
                            smiles = xml_subgroup.find(ns + 'descriptor').text
                if smiles is None:
                    continue  # pragma: no cover # element is always present

                # discard entries with coordinating metals
                mol = openbabel.OBMol()
                inchi_conv.ReadString(mol, smiles)
                inchi = inchi_conv.WriteString(mol)
                formula = inchi.split('/')[1]
                if '.' in formula:
                    continue

                # correct to residue
                mol = openbabel.OBMol()
                smiles_conv.ReadString(mol, smiles)

                _, i_n, i_c = self.get_termini(mol, residue=False)
                if i_n is None or i_c is None:
                    mol_2 = openbabel.OBMol()
                    smiles_conv.ReadString(mol_2, smiles)
                    _, i_n_2, i_c_2 = self.get_termini(mol_2)
                    if ((i_n_2 is not None) + (i_c_2 is not None)) > ((i_n is not None) + (i_c is not None)):
                        mol = mol_2
                        i_n = i_n_2
                        i_c = i_c_2

                smiles = smiles_conv.WriteString(mol).partition('\t')[0]

                # correct structure for pH
                if ph:
                    try:
                        smiles = get_major_micro_species(smiles, 'smiles', 'smiles', ph,
                                                         major_tautomer=major_tautomer, dearomatize=dearomatize)
                    except jnius.JavaException:
                        continue

                mol = openbabel.OBMol()
                smiles_conv.ReadString(mol, smiles)
                smiles = smiles_conv.WriteString(mol).partition('\t')[0]

                mol = openbabel.OBMol()
                smiles_conv.ReadString(mol, smiles)
                smiles = smiles_conv.WriteString(mol).partition('\t')[0]

                # exclude from alphabet because entry minus O- is equivalent to another entry
                if id in ['ASA', 'GND']:
                    continue

                # merge into alphabet
                monomer = None

                if id == 'ARG':
                    monomer = alphabet.monomers.R
                elif id == 'HIS':
                    monomer = alphabet.monomers.H
                elif id == 'LYS':
                    monomer = alphabet.monomers.K
                elif id not in ['ABA']:
                    monomer = smiles_to_monomer.get(smiles, None)

                if monomer is None:
                    monomer = pdb_ligand_to_monomer.get(id, None)
                    if monomer is None:
                        if id not in ['HSK', 'MTY', 'SUI']:
                            names = list(synonyms)
                            if name is not None:
                                names.append(name)

                            for n in names:
                                monomer = name_to_monomer.get(n.lower(), None)
                                if monomer is not None:
                                    same_names.append((id, monomer.id, n))
                                    break

                        if monomer is None:
                            new_monomers.append(id)
                    else:
                        same_pdb_ids.append((id, smiles, monomer.export('smiles', options=('c',))))
                else:
                    resid_id = ''
                    pdb_ligand_id = ''
                    for identifier in monomer.identifiers:
                        if identifier.ns == 'resid':
                            resid_id = identifier.id
                        if identifier.ns == 'pdb.ligand':
                            pdb_ligand_id = identifier.id
                    same_structures.append((id, resid_id, pdb_ligand_id))

                n_monomers += 1
                if monomer is not None:
                    monomer.synonyms.add(name)
                    monomer.synonyms.update(synonyms)
                    for identifier in monomer.identifiers:
                        if identifier.ns == 'pdb-ccd':
                            potential_incorrect_merges.append((id, identifier.id))
                    monomer.identifiers.update(identifiers)
                else:
                    monomer = Monomer(id=id, name=name, synonyms=synonyms,
                                      identifiers=identifiers, structure=smiles)

                    _, i_n, i_c = self.get_termini(monomer.structure, residue=False)
                    self.set_termini(mol, monomer, i_n, i_c)

                    assert id not in alphabet.monomers
                    alphabet.monomers[id] = monomer

                pdb_ccd_id_to_monomers[id] = monomer

                if n_monomers == self._max_monomers:
                    break

            # set base monomers
            for monomer_id, base_id in base_monomers.items():
                monomer = pdb_ccd_id_to_monomers.get(monomer_id, None)
                base = pdb_ccd_id_to_monomers.get(base_id, None)
                if monomer is not None and base is not None:
                    monomer.base_monomers.add(base)

            # save summary of merging
            filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'PDB', 'merge-report.txt'))
            with open(filename, 'w') as file:
                file.write(('{} monomers were merged into the alphabet because they have the same structures (SMILES):\n'
                            '  PDB-CCD ID\tRESID ID\tPDB Ligand RESID ID\n'
                            '  {}\n\n').format(
                    len(same_structures), '\n  '.join(sorted('\t'.join(n) for n in same_structures))))

                file.write(('{} monomers were merged into the alphabet because they have the same PDB ids (ligand / CCD), '
                            'but different structures (e.g., different stereochemistry or bond order):\n'
                            '  PDB-CCD ID\tPDB-CCD SMILES\tRESID SMILES\n'
                            '  {}\n\n').format(
                    len(same_pdb_ids), '\n  '.join(sorted('\t'.join(n) for n in same_pdb_ids))))

                file.write(('{} monomers were merged into the alphabet because they have the same names:\n'
                            '  PDB-CCD ID\tRESID ID\tName\n'
                            '  {}\n\n').format(
                    len(same_names), '\n  '.join(sorted('\t'.join(n) for n in same_names))))

                file.write(('{} monomers were potentially merged incorrectly:\n'
                            '  PDB-CCD ID-1\tPDB-CCD ID-2\n'
                            '  {}\n\n').format(
                    len(potential_incorrect_merges), '\n  '.join(sorted('\t'.join(n) for n in potential_incorrect_merges))))

                file.write('{} monomers were added to the alphabet:\n  {}\n'.format(
                    len(new_monomers), '\n  '.join(sorted(new_monomers))))

    def build_from_mod(self, alphabet, ph=None, major_tautomer=False, dearomatize=False):
        """ Build alphabet from PSI-MI ontology

        Args:
            alphabet (:obj:`Alphabet`): alphabet
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule
        """
        alphabet.monomers['SGergerOMeCys'] = Monomer(
            id='SGergerOMeCys',
            name='S-geranylgeranyl-L-cysteine methyl ester',
            identifiers=IdentifierSet([Identifier('mod', 'MOD:01119')]),
            base_monomers=[alphabet.monomers.C],
            structure='COC(=O)[C@H](CSC/C=C(/CC/C=C(/CC/C=C(/CCC=C(C)C)\C)\C)\C)[NH3+]',
            l_bond_atoms=[Atom(Monomer, element='N', position=29, charge=-1)],
            l_displaced_atoms=[Atom(Monomer, element='H', position=29),
                                  Atom(Monomer, element='H', position=29, charge=1)])

    def set_termini(self, mol, monomer, i_n, i_c):
        """ Set the C and N terminal bond atoms of a monomer

        Args:
            mol (:obj:`openbabel.OBMol`): molecule
            monomer (:obj:`Monomer`): monomer
            i_n (:obj:`int`): index of N terminus
            i_c (:obj:`int`): index of C terminus
        """
        if i_c:
            c_atom = mol.GetAtom(i_c)
            i_o = None
            for bond in openbabel.OBAtomBondIter(c_atom):
                if bond.GetBondOrder() != 1:
                    continue

                o_atom = bond.GetBeginAtom()
                if o_atom.GetIdx() == c_atom.GetIdx():
                    o_atom = bond.GetEndAtom()

                if o_atom.GetAtomicNum() == 8:
                    other_atoms = sorted([other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(o_atom)])
                    if len(other_atoms) == 1 or other_atoms[-2] == 1:
                        i_o = o_atom.GetIdx()
                        break

            if i_o:
                monomer.r_bond_atoms.append(Atom(Monomer, element='C', position=i_c))
                if o_atom.GetFormalCharge() == 0:
                    monomer.r_displaced_atoms.append(Atom(Monomer, element='O', position=i_o))
                    monomer.r_displaced_atoms.append(Atom(Monomer, element='H', position=i_o))
                else:
                    monomer.r_displaced_atoms.append(Atom(Monomer, element='O', position=i_o, charge=-1))

        if i_n:
            atom = mol.GetAtom(i_n)
            other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(atom)]
            tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(atom)])
            other_atoms += [1] * (3 + atom.GetFormalCharge() - tot_bond_order)
            n_charge = atom.GetFormalCharge()
            n_h = other_atoms.count(1)

            monomer.l_bond_atoms.append(Atom(Monomer, element='N', position=i_n, charge=-n_charge))
            if n_charge >= 0 and n_h >= 1:
                monomer.l_displaced_atoms.append(Atom(Monomer, element='H', position=i_n))
            if n_charge >= 1 and n_h >= 2:
                monomer.l_displaced_atoms.append(Atom(Monomer, element='H', position=i_n, charge=1))

    def get_termini(self, mol, residue=True):
        """ Get indices of atoms of N and C termini

        Args:
            mol (:obj:`openbabel.OBMol`): molecule
            residue (:obj:`bool`, optional): if :obj:`True`, search for a residue (H instead of O- at C terminus)

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
            if self.is_n_terminus(mol, atom):
                atom_ns.append(atom)
            elif self.is_c_terminus(mol, atom, residue=residue):
                atom_cs.append(atom)

        termini = []
        for atom_n in atom_ns:
            for atom_c in atom_cs:
                if self.is_terminus(atom_n, atom_c):
                    termini.append((atom_n, atom_c))

        if not termini:
            if atom_ns and atom_cs:
                termini.append((atom_ns[0], atom_cs[0]))
            elif atom_ns:
                termini.append((atom_ns[0], None))
            elif atom_cs:
                termini.append((None, atom_cs[0]))

        if termini:
            atom_n, atom_c = termini[0]

            if atom_n:
                idx_n = atom_n.GetIdx()
            else:
                idx_n = None

            if atom_c:
                idx_c = atom_c.GetIdx()
                if residue:
                    self.is_c_terminus(mol, atom_c, residue=residue, convert_to_aa=True)
            else:
                idx_c = None
        else:
            idx_n = None
            idx_c = None

        return (mol, idx_n, idx_c)

    def is_n_terminus(self, mol, atom):
        """ Determine if an atom is an N-terminus

        Args:
            mol (:obj:`openbabel.OBMol`): molecule
            atom (:obj:`openbabel.OBAtom`): atom

        Returns:
            :obj:`bool`: :obj:`True` if the atom is an N-terminus
        """
        if atom is None:
            return False

        # check atom is nitrogen
        if atom.GetAtomicNum() != 7:
            return False

        # check atom has charge 0 or +1
        if atom.GetFormalCharge() < 0:
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

        # check atom can form another bond
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(atom)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(atom)])
        other_atoms += [1] * (3 + atom.GetFormalCharge() - tot_bond_order)
        n_charge = atom.GetFormalCharge()
        n_h = other_atoms.count(1)
        if not ((n_charge == 1 and n_h >= 2) or
                (n_charge == 0 and (n_h >= 1 or len(other_atoms) <= 2))):
            return False

        # return True
        return True

    def is_c_terminus(self, mol, atom, residue=True, convert_to_aa=False):
        """ Determine if an atom is an C-terminus

        Args:
            mol (:obj:`openbabel.OBMol`): molecule
            atom (:obj:`openbabel.OBAtom`): atom
            residue (:obj:`bool`, optional): if :obj:`True`, search for a residue (H instead of O- at C terminus)
            convert_to_aa (:obj:`bool`, optional): if :obj:`True`, convert COH to COOH

        Returns:
            :obj:`bool`: :obj:`True` if the atom is an C-terminus
        """
        if atom is None:
            return False

        # check atom is carbon
        if atom.GetAtomicNum() != 6:
            return False

        # check atom has charge 0
        if atom.GetFormalCharge() != 0:
            return False

        # check atom is bonded to 1 oxygen, 1 carbon, and 1 hydrogen atom (residue=True) or 1 oxygen (residue=False) atom
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(atom)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(atom)])
        other_atoms += [1] * (4 + atom.GetFormalCharge() - tot_bond_order)
        other_atoms.sort()
        if set(other_atoms).difference(set([1, 6, 8])):
            return False
        if 6 not in other_atoms or 8 not in other_atoms:
            return False
        if len(other_atoms) != 3:
            return False
        if residue and other_atoms[0] != 1:
            return False
        if not residue and other_atoms[1] != 8:
            return False

        # check bonds to carbon and hydrogen are single bonds and bond to oxygen is a double bond
        o_single_bonds = []
        o_double_bonds = []
        for bond in openbabel.OBAtomBondIter(atom):
            other_atom = bond.GetBeginAtom()
            if other_atom == atom:
                other_atom = bond.GetEndAtom()
            if other_atom.GetAtomicNum() == 8:
                if bond.GetBondOrder() == 1:
                    o_single_bonds.append((bond, other_atom))
                elif bond.GetBondOrder() == 2:
                    o_double_bonds.append((bond, other_atom))
                else:
                    return False
            elif bond.GetBondOrder() != 1:
                return False

        if residue and (len(o_single_bonds) != 0 or len(o_double_bonds) != 1):
            return False
        if not residue and (len(o_single_bonds) != 1 or len(o_double_bonds) != 1):
            return False

        # if not residue, check that single oxygen is bound to H
        if not residue:
            oxygen_atom = o_single_bonds[0][1]
            for bond in openbabel.OBAtomBondIter(oxygen_atom):
                other_atom = bond.GetBeginAtom()
                if other_atom == oxygen_atom:
                    other_atom = bond.GetEndAtom()
                if other_atom.GetIdx() != atom.GetIdx() and other_atom.GetAtomicNum() != 1:
                    return False

        # correct residue to amino acid
        if convert_to_aa:
            h_atom = None
            for other_atom in openbabel.OBAtomAtomIter(atom):
                if other_atom.GetAtomicNum() == 1:
                    h_atom = other_atom
                    break

            if h_atom is not None:
                h_atom.SetAtomicNum(8)
            else:
                o_atom = mol.NewAtom()
                o_atom.SetAtomicNum(8)
                o_atom.SetFormalCharge(0)

                bond = openbabel.OBBond()
                bond.SetBegin(atom)
                bond.SetEnd(o_atom)
                bond.SetBondOrder(1)
                assert mol.AddBond(bond)

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
            identifiers.add(Identifier('resid', id))

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
            backbone=None,
            bond=Bond(
                r_bond_atoms=[Atom(Monomer, element='C', position=None)],
                l_bond_atoms=[Atom(Monomer, element='N', position=None)],
                r_displaced_atoms=[Atom(Monomer, element='O', position=None, charge=-1),
                                       Atom(Monomer, element='H', position=None)],
                l_displaced_atoms=[Atom(Monomer, element='H', position=None),
                                      Atom(Monomer, element='H', position=None, charge=1)]),
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
            backbone=None,
            bond=Bond(
                r_bond_atoms=[Atom(Monomer, element='C', position=None)],
                l_bond_atoms=[Atom(Monomer, element='N', position=None)],
                r_displaced_atoms=[Atom(Monomer, element='O', position=None, charge=-1),
                                       Atom(Monomer, element='H', position=None)],
                l_displaced_atoms=[Atom(Monomer, element='H', position=None),
                                      Atom(Monomer, element='H', position=None, charge=1)]),
            circular=circular)
