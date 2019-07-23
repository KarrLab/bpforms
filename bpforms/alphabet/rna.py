""" Alphabet and BpForm to represent modified RNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import (Alphabet, AlphabetBuilder, Monomer, MonomerSequence, Backbone,
                          Bond, Atom, BpForm, Identifier, IdentifierSet, SynonymSet,
                          BpFormsWarning)
from wc_utils.util.chem import EmpiricalFormula, get_major_micro_species
import bs4
import csv
import io
import openbabel
import os.path
import pkg_resources
import re
import requests
import requests_cache
import warnings

filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'rna.yml'))
rna_alphabet = Alphabet().from_yaml(filename)
# :obj:`Alphabet`: Alphabet for RNA nucleotide monophosphates

canonical_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'rna.canonical.yml'))
canonical_rna_alphabet = Alphabet().from_yaml(canonical_filename)
# :obj:`Alphabet`: Alphabet for canonical RNA nucleotide monophosphates


class RnaAlphabetBuilder(AlphabetBuilder):
    """ Build RNA alphabet from MODOMICS and the RNA Modification Database """

    MODOMICS_INDEX_ENDPOINT = 'http://modomics.genesilico.pl/modifications/'
    MODOMICS_INDEX_ASCII_ENDPOINT = 'http://modomics.genesilico.pl/modifications/?base=all&type=all&display_ascii=Display+as+ASCII'
    MODOMICS_ENTRY_ENDPOINT = 'http://modomics.genesilico.pl/modifications/{}/'
    RNA_MOD_DB_INDEX_ENDPOINT = 'https://mods.rna.albany.edu/mods/modifications/search'
    RNA_MOD_DB_ENTRY_ENDPOINT = 'https://mods.rna.albany.edu/mods/modifications/view/{}'
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
        return super(RnaAlphabetBuilder, self).run(ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize, path=path)

    def build(self, ph=None, major_tautomer=False, dearomatize=False):
        """ Build alphabet

        Args:
            ph (:obj:`float`, optional): pH at which to calculate the major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule

        Returns:
            :obj:`Alphabet`: alphabet
        """
        # initialize alphabet
        alphabet = Alphabet()
        alphabet.id = 'rna'
        alphabet.name = 'RNA nucleotide monophosphates'
        alphabet.description = ('The canonical RNA nucleotide monophosphates, '
                                'plus non-canonical RNA nucleotide monophosphates based on '
                                '<a href="http://modomics.genesilico.pl/modifications">MODOMICS</a> and '
                                'the <a href="https://mods.rna.albany.edu/mods/">RNA Modification Database</a>')

        # create requests session
        cache_name = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'rna'))
        session = requests_cache.core.CachedSession(cache_name, backend='sqlite', expire_after=None)
        session.mount('http://', requests.adapters.HTTPAdapter(max_retries=self.MAX_RETRIES))

        # build from databases
        self.build_modomics(alphabet, session, ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)
        self.build_rna_mod_db(alphabet, session, ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)

        # add missing parent/child information
        if '10G' in alphabet.monomers and 'G' in alphabet.monomers:
            alphabet.monomers['10G'].base_monomers.add(alphabet.monomers['G'])
        if 'manQ' in alphabet.monomers and '10G' in alphabet.monomers:
            alphabet.monomers['manQ'].base_monomers.add(alphabet.monomers['10G'])
        if '102G' in alphabet.monomers and '10G' in alphabet.monomers:
            alphabet.monomers['102G'].base_monomers.add(alphabet.monomers['10G'])
        if '103G' in alphabet.monomers and 'G' in alphabet.monomers:
            alphabet.monomers['103G'].base_monomers.add(alphabet.monomers['G'])
        if '104G' in alphabet.monomers and '10G' in alphabet.monomers:
            alphabet.monomers['104G'].base_monomers.add(alphabet.monomers['10G'])
        if '106G' in alphabet.monomers and '10G' in alphabet.monomers:
            alphabet.monomers['106G'].base_monomers.add(alphabet.monomers['10G'])

        # convert to nucleosides to nucleotide monophosphates
        conv = openbabel.OBConversion()
        assert conv.SetInFormat('smi'), 'Unable to set format to SMILES'
        assert conv.SetOutFormat('smi'), 'Unable to set format to SMILES'
        conv.SetOptions('c', conv.OUTOPTIONS)

        pi_smiles = 'OP([O-])([O-])=O'

        for monomer in alphabet.monomers.values():
            # add PI to molecule
            monomer_structure = openbabel.OBMol(monomer.structure)
            n_monomer_atoms = monomer_structure.NumAtoms()

            pi = openbabel.OBMol()
            conv.ReadString(pi, pi_smiles)

            mol = openbabel.OBMol()
            mol += monomer_structure
            mol += pi

            # get atom references
            monomer_backbone_o = mol.GetAtom(monomer.backbone_bond_atoms[0].position)
            assert monomer_backbone_o.GetAtomicNum() == 8

            monomer_backbone_h = None
            for other_atom in openbabel.OBAtomAtomIter(monomer_backbone_o):
                if other_atom.GetAtomicNum() == 1:
                    monomer_backbone_h = other_atom
                    break

            pi_p = mol.GetAtom(n_monomer_atoms + 2)
            pi_o = mol.GetAtom(n_monomer_atoms + 4)
            pi_oh = mol.GetAtom(n_monomer_atoms + 1)

            # form bond with PI
            bond = openbabel.OBBond()
            bond.SetBegin(monomer_backbone_o)
            bond.SetEnd(pi_p)
            bond.SetBondOrder(1)
            assert mol.AddBond(bond)

            # remove displaced atoms
            mol.DeleteAtom(pi_oh, True)
            if monomer_backbone_h is not None:
                mol.DeleteAtom(monomer_backbone_h, True)

            # canonicalize atom indices
            smiles = conv.WriteString(mol).partition('\t')[0]

            mol2 = openbabel.OBMol()
            conv.ReadString(mol2, smiles)
            smiles = conv.WriteString(mol2).partition('\t')[0]

            mol2 = openbabel.OBMol()
            conv.ReadString(mol2, smiles)
            smiles = conv.WriteString(mol2).partition('\t')[0]

            # get indices of termini
            atom_ls = []
            atom_rs = []
            for i_atom in range(1, mol2.NumAtoms() + 1):
                atom = mol2.GetAtom(i_atom)
                if self.is_l_atom(atom):
                    atom_ls.append(atom)
                elif self.is_r_bond_atom(atom):
                    atom_rs.append(atom)

            i_left_p = None
            i_right_o = None
            i_left_o = None

            for atom_l in atom_ls:
                for atom_r in atom_rs:
                    if self.is_nucleotide_terminus(atom_l, atom_r):
                        i_left_p = atom_l.GetIdx()
                        i_right_o = atom_r.GetIdx()
                        break

            for other_atom in openbabel.OBAtomAtomIter(atom_l):
                if other_atom.GetFormalCharge() == -1:
                    i_left_o = other_atom.GetIdx()

            # update properties of monomer
            monomer.structure = smiles
            monomer.backbone_bond_atoms = []
            monomer.backbone_displaced_atoms = []
            monomer.l_bond_atoms = [Atom(Monomer, element='P', position=i_left_p)]
            monomer.l_displaced_atoms = [Atom(Monomer, element='O', position=i_left_o, charge=-1)]
            monomer.r_bond_atoms = [Atom(Monomer, element='O', position=i_right_o)]
            monomer.r_displaced_atoms = [Atom(Monomer, element='H', position=i_right_o)]

            assert monomer.structure.GetAtom(i_right_o).GetAtomicNum() == 8
            assert monomer.structure.GetAtom(i_left_p).GetAtomicNum() == 15
            assert monomer.structure.GetAtom(i_left_o).GetAtomicNum() == 8

        # return alphabet
        return alphabet

    def build_modomics(self, alphabet, session, ph=None, major_tautomer=False, dearomatize=False):
        """ Build alphabet from MODOMICS

        Args:
            alphabet (:obj:`Alphabet`): alphabet
            session (:obj:`requests_cache.core.CachedSession`): request cache session
            ph (:obj:`float`, optional): pH at which to calculate the major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule
        """
        # get originating monomeric forms
        ascii_response = session.get(self.MODOMICS_INDEX_ASCII_ENDPOINT)
        ascii_response.raise_for_status()
        stream = io.StringIO(ascii_response.text)
        stream.readline()
        reader = csv.reader(stream, delimiter='\t', quoting=csv.QUOTE_NONE)
        base_monomer_short_names = {}
        for monomer in reader:
            for i in range(1, 4):
                short_name = monomer[i]
                if short_name:
                    break
            if short_name in ['A', 'C', 'G', 'U']:
                continue
            base_monomer_short_names[short_name] = monomer[i + 2]

        # get index of nucleosides
        response = session.get(self.MODOMICS_INDEX_ENDPOINT)
        response.raise_for_status()

        # get individual nucleosides and create monomeric forms
        doc = bs4.BeautifulSoup(response.text, 'html.parser')
        table = doc.find('table', {'class': 'datagrid'})
        tbody = table.find('tbody')
        mods = tbody.find_all('tr')
        monomer_short_names = {}
        invalid_nucleosides = []
        for i_mod, mod in enumerate(mods):
            if i_mod >= self._max_monomers:
                break

            cells = mod.find_all('td')
            new_nomenclature = cells[0].text
            name = cells[1].text
            short_name = cells[2].text
            abbrev = cells[3].text

            id = short_name

            chars = new_nomenclature
            if not chars:
                chars = id

            synonyms = SynonymSet()
            if abbrev:
                synonyms.add(abbrev)

            identifiers = IdentifierSet()
            identifiers.add(Identifier('modomics.short_name', short_name))
            if new_nomenclature:
                identifiers.add(Identifier('modomics.new_nomenclature', new_nomenclature))

            structure, more_identifiers = self.get_nucleoside_details_from_modomics(id, session)
            identifiers.update(more_identifiers)

            monomer = Monomer(
                id=chars,
                name=name,
                synonyms=synonyms,
                identifiers=identifiers,
                structure=structure,
            )

            if not monomer.structure:
                continue
            if '*' in monomer.export('smiles', options=('c',)):
                continue
            if ' (base)' in new_nomenclature:
                continue
            if ' (cap)' in name or ' cap' in name:
                continue
            if '-CoA)' in name:
                continue
            if new_nomenclature and new_nomenclature[-1] == 'N':
                continue

            conv = openbabel.OBConversion()
            assert conv.SetOutFormat('smi'), 'Unable to set format to SMILES'
            assert conv.SetInFormat('smi'), 'Unable to set format to SMILES'
            conv.SetOptions('c', conv.OUTOPTIONS)

            smiles = conv.WriteString(monomer.structure).partition('\t')[0]
            if ph is not None:
                smiles = get_major_micro_species(smiles, 'smiles', 'smiles', ph, major_tautomer=major_tautomer, dearomatize=dearomatize)
            smiles_mol = openbabel.OBMol()
            conv.ReadString(smiles_mol, smiles)

            smiles = conv.WriteString(smiles_mol).partition('\t')[0]
            smiles_mol2 = openbabel.OBMol()
            conv.ReadString(smiles_mol2, smiles)

            monomer.structure = smiles_mol2

            if not self.is_valid_nucleoside(monomer):
                invalid_nucleosides.append(chars)
                continue

            alphabet.monomers[chars] = monomer
            monomer_short_names[short_name] = monomer

        if invalid_nucleosides:
            warnings.warn('The following compounds were ignored because they do not appear to be nucleosides:\n- {}'.format(
                '\n- '.join(invalid_nucleosides)), BpFormsWarning)

        for short_name, monomer in monomer_short_names.items():
            base_monomer = monomer_short_names.get(base_monomer_short_names.get(short_name, None), None)
            if base_monomer:
                monomer.base_monomers.add(base_monomer)

        # return alphabet
        return alphabet

    def get_nucleoside_details_from_modomics(self, id, session):
        """ Get the structure of a nucleoside in the MODOMICS database

        Args:
            id (:obj:`str`): id of nucleoside in MODOMICS database

        Returns:
            :obj:`openbabel.OBMol`: structure
            :obj:`IdentifierSet`: identifiers
        """
        response = session.get(self.MODOMICS_ENTRY_ENDPOINT.format(id))
        response.raise_for_status()

        doc = bs4.BeautifulSoup(response.text, 'html.parser')

        table = doc.find(id='modification_details')
        tbody = table.find('tbody')
        rows = tbody.find_all('tr')
        mol = None
        identifiers = IdentifierSet()
        for row in rows:
            cells = row.find_all('td')
            if cells[0].text.startswith('SMILES'):
                smiles = cells[1].text
                if not smiles:
                    continue

                mol = openbabel.OBMol()
                conv = openbabel.OBConversion()
                assert conv.SetInFormat('smi')
                if not conv.ReadString(mol, smiles):
                    mol = None
                    continue

            elif cells[0].text.startswith('PubChem'):
                link = cells[1].find('a')
                identifiers.add(Identifier('pubchem.compound', link.text))

        return mol, identifiers

    def build_rna_mod_db(self, alphabet, session, ph=None, major_tautomer=False, dearomatize=False):
        """ Build alphabet from the RNA Modification Database

        Args:
            alphabet (:obj:`Alphabet`): alphabet
            session (:obj:`requests_cache.core.CachedSession`): request cache session
            ph (:obj:`float`, optional): pH at which to calculate the major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule
        """
        ########################################
        # get information from the RNA Modification Database

        # parse index
        response = session.get(self.RNA_MOD_DB_INDEX_ENDPOINT)
        response.raise_for_status()
        doc = bs4.BeautifulSoup(response.text, 'html.parser')
        table = doc.find('table', {'class': 'searchresult'})
        entry_ids = []
        for row in table.find_all('tr')[1:]:
            entry_ids.append(row.find('td').find('a').get('href').partition('/mods/modifications/view/')[2])

        # parse entries
        monomers = []
        for entry_id in entry_ids:
            response = session.get(self.RNA_MOD_DB_ENTRY_ENDPOINT.format(entry_id))
            response.raise_for_status()

            id = None
            name = None
            identifiers = IdentifierSet([Identifier('rnamods', str(entry_id))])
            comments = []

            doc = bs4.BeautifulSoup(response.text, 'html.parser')
            div = doc.find('div', {'class': 'modView'})

            dl = div.find('dl')
            key = None
            for child in dl.children:
                if isinstance(child, bs4.element.Tag):
                    if child.name == 'dt':
                        key = child.text[0:-1]
                    elif child.name == 'dd':
                        if key == 'Symbol':
                            id = child.text
                        elif key == 'Common name':
                            name = child.text
                        elif key == 'CA registry numbers' and child.text.startswith('ribonucleoside '):
                            cas_id = child.text[len('ribonucleoside '):].strip()
                            if cas_id:
                                identifiers.add(Identifier('cas', cas_id))

            # comment
            for p in div.find_all('p'):
                label = p.find('span', {'class': 'refLabel'}).text
                if label == 'Comment:':
                    comments.append('<p>{}</p>'.format(re.sub(' \[\d+\]', '', p.text.partition('Comment:')[2].strip())))

            # phylogenetic prevalance
            table = div.find('table', {'class': 'phylotable'})
            if table:
                rows = table.find_all('tr')

                domains = []
                for cell in rows[1].find_all('th'):
                    domains.append(cell.text.strip())

                domain_dist = []
                for row in rows[2:]:
                    rna_type_obj = row.find('th')
                    if rna_type_obj:
                        rna_type = rna_type_obj.text.strip()
                        assert rna_type.endswith('RNA')
                        for domain, cell in zip(domains, row.find_all('td')):
                            if cell.text:
                                domain_dist.append(domain + ' ' + rna_type)

                comments.append('<p>Phylogenetic distribution<ul><li>{}</li></ul></p>'.format(
                    '</li><li>'.join(domain_dist)))

            monomers.append(Monomer(id=id, name=name, identifiers=identifiers,
                                    comments=' '.join(comments)))

        ########################################
        # corrections to the RNA Modification Database: incorrect or missing CAS ids
        for monomer in monomers:
            # incorrect ids
            for identifier in monomer.identifiers:
                if identifier.ns == 'cas':
                    if identifier.id == '577773-09-02':
                        identifier.id = '577773-09-2'

                    elif identifier.id == '-56973-12-7':
                        identifier.id = '56973-12-7'

                    elif monomer.id == 'C+':
                        identifier.id = '1221169-70-5'

            # missing ids
            if monomer.id == 'cnm5U':
                monomer.identifiers.add(Identifier('cas', '58479-73-5'))

            elif monomer.id == 'gcmnm5s2U':
                monomer.identifiers.add(Identifier('cas', '1401688-11-6'))

            elif monomer.id == 'gmnm5s2U':
                monomer.identifiers.add(Identifier('cas', '1401688-10-5'))

        ########################################
        # get structures for entries from SciFinder

        # save list of CAS ids to search for in SciFinder
        rnamods_dirname = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'rnamods'))
        if not os.path.isdir(rnamods_dirname):
            os.makedirs(rnamods_dirname)  # pragma: no cover # already created to store .mol files

        with open(os.path.join(rnamods_dirname, 'index.tsv'), 'w') as file:
            csv_writer = csv.writer(file, dialect='excel-tab', quoting=csv.QUOTE_MINIMAL)
            csv_writer.writerow(['RNA Modification Database Id', 'RNA Modification Database Symbol', 'Name', 'CAS Id', 'Comments'])
            for monomer in monomers:
                rnamods_id = ''
                cas_id = ''
                for identifier in monomer.identifiers:
                    if identifier.ns == 'rnamods':
                        rnamods_id = identifier.id
                    elif identifier.ns == 'cas':
                        cas_id = identifier.id
                csv_writer.writerow([rnamods_id, monomer.id, monomer.name, cas_id, monomer.comments])

        # Get structures from SciFinder and save to file (done manually at https://scifinder.cas.org)
        # 1. Use `Substance Identifier` to find compounds by their CAS id
        # 2. Select `Explore by Structure` > `Chemical Structure` in the context menu for the compound
        # 3. Open the `Non-Java` structure viewer
        # 4. Click the save icon
        # 5. Save the structure in .mol format with the name equal to the CAS id

        # Read structures from file
        mol_conv = openbabel.OBConversion()
        assert mol_conv.SetInFormat('mol'), 'Unable to set format to MOL'

        smiles_conv = openbabel.OBConversion()
        assert smiles_conv.SetInFormat('smi'), 'Unable to set format to SMILES'
        assert smiles_conv.SetOutFormat('smi'), 'Unable to set format to SMILES'
        smiles_conv.SetOptions('c', smiles_conv.OUTOPTIONS)

        for monomer in monomers:
            for identifier in monomer.identifiers:
                if identifier.ns == 'cas':
                    cas_id = identifier.id
                    break

            # read .mol file
            mol_filename = os.path.join(rnamods_dirname, cas_id + '.mol')
            assert os.path.isfile(mol_filename)
            mol = openbabel.OBMol()
            mol_conv.ReadFile(mol, mol_filename)
            smiles = smiles_conv.WriteString(mol)

            # get major microspecies
            if ph:
                smiles = get_major_micro_species(smiles, 'smiles', 'smiles', ph, major_tautomer=major_tautomer, dearomatize=dearomatize)

            # sanitize SMILES
            mol = openbabel.OBMol()
            smiles_conv.ReadString(mol, smiles)
            smiles = smiles_conv.WriteString(mol)

            mol = openbabel.OBMol()
            smiles_conv.ReadString(mol, smiles)
            smiles = smiles_conv.WriteString(mol)

            # set structure
            monomer.structure = smiles

        ########################################
        # Determine valid monomers and find 3' and 5' sites
        for monomer in list(monomers):
            assert self.is_valid_nucleoside(monomer), "Monomer {} does not appear to be a valid nucleoside".format(monomer.id)

        ########################################
        # merge monomers with alphabet
        additional_monomers = []
        for monomer in monomers:
            inchi = monomer.export('inchi').partition('/h')[0]
            same_monomer = None
            for test_monomer in alphabet.monomers.values():
                if monomer.id == 'manQ':
                    break

                modomics_short_name = None
                for identifier in test_monomer.identifiers:
                    if identifier.ns == 'modomics.short_name':
                        modomics_short_name = identifier.id
                        break
                if monomer.id == modomics_short_name:
                    same_monomer = test_monomer
                    break

                if test_monomer.export('inchi').partition('/h')[0] == inchi:
                    same_monomer = test_monomer
                    break

            if same_monomer:
                if monomer.id not in [same_monomer.id, same_monomer.name]:
                    same_monomer.synonyms.add(monomer.id)

                if monomer.name not in [same_monomer.id, same_monomer.name]:
                    same_monomer.synonyms.add(monomer.name)

                same_monomer.identifiers.update(monomer.identifiers)

                assert not same_monomer.comments, "Comments must be merged, which isn't implemented"
                same_monomer.comments = monomer.comments
            else:
                code = monomer.id
                assert code not in alphabet.monomers, "Code already used. Another code must be chosen."
                alphabet.monomers[code] = monomer
                additional_monomers.append(code)

            if additional_monomers:
                print('{} monomeric forms were added to the alphabet:\n  {}'.format(
                    len(additional_monomers), '\n  '.join(additional_monomers)))

    def is_valid_nucleoside(self, monomer):
        """ Determine if nucleoside should be included in alphabet

        Args:
            monomer (:obj:`Monomer`): monomeric form

        Returns:
            :obj:`bool`: :obj:`True` if the monomeric form is a valid nucleoside
        """
        formula = monomer.get_formula()
        if formula.C < 9 or formula.O < 4 or formula.N < 2:
            return False

        atom_bs = []
        atom_rs = []
        for i_atom in range(1, monomer.structure.NumAtoms() + 1):
            atom = monomer.structure.GetAtom(i_atom)
            if self.is_backbone_atom(atom):
                atom_bs.append(atom)
            elif self.is_r_bond_atom(atom):
                atom_rs.append(atom)

        termini = []
        for atom_b in atom_bs:
            for atom_r in atom_rs:
                if self.is_terminus(atom_b, atom_r):
                    termini.append((atom_b, atom_r))

        if termini:
            atom_b_idx = termini[0][0].GetIdx()
            atom_r_idx = termini[0][1].GetIdx()
        else:  # pragma no cover: case not used by MODOMICS or the RNA Modification Database
            if atom_bs and not atom_rs:
                atom_b_idx = atom_bs[0].GetIdx()
                atom_r_idx = None
            elif atom_rs and not atom_bs:
                atom_b_idx = None
                atom_r_idx = atom_rs[0].GetIdx()
            else:
                atom_b_idx = None
                atom_r_idx = None

        if atom_b_idx:
            monomer.backbone_bond_atoms = [Atom(Monomer, 'O', atom_b_idx)]
            monomer.backbone_displaced_atoms = [Atom(Monomer, 'H', atom_b_idx)]
        if atom_r_idx:
            monomer.r_bond_atoms = [Atom(Monomer, 'O', atom_r_idx)]
            monomer.r_displaced_atoms = [Atom(Monomer, 'H', atom_r_idx)]

        return True

    def is_terminus(self, b_atom, r_atom):
        """ Determine if a pair of atoms is a valid pair of bonding sites

        Args:
            b_atom (:obj:`openbabel.OBAtom`): potential backbone atom
            r_atom (:obj:`openbabel.OBAtom`): potential right bond atom

        Returns:
            :obj:`bool`: :obj:`True`, if the atoms are a valid pair of bonding sites
        """
        r_atom_2 = self.is_backbone_atom(b_atom)
        if not r_atom_2:
            return False  # pragma no cover: case not used by MODOMICS or the RNA Modification Database

        if not self.is_r_bond_atom(r_atom):
            return False  # pragma no cover: case not used by MODOMICS or the RNA Modification Database

        return r_atom_2.GetIdx() == r_atom.GetIdx()

    def is_nucleotide_terminus(self, l_atom, r_atom):
        """ Determine if a pair of atoms is a valid pair of bonding sites

        Args:
            l_atom (:obj:`openbabel.OBAtom`): potential left atom
            r_atom (:obj:`openbabel.OBAtom`): potential right bond atom

        Returns:
            :obj:`bool`: :obj:`True`, if the atoms are a valid pair of bonding sites
        """
        r_atom_2 = self.is_l_atom(l_atom)
        if not r_atom_2:
            return False  # pragma no cover: case not used by MODOMICS or the RNA Modification Database

        if not self.is_r_bond_atom(r_atom):
            return False  # pragma no cover: case not used by MODOMICS or the RNA Modification Database

        return r_atom_2.GetIdx() == r_atom.GetIdx()

    def is_backbone_atom(self, b_atom):
        """ Determine if an atom is a valid backbone bonding site

        Args:
            b_atom (:obj:`openbabel.OBAtom`): potential backbone atom

        Returns:
            :obj:`bool`: :obj:`True`, if the atom is a valid backbone bonding site
        """
        if b_atom.GetAtomicNum() != 8:
            return False
        if b_atom.GetFormalCharge() != 0:
            return False

        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(b_atom)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(b_atom)])
        other_atoms += [1] * (2 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        if other_atoms != [1, 6]:
            return False

        # get first C
        for other_atom in openbabel.OBAtomAtomIter(b_atom):
            if other_atom.GetAtomicNum() == 6:
                c_1 = other_atom
        if c_1.GetFormalCharge() != 0:
            return False  # pragma no cover: case not used by MODOMICS or the RNA Modification Database
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(c_1)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(c_1)])
        other_atoms += [1] * (4 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        if other_atoms != [1, 1, 6, 8]:
            return False

        # get second C
        for other_atom in openbabel.OBAtomAtomIter(c_1):
            if other_atom.GetAtomicNum() == 6:
                c_2 = other_atom
        if c_2.GetFormalCharge() != 0:
            return False  # pragma no cover: case not used by MODOMICS or the RNA Modification Database
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(c_2)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(c_2)])
        other_atoms += [1] * (4 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        if other_atoms != [1, 6, 6, 8]:
            return False

        # get third C
        for other_atom in openbabel.OBAtomAtomIter(c_2):
            if other_atom.GetAtomicNum() == 6 and other_atom.GetIdx() != c_1.GetIdx():
                c_3 = other_atom
        if c_3.GetFormalCharge() != 0:
            return False  # pragma no cover: case not used by MODOMICS or the RNA Modification Database
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(c_3)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(c_3)])
        other_atoms += [1] * (4 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        if other_atoms != [1, 6, 6, 8]:
            return False

        # get second O
        for other_atom in openbabel.OBAtomAtomIter(c_3):
            if other_atom.GetAtomicNum() == 8:
                r_atom_2 = other_atom
        if r_atom_2.GetFormalCharge() != 0:
            return False  # pragma no cover: case not used by MODOMICS or the RNA Modification Database

        return r_atom_2

    def is_l_atom(self, l_atom):
        """ Determine if an atom is a valid left bonding site

        Args:
            l_atom (:obj:`openbabel.OBAtom`): potential left atom

        Returns:
            :obj:`bool`: :obj:`True`, if the atom is a valid left bonding site
        """
        if l_atom.GetAtomicNum() != 15:
            return False
        if l_atom.GetFormalCharge() != 0:
            return False

        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(l_atom)]
        if other_atoms != [8, 8, 8, 8]:
            return False
        bonds = sorted([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(l_atom)])
        if bonds != [1, 1, 1, 2]:
            return False
        c_1 = None
        o_minus = 0
        for bond in openbabel.OBAtomBondIter(l_atom):
            other_atom = bond.GetBeginAtom()
            if other_atom == l_atom:
                other_atom = bond.GetEndAtom()

            if other_atom.GetFormalCharge() == -1 and len(list(openbabel.OBAtomAtomIter(other_atom))) == 1:
                o_minus += 1
                continue

            for other_other_atom in openbabel.OBAtomAtomIter(other_atom):
                if other_other_atom.GetAtomicNum() == 6:
                    c_1 = other_other_atom
                    break
        if o_minus != 2:
            return False
        if c_1 is None:
            return False

        # get first C
        if c_1.GetFormalCharge() != 0:
            return False  # pragma no cover: case not used by MODOMICS or the RNA Modification Database
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(c_1)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(c_1)])
        other_atoms += [1] * (4 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        if other_atoms != [1, 1, 6, 8]:
            return False

        # get second C
        for other_atom in openbabel.OBAtomAtomIter(c_1):
            if other_atom.GetAtomicNum() == 6:
                c_2 = other_atom
        if c_2.GetFormalCharge() != 0:
            return False  # pragma no cover: case not used by MODOMICS or the RNA Modification Database
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(c_2)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(c_2)])
        other_atoms += [1] * (4 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        if other_atoms != [1, 6, 6, 8]:
            return False

        # get third C
        for other_atom in openbabel.OBAtomAtomIter(c_2):
            if other_atom.GetAtomicNum() == 6 and other_atom.GetIdx() != c_1.GetIdx():
                c_3 = other_atom
        if c_3.GetFormalCharge() != 0:
            return False  # pragma no cover: case not used by MODOMICS or the RNA Modification Database
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(c_3)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(c_3)])
        other_atoms += [1] * (4 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        if other_atoms != [1, 6, 6, 8]:
            return False

        # get second O
        for other_atom in openbabel.OBAtomAtomIter(c_3):
            if other_atom.GetAtomicNum() == 8:
                r_atom_2 = other_atom
        if r_atom_2.GetFormalCharge() != 0:
            return False  # pragma no cover: case not used by MODOMICS or the RNA Modification Database

        return r_atom_2

    def is_r_bond_atom(self, r_atom):
        """ Determine if an atom is a valid right bond bonding site

        Args:
            b_atom (:obj:`openbabel.OBAtom`): potential right bond atom

        Returns:
            :obj:`bool`: :obj:`True`, if the atom is a valid right bond bonding site
        """
        other_atoms = [other_atom.GetAtomicNum() for other_atom in openbabel.OBAtomAtomIter(r_atom)]
        tot_bond_order = sum([bond.GetBondOrder() for bond in openbabel.OBAtomBondIter(r_atom)])
        other_atoms += [1] * (2 - tot_bond_order)
        other_atoms = sorted(other_atoms)
        if other_atoms != [1, 6]:
            return False

        return True


class RnaForm(BpForm):
    """ RNA form """

    DEFAULT_FASTA_CODE = 'N'

    def __init__(self, seq=None, circular=False):
        """
        Args:
            seq (:obj:`MonomerSequence`, optional): sequence of monomeric forms of the DNA
            circular (:obj:`bool`, optional): if :obj:`True`, indicates that the biopolymer is circular
        """
        super(RnaForm, self).__init__(
            seq=seq, alphabet=rna_alphabet,
            backbone=None,
            bond=Bond(
                r_bond_atoms=[Atom(Monomer, element='O', position=None)],
                l_bond_atoms=[Atom(Monomer, element='P', position=None)],
                r_displaced_atoms=[Atom(Monomer, element='H', position=None)],
                l_displaced_atoms=[Atom(Monomer, element='O', position=None, charge=-1)]),
            circular=circular)


class CanonicalRnaForm(BpForm):
    """ Canonical RNA form """

    DEFAULT_FASTA_CODE = 'N'

    def __init__(self, seq=None, circular=False):
        """
        Args:
            seq (:obj:`MonomerSequence`, optional): sequence of monomeric forms of the DNA
            circular (:obj:`bool`, optional): if :obj:`True`, indicates that the biopolymer is circular
        """
        super(CanonicalRnaForm, self).__init__(
            seq=seq, alphabet=canonical_rna_alphabet,
            backbone=None,
            bond=Bond(
                r_bond_atoms=[Atom(Monomer, element='O', position=None)],
                l_bond_atoms=[Atom(Monomer, element='P', position=None)],
                r_displaced_atoms=[Atom(Monomer, element='H', position=None)],
                l_displaced_atoms=[Atom(Monomer, element='O', position=None, charge=-1)]),
            circular=circular)
