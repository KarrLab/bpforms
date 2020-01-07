""" Alphabet and BpForm to represent modified DNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import (Alphabet, AlphabetBuilder, Monomer, MonomerSequence, BpForm,
                          Bond, Atom, Identifier, IdentifierSet, SynonymSet,
                          BpFormsWarning)
from bpforms.alphabet.core import (download_pdb_ccd, parse_pdb_ccd, get_pdb_ccd_open_babel_mol,
                                   get_can_smiles)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from wc_utils.util.chem import EmpiricalFormula
from wc_utils.util.chem.marvin import get_major_micro_species
import csv
import math
import openbabel
import os
import pkg_resources
import requests
import sqlalchemy
import warnings


dna_mod_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'DNAmod.sqlite'))


def get_dnamod(filename):
    if not os.path.isfile(filename):
        response = requests.get('https://dnamod.hoffmanlab.org/DNAmod.sqlite')
        with open(filename, 'wb') as file:
            file.write(response.content)


get_dnamod(dna_mod_filename)

filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'dna.yml'))
dna_alphabet = Alphabet().from_yaml(filename)
# :obj:`Alphabet`: Alphabet for DNA nucleotide monophosphates

canonical_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'dna.canonical.yml'))
canonical_dna_alphabet = Alphabet().from_yaml(canonical_filename)
# :obj:`Alphabet`: Alphabet for canonical DNA nucleotide monophosphates

engine = sqlalchemy.create_engine('sqlite:///' + dna_mod_filename)
DeclarativeBase = declarative_base(engine)


class DnaAlphabetBuilder(AlphabetBuilder):
    """ Build DNA alphabet from MODOMICS """

    INVALID_IDS = (
        # nucleosides
        "2'-deoxyinosine",
        "2'-deoxyuridine",
        "2'-deoxy-5-(4,5-dihydroxypentyl)uridine",
        "5-(2-aminoethoxy)methyl-2'-deoxyuridine",
        "5,6-dihydroxy-2'-deoxyuridine",
        "5-(2-aminoethyl)-2'-deoxyuridine",
        "7-amido-7-deazaguanosine",
        "(beta-D-glucopyranosyloxymethyl)deoxyuridine 5'-monophosphate",
        "N(2),N(2)-dimethylguanosine",
    )

    class Names(DeclarativeBase):
        """"""
        __tablename__ = 'names'
        __table_args__ = {'autoload': True}

    class ExpandedAlphabet(DeclarativeBase):
        """"""
        __tablename__ = 'expanded_alphabet'
        __table_args__ = {'autoload': True}
        nameid = sqlalchemy.Column(primary_key=True)

    class ModBase(DeclarativeBase):
        """"""
        __tablename__ = 'modbase'
        __table_args__ = {'autoload': True, 'extend_existing': True}
        nameid = sqlalchemy.Column(primary_key=True)

    class ModBaseParents(DeclarativeBase):
        """"""
        __tablename__ = 'modbase_parents'
        __table_args__ = {'autoload': True, 'extend_existing': True}
        nameid = sqlalchemy.Column(primary_key=True)

    class CovMod(DeclarativeBase):
        """"""
        __tablename__ = 'covmod'
        __table_args__ = {'autoload': True, 'extend_existing': True}
        cmodid = sqlalchemy.Column(primary_key=True)

    def load_session(self):
        """ loads an SQLAlchemy session """
        metadata = DeclarativeBase.metadata
        Session = sessionmaker(bind=engine)
        session = Session()
        return session

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
        return super(DnaAlphabetBuilder, self).run(ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize, path=path)

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

        # create canonical monomeric forms
        alphabet.from_yaml(canonical_filename)
        alphabet.id = 'dna'
        alphabet.name = 'DNA nucleotide monophosphates'
        alphabet.description = ('The canonical DNA nucleotide monophosphates, '
                                'plus the non-canonical DNA nucleotide monophosphates based on '
                                '<a href="https://dnamod.hoffmanlab.org">DNAmod</a> and '
                                '<a href="http://repairtoire.genesilico.pl/damage/">REPAIRtoire</a>')

        # build from sources
        self.build_dnamod(alphabet, ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)
        self.build_repairtoire(alphabet, ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)
        self.build_pdb_ccd(alphabet, ph=ph, major_tautomer=major_tautomer, dearomatize=dearomatize)

        # corrections
        if 'dI' in alphabet.monomers and 'DI' in alphabet.monomers:
            dI = alphabet.monomers.dI
            alphabet.monomers.dI = alphabet.monomers.pop('DI')
            alphabet.monomers.dI.synonyms.add("2'-deoxyinosine, 9-[(2R,4S,5R)-4-hydroxy-5-"
                                              "(hydroxymethyl)tetrahydrofuran-2-yl]-9H-purin-6-ol")
            alphabet.monomers.dI.synonyms.add('deoxyinosine')
            alphabet.monomers.dI.synonyms.add('9-(2-deoxy-beta-d-erythro-pentofuranosyl)-9h-purin-6-ol')
            alphabet.monomers.dI.identifiers.add(Identifier('dnamod', "2'-deoxyinosine"))
            alphabet.monomers.dI.identifiers.add(Identifier('chebi', 'CHEBI:28997'))
            alphabet.monomers.dI.comments = ("A purine 2'-deoxyribonucleoside that is inosine in which the hydroxy "
                                             "group at position 2' is replaced by a hydrogen.")
            for monomer in alphabet.monomers.values():
                if dI in monomer.base_monomers:
                    monomer.base_monomers.remove(dI)
                    monomer.base_monomers.add(alphabet.monomers.dI)

        if 'dU' in alphabet.monomers:
            alphabet.monomers.dU.synonyms.add("2'-deoxyuridine")
            alphabet.monomers.dU.synonyms.add('1-(2-deoxy-beta-D-erythro-pentofuranosyl)uracil')
            alphabet.monomers.dU.synonyms.add("2'-deoxyuridine-5'-monophosphate")
            alphabet.monomers.dU.synonyms.add("dU")
            alphabet.monomers.dU.synonyms.add("deoxyuridine")
            alphabet.monomers.dU.synonyms.add("2-deoxyuridine")
            alphabet.monomers.dU.identifiers.add(Identifier('dnamod', "2'-deoxyuridine"))
            alphabet.monomers.dU.identifiers.add(Identifier('chebi', 'CHEBI:16450'))
            alphabet.monomers.dU.identifiers.add(Identifier('repairtoire', 'dU'))
            alphabet.monomers.dU.comments = (
                "A pyrimidine 2'-deoxyribonucleoside having uracil as the nucleobase. Cysteine "
                "spontaneously loses an amine group, which is replaced by a keto group at the "
                "C4-atom. This reaction occurs spontaneously and all the time, and is the reason "
                "why T is used in DNA. This way, dU and T can be distinguished from each other.")

        # return alphabet
        return alphabet

    def build_dnamod(self, alphabet, ph=None, major_tautomer=False, dearomatize=False):
        """ Build monomeric forms from DNAmod

        Args:
            alphabet (:obj:`Alphabet`): alphabet
            ph (:obj:`float`, optional): pH at which to calculate the major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule

        Returns:
            :obj:`Alphabet`: alphabet
        """
        filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'dnamod.csv'))
        monomers = {}
        with open(filename, 'r') as file:
            reader = csv.DictReader(file, dialect='excel')
            for row in reader:
                monomers[row['Id']] = {
                    'nucleobase': row['Nucleobase'],
                    'nucleotide': row['Nucleotide'],
                    'l_bond_atom': int(float(row['Left bond atom (P)'])),
                    'l_displaced_atom': int(float(row['Left displaced atom (O-)'])),
                    'r_bond_atom': int(float(row['Right bond atom (O)'])),
                }

        # get individual nucleobases and create monomeric forms
        session = self.load_session()
        with_smiles = session.query(self.Names).filter(self.Names.smiles != '[]')
        if not math.isinf(self._max_monomers):
            with_smiles = with_smiles.limit(self._max_monomers)

        monomer_nameids = {}
        all_base_monomers_codes = {}
        all_base_monomers_ids = {}
        invalid_nucleobases = []
        for item in with_smiles.all():
            row = session.query(self.ExpandedAlphabet).filter(self.ExpandedAlphabet.nameid == item.nameid).first()
            if row is None:
                continue
            if row.Symbol:
                chars = row.Symbol
            else:
                chars = row.Abbreviation

            id = item.chebiname

            name = item.iupacname[1:-1]

            synonyms = SynonymSet()
            other_names = item.othernames[1:-1].split(', ')
            if other_names != ['']:
                for other_name in other_names:
                    synonyms.add(other_name[1:-1])

            identifiers = IdentifierSet()
            identifiers.add(Identifier('chebi', item.nameid))
            identifiers.add(Identifier('dnamod', id))

            smiles = item.smiles.strip('[]')
            inchi = item.inchi.strip('[]')

            cmodid = None
            verified_status = False
            base_monomer_codes = set()
            for base_monomer_code in session.query(self.ModBase).filter(self.ModBase.nameid == item.nameid).all():
                cmodid = base_monomer_code.cmodid
                verified_status = base_monomer_code.verifiedstatus
                if self.ModBase.baseid != 'O':
                    base_monomer_codes.add(base_monomer_code.baseid)

            if not verified_status:
                continue

            base_monomer_ids = set()
            for base_monomer_id in session.query(self.ModBaseParents).filter(self.ModBaseParents.nameid == item.nameid).all():
                base_monomer_ids.add(base_monomer_id.parentid)

            comments = session.query(self.CovMod).filter(self.CovMod.cmodid == cmodid).first().definition

            mol = openbabel.OBMol()
            conv = openbabel.OBConversion()
            assert conv.SetInFormat('smi')
            if not conv.ReadString(mol, smiles):
                assert conv.SetInFormat('inchi')
                assert conv.ReadString(mol, inchi)
            assert conv.SetOutFormat('smiles')
            conv.SetOptions('c', conv.OUTOPTIONS)
            smiles = conv.WriteString(mol, True)

            assert conv.SetInFormat('smiles')

            if ph:
                smiles = get_major_micro_species(smiles, 'smiles', 'smiles', ph=ph,
                                                 major_tautomer=major_tautomer, dearomatize=dearomatize)

            mol = openbabel.OBMol()
            conv.ReadString(mol, smiles)
            smiles = conv.WriteString(mol, True)

            mol = openbabel.OBMol()
            conv.ReadString(mol, smiles)
            smiles = conv.WriteString(mol, True)

            monomer = alphabet.monomers.get(chars, Monomer())
            monomer.id = id
            monomer.name = name
            monomer.synonyms = synonyms
            monomer.identifiers = identifiers
            monomer.structure = smiles
            monomer.comments = comments

            if not self.is_nucleobase_valid(monomer):
                invalid_nucleobases.append(id)
                continue

            if monomer.id in monomers:
                if smiles != monomers[monomer.id]['nucleobase']:
                    raise Exception('Structure and atom indices may need to be updated for {}\n  {}\n  {}'.format(
                        monomer.id, smiles, monomers[monomer.id]['nucleobase']))

                smiles = monomers[monomer.id]['nucleotide']
                if ph:
                    smiles = get_major_micro_species(smiles, 'smiles', 'smiles', ph=ph,
                                                     major_tautomer=major_tautomer, dearomatize=dearomatize)
                mol = openbabel.OBMol()
                conv.ReadString(mol, smiles)
                smiles = conv.WriteString(mol, True)
                mol = openbabel.OBMol()
                conv.ReadString(mol, smiles)
                smiles = conv.WriteString(mol, True)

                if smiles != monomers[monomer.id]['nucleotide']:
                    raise Exception('Structure and atom indices may need to be updated for {}\n  {}\n  {}'.format(
                        monomer.id, smiles, monomers[monomer.id]['nucleotide']))

                monomer.structure = smiles
                monomer.l_bond_atoms = [Atom(Monomer, element='P', position=monomers[monomer.id]['l_bond_atom'])]
                monomer.l_displaced_atoms = [Atom(Monomer, element='O', position=monomers[monomer.id]
                                                  ['l_displaced_atom'], charge=-1)]
                monomer.r_bond_atoms = [Atom(Monomer, element='O', position=monomers[monomer.id]['r_bond_atom'])]
                monomer.r_displaced_atoms = [Atom(Monomer, element='H', position=monomers[monomer.id]['r_bond_atom'])]
            else:
                raise Exception('Structure and atom indices need to be cataloged for {}'.format(monomer.id))

            alphabet.monomers[chars] = monomer

            monomer_nameids[item.nameid] = monomer
            all_base_monomers_codes[monomer] = base_monomer_codes
            all_base_monomers_ids[monomer] = base_monomer_ids

        for monomer, base_monomer_codes in all_base_monomers_codes.items():
            for base_monomer_code in base_monomer_codes:
                base_monomer = alphabet.monomers.get(base_monomer_code, None)
                if base_monomer:
                    monomer.base_monomers.add(base_monomer)
        for monomer, base_monomer_ids in all_base_monomers_ids.items():
            for base_monomer_id in base_monomer_ids:
                base_monomer = monomer_nameids.get(base_monomer_id, None)
                if base_monomer:
                    monomer.base_monomers.add(base_monomer)

        if invalid_nucleobases:
            warnings.warn('The following compounds were ignored because they do not appear to be nucleobases:\n- {}'.format(
                '\n- '.join(invalid_nucleobases)), BpFormsWarning)

    def build_repairtoire(self, alphabet, ph=None, major_tautomer=False, dearomatize=False):
        """ Build monomeric forms from DNAmod

        Args:
            alphabet (:obj:`Alphabet`): alphabet
            ph (:obj:`float`, optional): pH at which to calculate the major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule

        Returns:
            :obj:`Alphabet`: alphabet
        """
        monomer_smiles = {}
        for monomer in alphabet.monomers.values():
            monomer_smiles[monomer.export('smiles', options=('c',))] = monomer

        filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'repairtoire.csv'))
        conv = openbabel.OBConversion()
        assert conv.SetInFormat('smiles')
        assert conv.SetOutFormat('smiles')
        conv.SetOptions('c', conv.OUTOPTIONS)
        merged_monomers = []
        new_monomers = []
        with open(filename, 'r') as file:
            reader = csv.DictReader(file, dialect='excel')
            for row in reader:
                id = row['Id']
                name = row['Name']
                smiles = row['Nucleotide monophosphate (2-) (cleaned)']
                i_l_bond_atom = int(float(row['Left bond atom (P)']))
                i_l_displaced_atom = int(float(row['Left displaced atom (O-)']))
                i_r_bond_atom = int(float(row['Right bond atom (O)']))
                comments = row['Comments']

                smiles2 = smiles
                if ph:
                    smiles2 = get_major_micro_species(smiles2, 'smiles', 'smiles', ph=ph,
                                                      major_tautomer=major_tautomer, dearomatize=dearomatize)
                mol = openbabel.OBMol()
                conv.ReadString(mol, smiles2)
                smiles2 = conv.WriteString(mol, True)
                mol = openbabel.OBMol()
                conv.ReadString(mol, smiles2)
                smiles2 = conv.WriteString(mol, True)
                if smiles2 != smiles:
                    raise Exception('Structure and atom indices may need to be updated for {}'.format(id))

                monomer = monomer_smiles.get(smiles, None)
                if monomer is None:
                    monomer = alphabet.monomers.get(id, None)

                if monomer is not None:
                    monomer.synonyms.add(id)
                    monomer.synonyms.add(name)
                    monomer.identifiers.add(Identifier('repairtoire', id))
                    if comments:
                        if monomer.comments:
                            monomer.comments += ' ' + comments
                        else:
                            monomer.comments = comments

                    merged_monomers.append((id, monomer.id))

                else:
                    monomer = Monomer(
                        id=id, name=name,
                        identifiers=[Identifier('repairtoire', id)],
                        structure=smiles,
                        l_bond_atoms=[Atom(Monomer, element='P', position=i_l_bond_atom)],
                        l_displaced_atoms=[Atom(Monomer, element='O', position=i_l_displaced_atom, charge=-1)],
                        r_bond_atoms=[Atom(Monomer, element='O', position=i_r_bond_atom)],
                        r_displaced_atoms=[Atom(Monomer, element='H', position=i_r_bond_atom)],
                        comments=comments or None)
                    alphabet.monomers[id] = monomer

                    new_monomers.append(id)

        print('Merged monomers\n  {}'.format('\n  '.join('\t'.join(ids) for ids in merged_monomers)))
        print('New monomers\n  {}'.format('\n  '.join(new_monomers)))

    def is_nucleobase_valid(self, monomer):
        """ Determine if monomeric form should be included in alphabet

        Args:
            monomer (:obj:`Monomer`): monomeric form

        Returns:
            :obj:`bool`: :obj:`True` if monomeric form should be included in alphabet
        """
        for invalid_id in self.INVALID_IDS:
            if invalid_id in monomer.id:
                return False

        return True

    def build_pdb_ccd(self, alphabet, ph=None, major_tautomer=False, dearomatize=False):
        """ Build monomeric forms from PDB CCD

        Args:
            alphabet (:obj:`Alphabet`): alphabet
            ph (:obj:`float`, optional): pH at which to calculate the major protonation state of each monomeric form
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            dearomatize (:obj:`bool`, optional): if :obj:`True`, dearomatize molecule

        Returns:
            :obj:`Alphabet`: alphabet
        """
        smiles_to_monomer = {get_can_smiles(monomer.structure): monomer for monomer in alphabet.monomers.values()}
        monomer_codes = {code: monomer for code, monomer in alphabet.monomers.items()}
        monomer_lc_codes = {code.lower(): code for code in alphabet.monomers.keys()}
        monomer_to_codes = {monomer: code for code, monomer in alphabet.monomers.items()}

        base_monomers = {}
        pdb_id_to_monomer = {}

        same_structures = []
        same_ids = []
        same_ids_case_insensitive = []

        replaced_monomers = []

        filename = download_pdb_ccd()
        valid_types = ('DNA linking',
                       'DNA OH 5 prime terminus',
                       'DNA OH 3 prime terminus')
        for pdb_monomer, base_monomer, smiles, pdb_structure, atoms in \
                parse_pdb_ccd(filename, valid_types, self._max_monomers):
            structure = get_pdb_ccd_open_babel_mol(pdb_structure)
            if structure is None:
                continue

            # set structure, atoms
            mol = openbabel.OBMol()
            mol += structure

            if "P" in atoms:
                atom = mol.GetAtom(atoms["P"]['position'])
                assert atom.GetAtomicNum() == 15
                atom.SetIsotope(1)

            ops = []
            for i_o in range(1, 4):
                id_o = 'OP' + str(i_o)
                id_h = 'H' + id_o
                if id_o in atoms and id_h in atoms:
                    atom_o = mol.GetAtom(atoms[id_o]['position'])
                    assert atom_o.GetAtomicNum() == 8
                    atom_o.SetIsotope(len(ops))

                    atom_h = mol.GetAtom(atoms[id_h]['position'])
                    assert atom_h.GetAtomicNum() == 1

                    ops.append((atom_o, atom_h))
                    if len(ops) == 2:
                        break

            if "O3'" in atoms:
                atom = mol.GetAtom(atoms["O3'"]['position'])
                assert atom.GetAtomicNum() == 8
                atom.SetIsotope(2)

            for op, hop in ops:
                assert mol.DeleteAtom(hop, True)
                op.SetFormalCharge(-1)

            conv = openbabel.OBConversion()
            assert conv.SetInFormat('smi')
            assert conv.SetOutFormat('smi')
            conv.SetOptions('c', conv.OUTOPTIONS)
            isotope_smiles = conv.WriteString(mol, True)

            mol = openbabel.OBMol()
            conv = openbabel.OBConversion()
            assert conv.SetInFormat('smi')
            assert conv.SetOutFormat('smi')
            conv.SetOptions('c', conv.OUTOPTIONS)
            conv.ReadString(mol, isotope_smiles)
            isotope_smiles = conv.WriteString(mol, True)

            mol = openbabel.OBMol()
            conv = openbabel.OBConversion()
            assert conv.SetInFormat('smi')
            assert conv.SetOutFormat('smi')
            conv.SetOptions('c', conv.OUTOPTIONS)
            conv.ReadString(mol, isotope_smiles)
            i_left_p = None
            i_right_o = None
            i_left_o = None
            for atom in openbabel.OBMolAtomIter(mol):
                if atom.GetAtomicNum() == 15 and atom.GetIsotope() == 1:
                    i_left_p = atom.GetIdx()
                    atom.SetIsotope(0)
                elif atom.GetAtomicNum() == 8 and atom.GetIsotope() == 1:
                    i_left_o = atom.GetIdx()
                    atom.SetIsotope(0)
                elif atom.GetAtomicNum() == 8 and atom.GetIsotope() == 2:
                    i_right_o = atom.GetIdx()
                    atom.SetIsotope(0)
            can_smiles = conv.WriteString(mol, True)

            atom_map = {}
            for atom in openbabel.OBMolAtomIter(mol):
                if atom.GetAtomicNum() > 1:
                    atom_map[atom.GetIdx()] = len(atom_map) + 1

            mol = openbabel.OBMol()
            conv = openbabel.OBConversion()
            assert conv.SetInFormat('smi')
            assert conv.SetOutFormat('smi')
            conv.SetOptions('c', conv.OUTOPTIONS)
            conv.ReadString(mol, can_smiles)
            can_smiles = conv.WriteString(mol, True)

            mol = openbabel.OBMol()
            conv = openbabel.OBConversion()
            assert conv.SetInFormat('smi')
            assert conv.SetOutFormat('smi')
            conv.SetOptions('c', conv.OUTOPTIONS)
            conv.ReadString(mol, can_smiles)

            atom_map_2 = {}
            for atom in openbabel.OBMolAtomIter(mol):
                if atom.GetAtomicNum() > 1:
                    atom_map_2[len(atom_map_2) + 1] = atom.GetIdx()

            if i_left_p is not None:
                i_left_p = atom_map_2[atom_map[i_left_p]]
                assert mol.GetAtom(i_left_p).GetAtomicNum() == 15
            if i_left_o is not None:
                i_left_o = atom_map_2[atom_map[i_left_o]]
                assert mol.GetAtom(i_left_o).GetAtomicNum() == 8
            if i_right_o is not None:
                i_right_o = atom_map_2[atom_map[i_right_o]]
                assert mol.GetAtom(i_right_o).GetAtomicNum() == 8

            pdb_monomer.structure = mol
            if i_left_p is not None and i_left_o is not None and \
                    ('HOP3' in atoms or mol.GetAtom(i_left_o).GetFormalCharge() == -1):
                pdb_monomer.l_bond_atoms = [Atom(Monomer, element='P', position=i_left_p)]
                pdb_monomer.l_displaced_atoms = [
                    Atom(Monomer, element='O', position=i_left_o, charge=-1)]
            if i_right_o is not None and "HO3'" in atoms:
                pdb_monomer.r_bond_atoms = [Atom(Monomer, element='O', position=i_right_o)]
                pdb_monomer.r_displaced_atoms = [Atom(Monomer, element='H', position=i_right_o)]

            # get canonical SMILES for entry
            can_smiles = smiles

            if ph:
                can_smiles = get_major_micro_species(
                    can_smiles, 'smiles', 'smiles', ph=ph,
                    major_tautomer=major_tautomer,
                    dearomatize=dearomatize)

            mol = openbabel.OBMol()
            conv = openbabel.OBConversion()
            assert conv.SetInFormat('smi')
            assert conv.SetOutFormat('smi')
            conv.SetOptions('c', conv.OUTOPTIONS)
            conv.ReadString(mol, can_smiles)
            can_smiles = conv.WriteString(mol, True)

            mol = openbabel.OBMol()
            conv = openbabel.OBConversion()
            assert conv.SetInFormat('smi')
            assert conv.SetOutFormat('smi')
            conv.SetOptions('c', conv.OUTOPTIONS)
            conv.ReadString(mol, can_smiles)
            can_smiles = conv.WriteString(mol, True)

            mol = openbabel.OBMol()
            conv = openbabel.OBConversion()
            assert conv.SetInFormat('smi')
            assert conv.SetOutFormat('smi')
            conv.SetOptions('c', conv.OUTOPTIONS)
            conv.ReadString(mol, can_smiles)
            can_smiles = get_can_smiles(mol)

            # determine if the entry is already represented in the alphabet
            merge_with_new = False

            monomer = None

            if pdb_monomer.id not in ['A3A', 'DNR']:
                monomer = smiles_to_monomer.get(can_smiles, None)
                if monomer is not None:
                    for identifier in monomer.identifiers:
                        if identifier.ns == 'pdb-ccd':
                            monomer = None
                            break
                if monomer is not None and (\
                    (monomer.l_bond_atoms and not pdb_monomer.l_bond_atoms) or \
                    (monomer.r_bond_atoms and not pdb_monomer.r_bond_atoms)):
                    monomer = None

            if monomer is not None:
                merge_with_new = True
                same_structures.append((pdb_monomer.id, monomer_to_codes[monomer]))
            if pdb_monomer.id in ['DA', 'DC', 'DG', 'DT']:
                merge_with_new = False
                monomer = alphabet.monomers.get(pdb_monomer.id[1], None)
            if monomer is None:
                monomer = monomer_codes.get(pdb_monomer.id, None)
                if monomer is not None:
                    same_ids.append(pdb_monomer.id)
            if monomer is None:
                monomer_diff_case = monomer_lc_codes.get(pdb_monomer.id.lower(), None)
                if monomer_diff_case is not None:
                    same_ids_case_insensitive.append((pdb_monomer.id, monomer_diff_case))

            # add/merge the entry with the alphabet
            if monomer is None:
                # add the entry to the alphabet
                alphabet.monomers[pdb_monomer.id] = pdb_monomer
                if base_monomer is not None:
                    base_monomers[pdb_monomer] = base_monomer
                pdb_id_to_monomer[pdb_monomer.id] = pdb_monomer

            elif merge_with_new:
                # merge an existing monomer with the PDB entry
                alphabet.monomers[monomer_to_codes[monomer]] = pdb_monomer

                pdb_monomer.synonyms.update(monomer.synonyms)
                pdb_monomer.synonyms.add(pdb_monomer.id)
                pdb_monomer.synonyms.add(pdb_monomer.name)
                pdb_monomer.id = monomer.id or pdb_monomer.id
                pdb_monomer.name = monomer.name or pdb_monomer.name
                pdb_monomer.synonyms.discard(pdb_monomer.id)
                pdb_monomer.synonyms.discard(pdb_monomer.name)

                pdb_monomer.identifiers.update(monomer.identifiers)
                pdb_monomer.base_monomers.update(monomer.base_monomers)
                pdb_monomer.comments = monomer.comments

                if base_monomer is not None:
                    base_monomers[pdb_monomer] = base_monomer
                pdb_id_to_monomer[pdb_monomer.id] = pdb_monomer

                replaced_monomers.append((monomer, pdb_monomer))

            else:
                # merge the entry with an existing monomer
                if pdb_monomer.id not in [monomer.id, monomer.name]:
                    monomer.synonyms.add(pdb_monomer.id)
                if pdb_monomer.name not in [monomer.id, monomer.name]:
                    monomer.synonyms.add(pdb_monomer.name)
                monomer.synonyms.update(pdb_monomer.synonyms)
                monomer.identifiers.update(pdb_monomer.identifiers)

                if base_monomer is not None:
                    base_monomers[monomer] = base_monomer
                pdb_id_to_monomer[pdb_monomer.id] = monomer

        # set base monomers
        for monomer, base_monomer in base_monomers.items():
            if base_monomer in pdb_id_to_monomer:
                monomer.base_monomers.add(pdb_id_to_monomer[base_monomer])

        for old_monomer, new_monomer in replaced_monomers:
            for o_monomer in alphabet.monomers.values():
                if old_monomer in o_monomer.base_monomers:
                    o_monomer.base_monomers.remove(old_monomer)
                    o_monomer.base_monomers.add(new_monomer)

        # print summary
        print('{} entries with the similar structures were joined:\n  {}\t{}\n  {}'.format(
            len(same_structures), 'PDB CCD id', 'Alphabet code',
            '\n  '.join(sorted('\t'.join(ids) for ids in same_structures))))
        print('{} entries with the same ids were joined:\n  {}'.format(
            len(same_ids), '\n  '.join(sorted(same_ids))))
        print('{} entries with similar ids potentially should be joined:\n  {}\t{}\n  {}'.format(
            len(same_ids_case_insensitive), 'PDB CCD id', 'Alphabet code',
            '\n  '.join(sorted('\t'.join(ids) for ids in same_ids_case_insensitive))))


class DnaForm(BpForm):
    """ DNA form """

    DEFAULT_FASTA_CODE = 'N'

    def __init__(self, seq=None, circular=False):
        """
        Args:
            seq (:obj:`MonomerSequence`, optional): monomeric forms of the DNA form
            circular (:obj:`bool`, optional): if :obj:`True`, indicates that the biopolymer is circular
        """
        super(DnaForm, self).__init__(
            seq=seq, alphabet=dna_alphabet,
            backbone=None,
            bond=Bond(
                r_bond_atoms=[Atom(Monomer, element='O', position=None)],
                l_bond_atoms=[Atom(Monomer, element='P', position=None)],
                r_displaced_atoms=[Atom(Monomer, element='H', position=None)],
                l_displaced_atoms=[Atom(Monomer, element='O', position=None, charge=-1)]),
            circular=circular)


class CanonicalDnaForm(BpForm):
    """ Canonical DNA form """

    DEFAULT_FASTA_CODE = 'N'

    def __init__(self, seq=None, circular=False):
        """
        Args:
            seq (:obj:`MonomerSequence`, optional): monomeric forms of the DNA form
            circular (:obj:`bool`, optional): if :obj:`True`, indicates that the biopolymer is circular
        """
        super(CanonicalDnaForm, self).__init__(
            seq=seq, alphabet=canonical_dna_alphabet,
            backbone=None,
            bond=Bond(
                r_bond_atoms=[Atom(Monomer, element='O', position=None)],
                l_bond_atoms=[Atom(Monomer, element='P', position=None)],
                r_displaced_atoms=[Atom(Monomer, element='H', position=None)],
                l_displaced_atoms=[Atom(Monomer, element='O', position=None, charge=-1)]),
            circular=circular)
