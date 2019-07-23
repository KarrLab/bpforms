""" Alphabet and BpForm to represent modified DNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import (Alphabet, AlphabetBuilder, Monomer, MonomerSequence, BpForm,
                          Bond, Atom, Identifier, IdentifierSet, SynonymSet,
                          BpFormsWarning)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from wc_utils.util.chem import EmpiricalFormula, get_major_micro_species
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

    INVALID_NAMES = ()

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
            smiles = conv.WriteString(mol).partition('\t')[0]

            assert conv.SetInFormat('smiles')

            if ph:
                smiles = get_major_micro_species(smiles, 'smiles', 'smiles', ph=ph,
                                                 major_tautomer=major_tautomer, dearomatize=dearomatize)

            mol = openbabel.OBMol()
            conv.ReadString(mol, smiles)
            smiles = conv.WriteString(mol).partition('\t')[0]

            mol = openbabel.OBMol()
            conv.ReadString(mol, smiles)
            smiles = conv.WriteString(mol).partition('\t')[0]

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
                smiles = conv.WriteString(mol).partition('\t')[0]
                mol = openbabel.OBMol()
                conv.ReadString(mol, smiles)
                smiles = conv.WriteString(mol).partition('\t')[0]

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
                smiles2 = conv.WriteString(mol).partition('\t')[0]
                mol = openbabel.OBMol()
                conv.ReadString(mol, smiles2)
                smiles2 = conv.WriteString(mol).partition('\t')[0]
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
        for invalid_name in self.INVALID_NAMES:
            if invalid_name in monomer.name:
                return False
            for synonym in monomer.synonyms:
                if invalid_name in synonym:
                    return False

        return True


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
