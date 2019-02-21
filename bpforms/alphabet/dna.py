""" Alphabet and BpForm to represent modified DNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import (Alphabet, AlphabetBuilder, Monomer, MonomerSequence, BpForm,
                          Backbone, Bond, Atom, Identifier, IdentifierSet, SynonymSet)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from wc_utils.util.chem import EmpiricalFormula
import math
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
# :obj:`Alphabet`: Alphabet for DNA nucleobases

canonical_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'dna.canonical.yml'))
canonical_dna_alphabet = Alphabet().from_yaml(canonical_filename)
# :obj:`Alphabet`: Alphabet for canonical DNA nucleobases

engine = sqlalchemy.create_engine('sqlite:///' + dna_mod_filename)
DeclarativeBase = declarative_base(engine)


class DnaAlphabetBuilder(AlphabetBuilder):
    """ Build DNA alphabet from MODOMICS """

    INVALID_NAMES = (
        'adenosine', 'cytidine', 'guanosine', 'uridine', 'side',
        'ribose', 'ribulose', 'phosphate',
        'antelmycin', 'cytosylglucuronic acid', '1-(beta-D-xylopyranosyl)cytosine',
        'telbivudine', 'stavudine', 'tenofovir',
        '9-{2,5-anhydro-4-[(phosphonooxy)methyl]-alpha-L-lyxofuranosyl}-9H-purin-6-amine',
        'oxetanocin',
        'adefovir', '5-amino-6-(5-phospho-beta-D-ribosylamino)uracil',
        "2'-deoxyinosine",
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

    def run(self, ph=None, major_tautomer=False, path=filename):
        """ Build alphabet and, optionally, save to YAML file

        Args:
            ph (:obj:`float`, optional): pH at which calculate major protonation state of each monomer
            major_tautomer (:obj:`bool`, optional): if :obj:`True`, calculate the major tautomer
            path (:obj:`str`, optional): path to save alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        return super(DnaAlphabetBuilder, self).run(ph=ph, major_tautomer=major_tautomer, path=path)

    def build(self):
        """ Build alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        # initialize alphabet
        alphabet = Alphabet()

        # create canonical monomers
        alphabet.from_yaml(canonical_filename)
        alphabet.id = 'dna'
        alphabet.name = 'DNAmod DNA nucleobases'
        alphabet.description = ('The four canonical DNA nucleobases, plus the non-canonical DNA nucleobases in '
                                '<a href="https://dnamod.hoffmanlab.org/modifications">DNAmod</a>')

        # get individual nucleobases and create monomers
        session = self.load_session()
        with_inchi = session.query(self.Names).filter(self.Names.inchi != '[]')
        if not math.isinf(self._max_monomers):
            with_inchi = with_inchi.limit(self._max_monomers)

        monomer_nameids = {}
        all_base_monomers_codes = {}
        all_base_monomers_ids = {}
        invalid_nucleobases = []
        for item in with_inchi.all():
            row = session.query(self.ExpandedAlphabet).filter(self.ExpandedAlphabet.nameid == item.nameid).first()
            if row is None:
                chars = 'base'
            elif row.Symbol:
                chars = row.Symbol
            else:
                chars = row.Abbreviation
            idx = 0
            tmp = chars
            while chars in alphabet.monomers:
                idx += 1
                chars = tmp + '-' + str(idx)

            id = item.chebiname

            name = item.iupacname[1:-1]

            synonyms = SynonymSet()
            other_names = item.othernames[1:-1].split(', ')
            if other_names != ['']:
                for other_name in other_names:
                    synonyms.add(other_name[1:-1])

            identifiers = IdentifierSet()
            identifiers.add(Identifier('chebi', item.nameid))

            structure = item.inchi.strip('[]')

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

            monomer = Monomer(
                id=id,
                name=name,
                synonyms=synonyms,
                identifiers=identifiers,
                structure=structure,
                comments=comments,
            )

            if not self.is_nucleobase_valid(monomer):
                invalid_nucleobases.append(id)
                continue

            alphabet.monomers[chars] = monomer

            monomer_nameids[item.nameid] = monomer
            all_base_monomers_codes[monomer] = base_monomer_codes
            all_base_monomers_ids[monomer] = base_monomer_ids

        for monomer, base_monomer_codes in all_base_monomers_codes.items():
            for base_monomer_code in base_monomer_codes:
                monomer.base_monomers.add(alphabet.monomers[base_monomer_code])
        for monomer, base_monomer_ids in all_base_monomers_ids.items():
            for base_monomer_id in base_monomer_ids:
                base_monomer = monomer_nameids.get(base_monomer_id, None)
                if base_monomer:
                    monomer.base_monomers.add(base_monomer)

        warnings.warn('The following compounds were ignored because they do not appear to be nucleobases:\n- {}'.format(
            '\n- '.join(invalid_nucleobases)), UserWarning)

        return alphabet

    def is_nucleobase_valid(self, monomer):
        """ Determine if monomer should be included in alphabet

        Args:
            monomer (:obj:`Monomer`): monomer

        Returns:
            :obj:`bool`: :obj:`True` if monomer should be included in alphabet
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

    def __init__(self, monomer_seq=None):
        """
        Args:
            monomer_seq (:obj:`MonomerSequence`, optional): monomers of the DNA form
        """
        super(DnaForm, self).__init__(
            monomer_seq=monomer_seq, alphabet=dna_alphabet,
            backbone=Backbone(
                structure='InChI=1/C5H11O6P/c6-4-1-2-10-5(4)3-11-12(7,8)9/h4-6H,1-3H2,(H2,7,8,9)/p-2',
                monomer_bond_atoms=[Atom(element='N', position=None)],
                backbone_bond_atoms=[Atom(element='C', position=2)],
                monomer_displaced_atoms=[Atom(element='H', position=None)],
                backbone_displaced_atoms=[Atom(element='H', position=2)]),
            bond=Bond(
                left_participant=Monomer,
                right_participant=Monomer,
                left_bond_atoms=[Atom(element='O', position=6)],
                right_bond_atoms=[Atom(element='P', position=12)],
                left_displaced_atoms=[Atom(element='H', position=6)],
                right_displaced_atoms=[Atom(element='O', position=8, charge=-1)]))


class CanonicalDnaForm(BpForm):
    """ Canonical DNA form """

    DEFAULT_FASTA_CODE = 'N'

    def __init__(self, monomer_seq=None):
        """
        Args:
            monomer_seq (:obj:`MonomerSequence`, optional): monomers of the DNA form
        """
        super(CanonicalDnaForm, self).__init__(
            monomer_seq=monomer_seq, alphabet=canonical_dna_alphabet,
            backbone=Backbone(
                structure='InChI=1/C5H11O6P/c6-4-1-2-10-5(4)3-11-12(7,8)9/h4-6H,1-3H2,(H2,7,8,9)/p-2',
                monomer_bond_atoms=[Atom(element='N', position=None)],
                backbone_bond_atoms=[Atom(element='C', position=2)],
                monomer_displaced_atoms=[Atom(element='H', position=None)],
                backbone_displaced_atoms=[Atom(element='H', position=2)]),
            bond=Bond(
                left_participant=Monomer,
                right_participant=Monomer,
                left_bond_atoms=[Atom(element='O', position=6)],
                right_bond_atoms=[Atom(element='P', position=12)],
                left_displaced_atoms=[Atom(element='H', position=6)],
                right_displaced_atoms=[Atom(element='O', position=8, charge=-1)]))
