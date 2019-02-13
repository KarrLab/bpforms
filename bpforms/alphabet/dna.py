""" Alphabet and BpForm to represent modified DNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Alphabet, AlphabetBuilder, Monomer, MonomerSequence, BpForm, Identifier, IdentifierSet, SynonymSet
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from wc_utils.util.chem import EmpiricalFormula
import math
import os
import pkg_resources
import sqlalchemy
import warnings


dna_mod_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'DNAmod.sqlite'))

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

    def load_session(self):
        """ loads an SQLAlchemy session """
        metadata = DeclarativeBase.metadata
        Session = sessionmaker(bind=engine)
        session = Session()
        return session

    def run(self, path=filename):
        """ Build alphabet and, optionally, save to YAML file

        Args:
            path (:obj:`str`, optional): path to save alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        return super(DnaAlphabetBuilder, self).run(path)

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

        invalid_nucleobases = []
        for item in with_inchi.all():
            if item.nameid:
                row = session.query(self.ExpandedAlphabet).filter(self.ExpandedAlphabet.nameid == item.nameid).first()
                if row is None:
                    chars = 'dNMP'
                else:
                    chars = row.Abbreviation
                idx = 0
                tmp = chars
                while chars in alphabet.monomers:
                    idx += 1
                    chars = tmp+'_'+str(idx)

            id = item.chebiname

            name = item.iupacname[1:-1]

            synonyms = SynonymSet()
            # if item.chebiname:
            #     synonyms.add(item.chebiname)

            listOfNames = item.othernames[1:-1].split(', ')
            if listOfNames != ['']:
                for otherName in listOfNames:
                    synonyms.add(otherName[1:-1])

            identifiers = IdentifierSet()
            if item.nameid:
                identifiers.add(Identifier('chebi', item.nameid))

            structure = item.inchi.strip('[]')

            monomer = Monomer(
                id=id,
                name=name,
                synonyms=synonyms,
                identifiers=identifiers,
                structure=structure,
            )

            if not self.is_nucleobase_valid(monomer):
                invalid_nucleobases.append(id)
                continue

            alphabet.monomers[chars] = monomer

            if invalid_nucleobases:
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

    def __init__(self, monomer_seq=None):
        """
        Args:
            monomer_seq (:obj:`MonomerSequence`, optional): monomers of the DNA form
        """
        super(DnaForm, self).__init__(monomer_seq=monomer_seq, alphabet=dna_alphabet,
                                      backbone_formula=EmpiricalFormula('C5H7O6P'), backbone_charge=-2,
                                      bond_formula=EmpiricalFormula('H') * -1, bond_charge=1)


class CanonicalDnaForm(BpForm):
    """ Canonical DNA form """

    def __init__(self, monomer_seq=None):
        """
        Args:
            monomer_seq (:obj:`MonomerSequence`, optional): monomers of the DNA form
        """
        super(CanonicalDnaForm, self).__init__(monomer_seq=monomer_seq, alphabet=canonical_dna_alphabet,
                                               backbone_formula=EmpiricalFormula('C5H7O6P'), backbone_charge=-2,
                                               bond_formula=EmpiricalFormula('H') * -1, bond_charge=1)
