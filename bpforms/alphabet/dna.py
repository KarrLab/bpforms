""" Alphabet and BpForm to represent modified DNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Alphabet, AlphabetBuilder, Base, BpForm, Identifier, IdentifierSet, SynonymSet
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
# :obj:`Alphabet`: Alphabet for DNA nucleotides

canonical_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'dna.canonical.yml'))
canonical_dna_alphabet = Alphabet().from_yaml(canonical_filename)
# :obj:`Alphabet`: Alphabet for canonical DNA nucleotides

engine = sqlalchemy.create_engine('sqlite:///' + dna_mod_filename, echo=True)
DeclarativeBase = declarative_base(engine)


class DnaAlphabetBuilder(AlphabetBuilder):
    """ Build DNA alphabet from MODOMICS """

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

        # create canonical bases
        alphabet.from_yaml(canonical_filename)
        alphabet.id = 'dna'
        alphabet.name = 'DNA'
        alphabet.description = ('The four canonical bases, plus the modified bases in '
                                '<a href="https://dnamod.hoffmanlab.org/modifications">DNAmod</a>')

        # get individual modifications and create bases
        session = self.load_session()
        with_inchi = session.query(self.Names).filter(self.Names.inchi != '[]')
        if not math.isinf(self._max_bases):
            with_inchi = with_inchi.limit(self._max_bases)

        for item in with_inchi.all():
            if item.nameid:
                row = session.query(self.ExpandedAlphabet).filter(self.ExpandedAlphabet.nameid == item.nameid).first()
                if row is None:
                    chars = 'dNMP'
                else:
                    chars = row.Abbreviation
                idx = 0
                tmp = chars
                while chars in alphabet.bases:
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
                identifiers.add(Identifier('ChEBI ID', item.nameid))

            structure = item.inchi.strip('[]')

            alphabet.bases[chars] = Base(
                id=id,
                name=name,
                synonyms=synonyms,
                identifiers=identifiers,
                structure=structure,
                # comments="Modification of {}.".format(mod['originating_base'])
            )

        return alphabet


class DnaForm(BpForm):
    """ DNA form """

    def __init__(self, base_seq=None):
        """
        Args:
            base_seq (:obj:`BaseSequence`, optional): bases of the DNA form
        """
        super(DnaForm, self).__init__(base_seq=base_seq, alphabet=dna_alphabet,
                                      bond_formula=EmpiricalFormula('H') * -1, bond_charge=1)


class CanonicalDnaForm(BpForm):
    """ Canonical DNA form """

    def __init__(self, base_seq=None):
        """
        Args:
            base_seq (:obj:`BaseSequence`, optional): bases of the DNA form
        """
        super(CanonicalDnaForm, self).__init__(base_seq=base_seq, alphabet=canonical_dna_alphabet,
                                               bond_formula=EmpiricalFormula('H') * -1, bond_charge=1)
