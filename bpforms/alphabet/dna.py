""" Alphabet and BpForm to represent modified DNA

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import Alphabet, AlphabetBuilder, Base, BpForm, Identifier, IdentifierSet, SynonymSet
from wc_utils.util.chem import EmpiricalFormula
import os
import pkg_resources
import warnings
import sqlalchemy
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker


dbFileName = 'DNAmod.sqlite'

filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'dna.yml'))
dna_alphabet = Alphabet().from_yaml(filename)
# :obj:`Alphabet`: Alphabet for DNA nucleotides

canonical_filename = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'dna.canonical.yml'))
canonical_dna_alphabet = Alphabet().from_yaml(canonical_filename)
# :obj:`Alphabet`: Alphabet for canonical DNA nucleotides

engine = sqlalchemy.create_engine('sqlite:///'+dbFileName, echo=True)
declBase = declarative_base(engine)

class DnaAlphabetBuilder(AlphabetBuilder):
    """ Build DNA alphabet from MODOMICS """

    class Names(declBase):
        """"""
        __tablename__ = 'names'
        __table_args__ = {'autoload': True}

    class Expanded_Alphabet(declBase):
        """"""
        __tablename__ = 'expanded_alphabet'
        __table_args__ = {'autoload': True}
        nameid = sqlalchemy.Column(primary_key=True)

    def load_session(self):
        """ loads an SQLAlchemy session """
        metadata = declBase.metadata
        print(32)
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

    def __str__(self):
        return str(self)

    def build(self):
        """ Build alphabet

        Returns:
            :obj:`Alphabet`: alphabet
        """
        # initialize alphabet
        alphabet = Alphabet()

        # create canonical bases
        alphabet.from_yaml(canonical_filename)

        # get individual modifications and create bases
        session = self.load_session()
        with_inchi = session.query(self.Names).filter(self.Names.inchi != '[]').all()
        
        for item in with_inchi:
            if item.nameid:
                row = session.query(self.Expanded_Alphabet).filter(self.Expanded_Alphabet.nameid == item.nameid).first()
                if row is None:
                    chars = 'dNMP'
                else:
                    chars = row.Abbreviation
                print('The chars is {}'.format(chars))
                print(type(chars))
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
                print(type(listOfNames))
                print(102)
                for otherName in listOfNames:
                    print(otherName[1:-1])
                    print(106)
                    synonyms.add(otherName[1:-1])

            identifiers = IdentifierSet()
            if item.nameid:
                identifiers.add(Identifier('ChEBI ID', item.nameid))

            structure = item.inchi.strip('[]')
            print(structure)

            if chars in alphabet.bases:
                warnings.warn('Ignoring canonical base {}'.format(chars), UserWarning)
                continue

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
                                      bond_formula=EmpiricalFormula('H2O') * -1, bond_charge=0) # I need to check with Jonathan if this change is OK!


class CanonicalDnaForm(BpForm):
    """ Canonical DNA form """

    def __init__(self, base_seq=None):
        """
        Args:
            base_seq (:obj:`BaseSequence`, optional): bases of the DNA form
        """
        super(CanonicalDnaForm, self).__init__(base_seq=base_seq, alphabet=canonical_dna_alphabet,
                                               bond_formula=EmpiricalFormula('H2O') * -1, bond_charge=O)




if __name__ == "__main__":
    builder = DnaAlphabetBuilder()
    builder.run()
