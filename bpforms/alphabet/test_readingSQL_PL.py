from bpforms.core import Alphabet, AlphabetBuilder, Base, BpForm, Identifier, IdentifierSet, SynonymSet
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import sqlalchemy

engine = create_engine('sqlite:///DNAmod.sqlite', echo=True)
declBase = declarative_base(engine)
########################################################################


class Names(declBase):
	""""""
	__tablename__ = 'names'
	__table_args__ = {'autoload': True}

print(16)
class Expanded_Alphabet(declBase):
	""""""
	__tablename__ = 'expanded_alphabet'
	__table_args__ = {'autoload': True}
	nameid = sqlalchemy.Column(primary_key=True)

#----------------------------------------------------------------------


def load_session():
	""""""
	metadata = declBase.metadata
	print(28)
	#expanded_alphabet = Table('expanded_alphabet', metadata, 
		#              sqlalchemy.Column("nameid", sqlalchemy.String, primary_key=True),
		#              autoload=True)
	print(32)
	Session = sessionmaker(bind=engine)
	session = Session()
	return session

# get individual modifications and create bases

alphabet = Alphabet()
def gen_alphabet():
	session = load_session()
	with_inchi = session.query(Names).filter(Names.inchi != '[]').all()
	for item in with_inchi:

		if item.nameid:
			id = session.query(Expanded_Alphabet).filter(Expanded_Alphabet.nameid == item.nameid).first()
			id = str(id)
			print(type(id))
			idx = 0
			while id in alphabet.bases:
				idx += 1
				id = id+'_'+str(idx)

		chars = id

		name = item.iupacname

		synonyms = SynonymSet()
		if item.chebiname:
			synonyms.add(item.chebiname)

		if item.othernames:
			listOfNames = item.othernames.split(', ')
			print(type(listOfNames))
			for otherName in listOfNames:
				print(type(otherName))
				synonyms.add(otherName)

		identifiers = IdentifierSet()
		if item.nameid:
			identifiers.add(Identifier('ChEBI ID', item.nameid))

		structure = item.inchi.strip('[]')
		print(structure)

		if chars in alphabet.bases:
			warnings.warn('Ignoring canonical base {}'.format(chars), UserWarning)
			continue

		alphabet.bases[chars] = Base(
			id=chars,
			name=name,
			synonyms=synonyms,
			identifiers=identifiers,
			structure=structure,
			# comments="Modification of {}.".format(mod['originating_base'])
		)

	# return alphabet
	return alphabet


if __name__ == "__main__":
	gen_alphabet()
