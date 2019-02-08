from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
 
engine = create_engine('sqlite:///DNAmod.sqlite', echo=True)
declBase = declarative_base(engine)
########################################################################
class Names(declBase):
    """"""
    __tablename__ = 'names'
    __table_args__ = {'autoload':True}

class Expanded_alphabet(declBase):
    """"""
    __tablename__ = 'expanded_alphabet'
    __table_args__ = {'autoload':True}
 
#----------------------------------------------------------------------
def loadSession():
    """"""
    metadata = declBase.metadata
    Session = sessionmaker(bind=engine)
    session = Session()
    return session
 
if __name__ == "__main__":
    session = loadSession()
    with_inchi = session.query(Names).filter(Names.inchi != '[]').all()
    print(res[1].othernames)


# get individual modifications and create bases
for item in with_inchi:

    if item.nameid:
        id = session.query(Expanded_alphabet).filter(Expanded_alphabet.nameid == item.nameid).one()
        idx = 0
        while id in alphabet.bases:
            idx+=1
            id = id+'_'+str(idx)

    chars = id

    name = item.iupacname

    synonyms = SynonymSet()
    if item.chebiname:
        synonyms.add(item.chebiname)
    if item.othernames:
        for name in item.othernemes:
            synonyms.add(item.othernames[name])

    identifiers = IdentifierSet()
    if item.nameid:
        identifiers.add(Identifier('ChEBI ID', item.namid))
    
    structure = item.inchi

    if chars in alphabet.bases:
        warnings.warn('Ignoring canonical base {}'.format(chars), UserWarning)
        continue

    alphabet.bases[chars] = Base(
        id=chars,
        name=mod['name'],
        synonyms=synonyms,
        identifiers=identifiers,
        structure=structure,
        # comments="Modification of {}.".format(mod['originating_base'])
    )

# return alphabet
return alphabet