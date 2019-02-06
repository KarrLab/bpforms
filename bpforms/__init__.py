import pkg_resources

# read version
with open(pkg_resources.resource_filename('bpforms', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()

from .core import (Identifier, IdentifierSet, SynonymSet, Base, BaseSequence,
                   Alphabet, BpForm)
from . import alphabet
from .alphabet.dna import dna_alphabet, DnaForm
from .alphabet.rna import rna_alphabet, RnaForm
from .alphabet.protein import protein_alphabet, ProteinForm
from . import util
from .util import get_form
