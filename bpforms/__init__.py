import pkg_resources

# read version
with open(pkg_resources.resource_filename('bpforms', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()

from .core import (Identifier, IdentifierSet, SynonymSet, Base, BaseSequence,
                   Alphabet, BpForm)
from .dna import dna_alphabet, DnaForm
from .rna import rna_alphabet, RnaForm
from .protein import protein_alphabet, ProteinForm
from .util import get_form
