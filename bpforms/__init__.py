import pkg_resources

# read version
with open(pkg_resources.resource_filename('bpforms', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()

from .core import (Identifier, IdentifierSet, SynonymSet, Base, BaseSequence,
                   Alphabet, dna_alphabet, dna_alphabet, protein_alphabet,
                   BpForm, DnaForm, RnaForm, ProteinForm, get_form)
