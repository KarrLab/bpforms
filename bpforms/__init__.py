import pkg_resources

# read version
with open(pkg_resources.resource_filename('bpforms', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()

from .core import (Identifier, IdentifierSet, SynonymSet, Monomer, MonomerSequence, MonomerDict,
                   Backbone, BondBase, Bond, BondOrder, BondStereo, OntoBond, BondSet, Atom, AtomList,
                   Alphabet, AlphabetBuilder, BpForm, BpFormFeature, BpFormFeatureSet, BpFormsWarning,
                   Nick, NickSet)
from .alphabet.dna import dna_alphabet, canonical_dna_alphabet, DnaForm, CanonicalDnaForm
from .alphabet.rna import rna_alphabet, canonical_rna_alphabet, RnaForm, CanonicalRnaForm
from .alphabet.protein import protein_alphabet, canonical_protein_alphabet, ProteinForm, CanonicalProteinForm
from .util import get_form
from .xlink.core import crosslinks_onto
