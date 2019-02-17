""" Test of bpforms.util

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.alphabet import dna
from bpforms.alphabet import rna
from bpforms.alphabet import protein
from bpforms import core
from bpforms import util
import os
import shutil
import tempfile
import unittest


class UtilTestCase(unittest.TestCase):
    def setUp(self):
        os.rename(dna.filename, dna.filename + '.save')
        os.rename(rna.filename, rna.filename + '.save')
        os.rename(protein.filename, protein.filename + '.save')
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        os.rename(dna.filename + '.save', dna.filename)
        os.rename(rna.filename + '.save', rna.filename)
        os.rename(protein.filename + '.save', protein.filename)
        shutil.rmtree(self.tempdir)

    def test_get_alphabets(self):
        self.assertEqual(util.get_alphabets()['dna'], dna.dna_alphabet)

    def test_get_alphabet(self):
        self.assertEqual(util.get_alphabet('dna'), dna.dna_alphabet)
        self.assertEqual(util.get_alphabet('rna'), rna.rna_alphabet)
        self.assertEqual(util.get_alphabet('protein'), protein.protein_alphabet)
        with self.assertRaises(ValueError):
            util.get_alphabet('lipid')

    def test_get_form(self):
        self.assertEqual(util.get_form('dna'), dna.DnaForm)
        self.assertEqual(util.get_form('rna'), rna.RnaForm)
        self.assertEqual(util.get_form('protein'), protein.ProteinForm)
        self.assertEqual(util.get_form('canonical_dna'), dna.CanonicalDnaForm)
        self.assertEqual(util.get_form('canonical_rna'), rna.CanonicalRnaForm)
        self.assertEqual(util.get_form('canonical_protein'), protein.CanonicalProteinForm)
        self.assertEqual(util.get_form('curated_protein'), protein.CuratedProteinForm)
        with self.assertRaises(ValueError):
            util.get_form('lipid')

    def test_build_alphabets(self):
        self.assertFalse(os.path.isfile(dna.filename))
        self.assertFalse(os.path.isfile(rna.filename))
        self.assertFalse(os.path.isfile(protein.filename))

        util.build_alphabets(_max_monomers=3)

        self.assertTrue(os.path.isfile(dna.filename))
        self.assertTrue(os.path.isfile(rna.filename))
        self.assertTrue(os.path.isfile(protein.filename))

    def test_gen_html_viz_alphabet(self):
        filename = os.path.join(self.tempdir, 'alphabet.html')
        util.gen_html_viz_alphabet(dna.canonical_dna_alphabet, filename)
        self.assertTrue(os.path.isfile(filename))

        util.gen_html_viz_alphabet(core.Alphabet(monomers={'A': core.Monomer()}), filename)
