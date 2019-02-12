""" Tests of bpforms command line interface (bpforms.__main__)

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms import __main__
import bpforms
import bpforms.alphabet.dna
import bpforms.alphabet.rna
import bpforms.alphabet.protein
import capturer
import mock
import os
import shutil
import tempfile
import unittest

ala_inchi = 'InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1'
ala_inchi_ph_14 = 'InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/p-1/t2-/m0/s1'


class CliTestCase(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_cli(self):
        with mock.patch('sys.argv', ['bpforms', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegex(context.Exception, 'usage: bpforms')

    def test_help(self):
        with self.assertRaises(SystemExit):
            with __main__.App(argv=[]) as app:
                app.run()

        with self.assertRaises(SystemExit):
            with __main__.App(argv=['--help']) as app:
                app.run()

    def test_version(self):
        with __main__.App(argv=['-v']) as app:
            with capturer.CaptureOutput(merged=False, relay=False) as captured:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(captured.stdout.get_text(), bpforms.__version__)
                self.assertEqual(captured.stderr.get_text(), '')

        with __main__.App(argv=['--version']) as app:
            with capturer.CaptureOutput(merged=False, relay=False) as captured:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(captured.stdout.get_text(), bpforms.__version__)
                self.assertEqual(captured.stderr.get_text(), '')

    def test_validate(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['validate', 'canonical_dna', 'ACGT']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), 'Form is valid')
                self.assertEqual(captured.stderr.get_text(), '')

        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['validate', 'canonical_dna', 'ACG[id: "ala" | structure: {}]T'.format(ala_inchi)]) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), 'Form is valid')
                self.assertEqual(captured.stderr.get_text(), '')

        with self.assertRaisesRegex(SystemExit, '^Form is invalid'):
            with __main__.App(argv=['validate', 'canonical_dna', 'ACGT[']) as app:
                # run app
                app.run()

    def test_get_properties(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['get-properties', 'canonical_dna', 'ACGT']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                text = captured.stdout.get_text()
                self.assertIn('Length: 4', text)
                self.assertIn('Formula: C39', text)
                self.assertIn('Molecular weight: ', text)
                self.assertIn('Charge: -', text)
                self.assertEqual(captured.stderr.get_text(), '')

        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            monomer_seq = ''.join([
                '[structure: {}]'.format(bpforms.alphabet.dna.dna_alphabet.monomers.A.get_inchi()),
                '[structure: {}]'.format(bpforms.alphabet.dna.dna_alphabet.monomers.C.get_inchi()),
                '[structure: {}]'.format(bpforms.alphabet.dna.dna_alphabet.monomers.G.get_inchi()),
                '[structure: {}]'.format(bpforms.alphabet.dna.dna_alphabet.monomers.T.get_inchi()),
            ])
            with __main__.App(argv=['get-properties', 'canonical_dna', monomer_seq, '--ph', '7.0']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                text = captured.stdout.get_text()
                self.assertIn('Length: 4', text)
                self.assertIn('Formula: C39', text)
                self.assertIn('Molecular weight: ', text)
                self.assertIn('Charge: -', text)

        with self.assertRaises(SystemExit):
            with __main__.App(argv=['get-properties', 'canonical_dna', 'ACGT[']) as app:
                # run app
                app.run()

    def test_protonate(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['protonate', 'canonical_dna', '[id: "ala" | structure: {0}][id: "ala" | structure: {0}]'.format(
                    ala_inchi), '14.']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), '[id: "ala" | structure: {0}][id: "ala" | structure: {0}]'.format(
                    ala_inchi_ph_14))

        with self.assertRaises(SystemExit):
            with __main__.App(argv=['protonate', 'canonical_dna', 'ACGT[', '7.']) as app:
                # run app
                app.run()


class BuildAlphabetsCliTestCase(unittest.TestCase):
    def setUp(self):
        os.rename(bpforms.alphabet.dna.filename, bpforms.alphabet.dna.filename + '.save')
        os.rename(bpforms.alphabet.rna.filename, bpforms.alphabet.rna.filename + '.save')
        os.rename(bpforms.alphabet.protein.filename, bpforms.alphabet.protein.filename + '.save')

    def tearDown(self):
        os.rename(bpforms.alphabet.dna.filename + '.save', bpforms.alphabet.dna.filename)
        os.rename(bpforms.alphabet.rna.filename + '.save', bpforms.alphabet.rna.filename)
        os.rename(bpforms.alphabet.protein.filename + '.save', bpforms.alphabet.protein.filename)

    def test_build_alphabets(self):
        self.assertFalse(os.path.isfile(bpforms.alphabet.dna.filename))
        self.assertFalse(os.path.isfile(bpforms.alphabet.rna.filename))
        self.assertFalse(os.path.isfile(bpforms.alphabet.protein.filename))

        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['build-alphabets', '--max-monomers', '3']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertIn('Alphabets successfully built', captured.stdout.get_text())

        self.assertTrue(os.path.isfile(bpforms.alphabet.dna.filename))
        self.assertTrue(os.path.isfile(bpforms.alphabet.rna.filename))
        self.assertTrue(os.path.isfile(bpforms.alphabet.protein.filename))
