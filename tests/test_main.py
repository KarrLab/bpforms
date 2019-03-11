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

dI_smiles = 'O=C1NC=NC2=C1N=CN2'
dI_smiles_ph_14 = 'O=c1[n-]cnc2c1nc[n-]2'


class CliTestCase(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempdir)
        bpforms.alphabet.dna.canonical_dna_alphabet.from_yaml(bpforms.alphabet.dna.canonical_filename)

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
            with __main__.App(argv=['validate', 'canonical_dna', 'ACG[id: "dI" | structure: "{}"]T'.format(dI_smiles)]) as app:
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
                self.assertNotIn('Structure: None', text)
                self.assertIn('Formula: C39', text)
                self.assertIn('Molecular weight: ', text)
                self.assertIn('Charge: -', text)
                self.assertEqual(captured.stderr.get_text(), '')

        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['get-properties', 'canonical_dna', 'ACGT', '--ph', '7.0']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                text = captured.stdout.get_text()
                self.assertIn('Length: 4', text)
                self.assertIn('Structure: ' + (
                    'Cc1cn(C2CC(O)C(COP(=O)([O-])OC3CC(OC3COP(=O)([O-])OC3CC(OC3COP(=O)([O-])'
                    'OC3CC(OC3COP(=O)([O-])[O-])n3cnc4c(N)ncnc34)n3ccc(N)nc3=O)n3cnc4c3nc(N)[nH]'
                    'c4=O)O2)c(=O)[nH]c1=O'), text)
                self.assertIn('Formula: C39', text)
                self.assertIn('Molecular weight: 1248.772047992', text)
                self.assertIn('Charge: -5', text)

        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['get-properties', 'canonical_dna', 'ACGT', '--ph', '7.0']) as app:
                # run app
                with mock.patch.object(bpforms.BpForm, 'export', side_effect=Exception('error')):
                    app.run()

                # test that the CLI produced the correct output
                text = captured.stdout.get_text()
                self.assertIn('Length: 4', text)
                self.assertIn('Structure: None', text)
                self.assertIn('Formula: None', text)
                self.assertIn('Molecular weight: None', text)
                self.assertIn('Charge: None', text)

        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['get-properties', 'canonical_dna', 'ACGT', '--circular']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                text = captured.stdout.get_text()
                self.assertIn('Length: 4', text)
                self.assertNotIn('Structure: None', text)
                self.assertIn('Formula: C39', text)
                self.assertIn('Molecular weight: ', text)
                self.assertIn('Charge: -', text)
                self.assertEqual(captured.stderr.get_text(), '')

        with self.assertRaises(SystemExit):
            with __main__.App(argv=['get-properties', 'canonical_dna', 'ACGT[']) as app:
                # run app
                app.run()

    def test_get_major_micro_species(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['get-major-micro-species', 'canonical_dna',
                                    '[id: "dI" | structure: "{0}"][id: "dI" | structure: "{0}"]'.format(
                    dI_smiles), '14.']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), 'OC1CCOC1COP(=O)([O-])OC1CCOC1COP(=O)([O-])[O-]')

        with self.assertRaises(SystemExit):
            with __main__.App(argv=['get-major-micro-species', 'canonical_dna', 'ACGT[', '7.']) as app:
                # run app
                app.run()

    def test_viz_alphabet(self):
        path = os.path.join(self.tempdir, 'alphabet.html')

        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['viz-alphabet', 'canonical_dna', path]) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), 'Visualization saved to {}'.format(path))
        self.assertTrue(os.path.isfile(path))


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
