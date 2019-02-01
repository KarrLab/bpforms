""" Tests of bpforms command line interface (bpforms.__main__)

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms import __main__
import bpforms
import capturer
import mock
import os
import shutil
import tempfile
import unittest


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

    @unittest.skip('Code not implemented yet')
    def test_validate(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['validate', 'rna', 'ACGU']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), 'Form is valid')
                self.assertEqual(captured.stderr.get_text(), '')

        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['validate', 'rna', 'ACGU[']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), '')
                self.assertRegex(captured.stderr.get_text(), '^Form is invalid')

    @unittest.skip('Code not implemented yet')
    def test_get_properties(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['get-formula', 'rna', 'ACGU', '--ph', '7.0']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), 'C H O N P')
                self.assertEqual(captured.stderr.get_text(), '')

    @unittest.skip('Code not implemented yet')
    def test_protonate(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['protonate', 'rna', 'ACGU[structure: InChI=1/]']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), 'ACGU[structure: InChI=1/]')
                self.assertEqual(captured.stderr.get_text(), '')
