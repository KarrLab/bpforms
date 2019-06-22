""" Test Jupyter notebooks with examples

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-20
:Copyright: 2019, Karr Lab
:License: MIT
"""

import glob
import itertools
import json
import nbconvert.preprocessors
import nbformat
import os
import sys
import unittest


@unittest.skipIf(os.getenv('CIRCLECI', '0') in ['1', 'true'], 'Jupyter server not setup in CircleCI')
class ExamplesTestCase(unittest.TestCase):
    TIMEOUT = 600

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, 'examples')

    @classmethod
    def tearDownClass(cls):
        sys.path.remove('examples')

    def test_jupyter(self):
        for filename in itertools.chain(glob.glob('examples/*.ipynb'), glob.glob('examples/**/*.ipynb')):
            with open(filename) as file:
                version = json.load(file)['nbformat']
            with open(filename) as file:
                notebook = nbformat.read(file, as_version=version)
            execute_preprocessor = nbconvert.preprocessors.ExecutePreprocessor(timeout=self.TIMEOUT)
            execute_preprocessor.preprocess(notebook, {'metadata': {'path': 'examples/'}})

    def test_Bouhaddou_model(self):
        import bouhaddou_et_al_plos_comput_biol_2018
        if os.path.isfile(bouhaddou_et_al_plos_comput_biol_2018.OUT_FILENAME):
            os.remove(bouhaddou_et_al_plos_comput_biol_2018.OUT_FILENAME)
        bouhaddou_et_al_plos_comput_biol_2018.run()
        self.assertTrue(os.path.isfile(bouhaddou_et_al_plos_comput_biol_2018.OUT_FILENAME))

    def test_modomics(self):
        import modomics
        modomics.run()
        out_filename = os.path.join('examples', 'modomics.ssu-rrna.tsv')
        self.assertTrue(os.path.isfile(out_filename))
