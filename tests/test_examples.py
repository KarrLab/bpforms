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
import shutil
import sys
import tempfile
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

    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

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
        filename = os.path.join('examples', 'modomics.ssu-rrna.tsv')
        self.assertTrue(os.path.isfile(filename))

    def test_pro(self):
        import pro

        in_pkl_filename = os.path.join(self.dirname, 'in.pkl')
        out_pickle_filename = os.path.join(self.dirname, 'out.pkl')
        out_pickle_filename_2 = os.path.join(self.dirname, 'out.2.pkl')
        out_tsv_filename = os.path.join(self.dirname, 'out.tsv')
        out_fasta_filename = os.path.join(self.dirname, 'out.fasta')
        out_fig_filename = os.path.join(self.dirname, 'out.svg')

        pro.run(in_pkl_filename=in_pkl_filename, max_num_proteins=100,
                out_pickle_filename=out_pickle_filename,
                out_pickle_filename_2=out_pickle_filename_2,
                out_tsv_filename=out_tsv_filename, out_fasta_filename=out_fasta_filename,
                out_fig_filename=out_fig_filename,
                )

        self.assertTrue(os.path.isfile(out_tsv_filename))
