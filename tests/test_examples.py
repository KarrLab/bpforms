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
import unittest


class ExamplesTestCase(unittest.TestCase):
    TIMEOUT = 600

    def test(self):
        for filename in itertools.chain(glob.glob('examples/*.ipynb'), glob.glob('examples/**/*.ipynb')):
            with open(filename) as file:
                version = json.load(file)['nbformat']
            with open(filename) as file:
                notebook = nbformat.read(file, as_version=version)
            execute_preprocessor = nbconvert.preprocessors.ExecutePreprocessor(timeout=self.TIMEOUT)
            execute_preprocessor.preprocess(notebook, {'metadata': {'path': 'examples/'}})
