""" Tests of building examples

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-14
:Copyright: 2019, Karr Lab
:License: MIT
"""

import shutil
import sys
import tempfile
import unittest


class BuildExamplesTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dirname = tempfile.mkdtemp()
        sys.path.append('docs')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dirname)
        sys.path.remove('docs')

    def test(self):
        import build_examples
        build_examples.build(self.dirname)
