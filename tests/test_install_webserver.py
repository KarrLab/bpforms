""" Tests of installing webserver

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-14
:Copyright: 2019, Karr Lab
:License: MIT
"""

import sys
import unittest


class InstallServerTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.append('.')

    @classmethod
    def tearDownClass(cls):
        sys.path.remove('.')

    def test(self):
        import install_webserver
        install_webserver.build(['dna'])
