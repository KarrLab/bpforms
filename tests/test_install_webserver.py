""" Tests of installing webserver

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-14
:Copyright: 2019, Karr Lab
:License: MIT
"""
import requests
import sys
import unittest

try:
    response = requests.get('http://modomics.genesilico.pl/modifications/')
    modomics_available = response.status_code == 200 and response.elapsed.total_seconds() < 2.0
except requests.exceptions.ConnectionError:
    modomics_available = False


class InstallServerTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.append('.')

    @classmethod
    def tearDownClass(cls):
        sys.path.remove('.')

    def test(self):
        import install_webserver
        install_webserver.build(['dna'], build_modomics=modomics_available, pro_max_num_proteins=0)
