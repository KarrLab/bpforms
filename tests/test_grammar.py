""" Test that grammar is compatible with EBNF

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-07-18
:Copyright: 2019, Karr Lab
:License: MIT
"""

import pkg_resources
import plyplus
import unittest


class GrammarTestCase(unittest.TestCase):
    def test(self):
        filename = pkg_resources.resource_filename('bpforms', 'grammar.ebnf')

        with open(filename, 'r') as file:
            plyplus.Grammar(file)
