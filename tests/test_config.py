""" Tests of config module

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-03-11
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms import config
import unittest


class ConfigTestCase(unittest.TestCase):
    def test_get_config(self):
        config.get_config()
