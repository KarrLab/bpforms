""" Test of calculating structures of biopolymers of various sizes

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-03-11
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms import __main__
from bpforms import core
from bpforms import rest
from bpforms.alphabet import dna
from bpforms.alphabet import protein
from bpforms.alphabet import rna
from wc_utils.util.chem import EmpiricalFormula, OpenBabelUtils
import hurry.filesize
import openbabel
import os
import psutil
import time
import unittest


class LargeBpFormsTestCase(unittest.TestCase):
    def test_dna(self):
        self.verify_large_polymers(dna.DnaForm, 'ACGT')

    def test_rna(self):
        self.verify_large_polymers(rna.RnaForm, 'ACGU')

    def test_protein(self):
        self.verify_large_polymers(protein.ProteinForm, 'ACDE')

    def test_get_structure(self):
        form = protein.ProteinForm().from_str('ARCGY' * 100)
        structure, _ = form.get_structure()
        self.assertIsInstance(structure, openbabel.OBMol)
        cml = OpenBabelUtils.export(structure, 'cml')
        self.assertTrue(cml.startswith('<molecule'))

    def verify_large_polymers(self, form_type, alphabet):
        # test Python API
        form = form_type()
        for i_trial in range(6):
            self.verify_large_polymer(form, alphabet, i_trial)

        # test CLI
        with __main__.App(argv=['get-properties', form.alphabet.id, alphabet * pow(2, 5)]) as app:
            app.run()
        with __main__.App(argv=['get-properties', form.alphabet.id, alphabet * pow(2, 5), '--ph', '7.4']) as app:
            app.run()

        # test REST API
        client = rest.app.test_client()
        rv = client.post('/api/bpform/', json=dict(alphabet=form.alphabet.id, seq=alphabet * pow(2, 5),
                                                   ph=7.4, major_tautomer=True))
        self.assertEqual(rv.status_code, 200)

    def verify_large_polymer(self, form, alphabet, i_trial):
        start = time.time()
        process = psutil.Process(os.getpid())
        uss0 = process.memory_full_info().uss

        form.from_str(alphabet * pow(2, i_trial))
        length = 4 * pow(2, i_trial)
        self.assertEqual(len(form), length)
        formula = form.get_formula()
        charge = form.get_charge()
        end1 = time.time()
        uss1 = process.memory_full_info().uss

        if length <= 50:
            structure = form.get_structure()[0]
            self.assertEqual(OpenBabelUtils.get_formula(structure), formula)
            self.assertEqual(structure.GetTotalCharge(), charge)
        else:
            structure = None
        end2 = time.time()
        uss2 = process.memory_full_info().uss

        if length <= 20:
            form.get_major_micro_species(7.4, major_tautomer=False)
        end3 = time.time()
        uss3 = process.memory_full_info().uss

        if length <= 5:
            form.get_major_micro_species(7.4, major_tautomer=True)
        end4 = time.time()
        uss4 = process.memory_full_info().uss

        if structure is not None:
            OpenBabelUtils.export(structure, 'smiles')
        end5 = time.time()
        uss5 = process.memory_full_info().uss

        print(('Calculating polymer of length {} took {:.3f} s'
               '\n  Parsing, length, formula, charge: {:.3f} s, {}'
               '\n  Structure: {:.3f} s, {}'
               '\n  Major microspecies: {:.3f} s, {}'
               '\n  Major tautomer: {:.3f} s, {}'
               '\n  Canonical SMILES: {:.3f} s, {}'
               ).format(
            length, end3 - start,
            end1 - start, hurry.filesize.size(uss1 - uss0),
            end2 - end1, hurry.filesize.size(uss2 - uss1),
            end3 - end2, hurry.filesize.size(uss3 - uss2),
            end4 - end3, hurry.filesize.size(uss4 - uss3),
            end5 - end4, hurry.filesize.size(uss5 - uss4),
        ))
