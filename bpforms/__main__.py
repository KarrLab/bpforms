""" bpforms command line interface

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

import bpforms
import bpforms.util
import cement


class BaseController(cement.Controller):
    """ Base controller for command line application """

    class Meta:
        label = 'base'
        description = "bpforms"
        arguments = [
            (['-v', '--version'], dict(action='version', version=bpforms.__version__)),
        ]

    @cement.ex(hide=True)
    def _default(self):
        raise SystemExit(self._parser.print_help())


class ValidateController(cement.Controller):
    """ Validate a biopolymer form """

    class Meta:
        label = 'validate'
        description = 'Validate a biopolymer form'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['alphabet'], dict(type=str, help='Biopolymer alphabet')),
            (['monomer_seq'], dict(type=str, help='Monomer sequence')),
            (['--circular'], dict(action='store_true', default=False, help='Biopolymer circularity')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        type = bpforms.util.get_form(args.alphabet)
        try:
            form = type(circular=args.circular).from_str(args.monomer_seq)
            print('Form is valid')
        except Exception as error:
            raise SystemExit('Form is invalid: {}'.format(str(error)))


class GetPropertiesController(cement.Controller):
    """ Calculate physical properties such as length, chemical formula, molecular weight, and charge """

    class Meta:
        label = 'get-properties'
        description = 'Calculate physical properties such as length, chemical formula, molecular weight, and charge'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['alphabet'], dict(type=str, help='Biopolymer alphabet')),
            (['monomer_seq'], dict(type=str, help='Monomer sequence')),
            (['--circular'], dict(action='store_true', default=False, help='Biopolymer circularity')),
            (['--ph'], dict(default=None, type=float,
                            help='pH at which calculate major protonation state of each monomer')),
            (['--not-major-tautomer'], dict(action='store_true', default=False,
                                        help='If set, do not calculate the major tautomer')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        type = bpforms.util.get_form(args.alphabet)
        try:
            form = type(circular=args.circular).from_str(args.monomer_seq)
        except Exception as error:
            raise SystemExit('Form is invalid: {}'.format(str(error)))
        if args.ph is not None:
            form.get_major_micro_species(args.ph, major_tautomer=not args.not_major_tautomer)
        print('Length: {}'.format(len(form)))

        try:
            structure = form.export('smiles')
        except Exception:
            structure = None
        print('Structure: {}'.format(structure))
        print('Formula: {}'.format(form.get_formula()))
        print('Molecular weight: {}'.format(form.get_mol_wt()))
        print('Charge: {}'.format(form.get_charge()))


class GetMajorMicroSpeciesController(cement.Controller):
    """ Calculate the major protonation and tautomerization """

    class Meta:
        label = 'get-major-micro-species'
        description = 'Calculate the major protonation and tautomerization state of a biopolymer form to a specific pH'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['alphabet'], dict(type=str, help='Biopolymer alphabet')),
            (['monomer_seq'], dict(type=str, help='Monomer sequence')),
            (['--circular'], dict(action='store_true', default=False, help='Biopolymer circularity')),
            (['ph'], dict(type=float, help='pH')),
            (['--not-major-tautomer'], dict(action='store_true', default=False,
                                        help='If set, do not calculate the major tautomer')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        type = bpforms.util.get_form(args.alphabet)
        try:
            form = type(circular=args.circular).from_str(args.monomer_seq)
        except Exception as error:
            raise SystemExit('Form is invalid: {}'.format(str(error)))
        form.get_major_micro_species(args.ph, major_tautomer=not args.not_major_tautomer)
        print(str(form))


class BuildAlphabetsController(cement.Controller):
    """ Build DNA, RNA, and protein alphabets from DNAmod, MODOMICS, and RESID """

    class Meta:
        label = 'build-alphabets'
        description = 'Build DNA, RNA, and protein alphabets from DNAmod, MODOMICS, and RESID'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['--ph'], dict(type=float, default=7.4,
                            help='pH at which calculate major protonation state of each monomer')),
            (['--not-major-tautomer'], dict(action='store_true', default=False,
                                            help='If set, do not calculate the major tautomer')),
            (['--max-monomers'], dict(type=float, default=float('inf'),
                                      help='Maximum number of monomers to build. Used for testing')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        bpforms.util.build_alphabets(ph=args.ph, major_tautomer=not args.not_major_tautomer, _max_monomers=args.max_monomers)
        print('Alphabets successfully built')


class VizAlphabetController(cement.Controller):
    """ Visualize an alphabet """

    class Meta:
        label = 'viz-alphabet'
        description = 'Visualize an alphabet'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['alphabet'], dict(type=str, help='Biopolymer alphabet')),
            (['path'], dict(type=str, help='Path to save visualization of alphabet')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        alphabet = bpforms.util.get_alphabet(args.alphabet)
        bpforms.util.gen_html_viz_alphabet(alphabet, args.path)
        print('Visualization saved to {}'.format(args.path))


class App(cement.App):
    """ Command line application """
    class Meta:
        label = 'bpforms'
        base_controller = 'base'
        handlers = [
            BaseController,
            ValidateController,
            GetPropertiesController,
            GetMajorMicroSpeciesController,
            BuildAlphabetsController,
            VizAlphabetController,
        ]


def main():
    with App() as app:
        app.run()
