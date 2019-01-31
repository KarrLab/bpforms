""" bpforms command line interface

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

import cement
import bpforms
import bpforms.core


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
        self._parser.print_help()


class ValidateController(cement.Controller):
    """ Validate a biopolymer form """

    class Meta:
        label = 'validate'
        description = 'Validate a biopolymer form'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['type'], dict(type=str, help='Type of biopolymer')),
            (['structure'], dict(type=str, help='Biopolymer structure')),
            (['--ph'], dict(default=7., type=float, help='pH')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        type = bpforms.core.get_form(args.type)
        try:
            form = type(args.structure, ph=args.ph)
            print('Form is valid')
        except Exception as error:
            print('Form is invalid: '.format(str(error)))


class CalcFormulaController(cement.Controller):
    """ Calculate chemical formula """

    class Meta:
        label = 'get-formula'
        description = 'Calculate chemical formula'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['type'], dict(type=str, help='Type of biopolymer')),
            (['structure'], dict(type=str, help='Biopolymer structure')),
            (['--ph'], dict(default=7., type=float, help='pH')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        type = bpforms.core.get_form(args.type)
        form = type(args.structure, ph=args.ph)
        print(form.get_formula())


class ProtonateController(cement.Controller):
    """ Protonate a biopolymer form to a specific pH """

    class Meta:
        label = 'protonate'
        description = 'Protonate a biopolymer form to a specific pH'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['type'], dict(type=str, help='Type of biopolymer')),
            (['structure'], dict(type=str, help='Biopolymer structure')),
            (['--ph'], dict(default=7., type=float, help='pH')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        type = bpforms.core.get_form(args.type)
        form = type(args.structure, ph=args.ph)
        form.protonate()
        print(form.structure)


class VisualizeController(cement.Controller):
    """ Generate a graphical depiction of the chemical structure of a biopolymer form """

    class Meta:
        label = 'visualize'
        description = 'Generate a graphical depiction of the chemical structure of a biopolymer form'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['type'], dict(type=str, help='Type of biopolymer')),
            (['structure'], dict(type=str, help='Biopolymer structure')),
            (['filename'], dict(type=str, help='Path to save graphical depiction of the biopolymer form')),
            (['--ph'], dict(default=7., type=float, help='pH')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        type = bpforms.core.get_form(args.type)
        form = type(args.structure, ph=args.ph)
        form.gen_visualization(args.filename)


class App(cement.App):
    """ Command line application """
    class Meta:
        label = 'bpforms'
        base_controller = 'base'
        handlers = [
            BaseController,
            ValidateController,
            CalcFormulaController,
            ProtonateController,
            VisualizeController,
        ]


def main():
    with App() as app:
        app.run()
