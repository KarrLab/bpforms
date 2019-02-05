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
            (['type'], dict(type=str, help='Type of biopolymer')),
            (['structure'], dict(type=str, help='Biopolymer structure')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        type = bpforms.util.get_form(args.type)
        try:
            form = type.from_str(args.structure)
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
            (['type'], dict(type=str, help='Type of biopolymer')),
            (['structure'], dict(type=str, help='Biopolymer structure')),
            (['--ph'], dict(default=None, type=float, help='pH')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        type = bpforms.util.get_form(args.type)
        try:
            form = type.from_str(args.structure)
        except Exception as error:
            raise SystemExit('Form is invalid: {}'.format(str(error)))
        if args.ph is not None:
            form.protonate(args.ph)
        print('Length: {}'.format(len(form)))
        print('Formula: {}'.format(form.get_formula()))
        print('Molecular weight: {}'.format(form.get_mol_wt()))
        print('Charge: {}'.format(form.get_charge()))


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
            (['ph'], dict(type=float, help='pH')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        type = bpforms.util.get_form(args.type)
        try:
            form = type.from_str(args.structure)
        except Exception as error:
            raise SystemExit('Form is invalid: {}'.format(str(error)))
        form.protonate(args.ph)
        print(str(form))


class App(cement.App):
    """ Command line application """
    class Meta:
        label = 'bpforms'
        base_controller = 'base'
        handlers = [
            BaseController,
            ValidateController,
            GetPropertiesController,
            ProtonateController,
        ]


def main():
    with App() as app:
        app.run()
