""" bpforms command line interface

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_utils.util.chem import OpenBabelUtils
import bpforms
import bpforms.util
import cement


class BaseController(cement.Controller):
    """ Base controller for command line application """

    class Meta:
        label = 'base'
        description = "bpforms"
        help = "bpforms"
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
        help = 'Validate a biopolymer form'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['alphabet'], dict(type=str, help='Biopolymer alphabet')),
            (['seq'], dict(type=str, help='Sequence of monomeric forms')),
            (['--circular'], dict(action='store_true', default=False, help='Biopolymer circularity')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        type = bpforms.util.get_form(args.alphabet)

        try:
            form = type().from_str(args.seq)
        except Exception as error:
            raise SystemExit('Form is invalid: {}'.format(str(error)))
        form.circular = args.circular

        errors = form.validate()
        if errors:
            raise SystemExit('Form is invalid:\n  {}'.format('\n  '.join(errors)))

        print('Form is valid')


class GetPropertiesController(cement.Controller):
    """ Calculate physical properties such as length, chemical formula, molecular weight, and charge """

    class Meta:
        label = 'get-properties'
        description = 'Calculate physical properties such as length, chemical formula, molecular weight, and charge'
        help = 'Calculate physical properties such as length, chemical formula, molecular weight, and charge'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['alphabet'], dict(type=str, help='Biopolymer alphabet')),
            (['seq'], dict(type=str, help='Sequence of monomeric forms')),
            (['--circular'], dict(action='store_true', default=False, help='Biopolymer circularity')),
            (['--ph'], dict(default=None, type=float,
                            help='pH at which calculate major protonation state of each monomeric form')),
            (['--major-tautomer'], dict(action='store_true', default=False,
                                        help='If set, calculate the major tautomer')),
            (['--dearomatize'], dict(action='store_true', default=False,
                                     help='If set, dearomatize molecule')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        type = bpforms.util.get_form(args.alphabet)

        try:
            form = type().from_str(args.seq)
        except Exception as error:
            raise SystemExit('Form is invalid: {}'.format(str(error)))
        form.circular = args.circular

        errors = form.validate()
        if errors:
            raise SystemExit('Form is invalid:\n  {}'.format('\n  '.join(errors)))

        smiles = None
        formula = None
        mol_wt = None
        charge = None
        try:
            if args.ph is None:
                formula = form.get_formula()
                mol_wt = form.get_mol_wt()
                charge = form.get_charge()
                structure = form.get_structure()[0]
            else:
                structure = form.get_major_micro_species(args.ph, major_tautomer=args.major_tautomer, dearomatize=args.dearomatize)
                if structure is not None:
                    formula = OpenBabelUtils.get_formula(structure)
                    mol_wt = formula.get_molecular_weight()
                    charge = structure.GetTotalCharge()

            if structure is None:
                smiles = None
            else:
                smiles = OpenBabelUtils.export(structure, 'smiles')
        except Exception:
            pass

        print('Length: {}'.format(len(form)))
        print('Structure: {}'.format(smiles))
        print('Formula: {}'.format(formula))
        print('Molecular weight: {}'.format(mol_wt))
        print('Charge: {}'.format(charge))


class GetMajorMicroSpeciesController(cement.Controller):
    """ Calculate the major protonation and tautomerization """

    class Meta:
        label = 'get-major-micro-species'
        description = 'Calculate the major protonation and tautomerization state of a biopolymer form to a specific pH'
        help = 'Calculate the major protonation and tautomerization state of a biopolymer form to a specific pH'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['alphabet'], dict(type=str, help='Biopolymer alphabet')),
            (['seq'], dict(type=str, help='Sequence of monomeric forms')),
            (['--circular'], dict(action='store_true', default=False, help='Biopolymer circularity')),
            (['ph'], dict(type=float, help='pH')),
            (['--major-tautomer'], dict(action='store_true', default=False,
                                        help='If set, calculate the major tautomer')),
            (['--dearomatize'], dict(action='store_true', default=False,
                                     help='If set, dearomatize molecule')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        type = bpforms.util.get_form(args.alphabet)

        try:
            form = type().from_str(args.seq)
        except Exception as error:
            raise SystemExit('Form is invalid: {}'.format(str(error)))
        form.circular = args.circular

        errors = form.validate()
        if errors:
            raise SystemExit('Form is invalid:\n  {}'.format('\n  '.join(errors)))

        structure = form.get_major_micro_species(args.ph, major_tautomer=args.major_tautomer, dearomatize=args.dearomatize)
        print(OpenBabelUtils.export(structure, 'smiles'))


class BuildAlphabetsController(cement.Controller):
    """ Build DNA, RNA, and protein alphabets from DNAmod, MODOMICS, the PDB Chemical Component Dictionary, RESID,
    and the RNA Modification Database """

    class Meta:
        label = 'build-alphabets'
        description = ('Build DNA, RNA, and protein alphabets from DNAmod, MODOMICS, '
                       'the PDB Chemical Component Dictionary, RESID, and the RNA Modification Database')
        help = ('Build DNA, RNA, and protein alphabets from DNAmod, MODOMICS, '
                       'the PDB Chemical Component Dictionary, RESID, and the RNA Modification Database')
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['--ph'], dict(type=float, default=7.4,
                            help='pH at which calculate major protonation state of each monomeric form')),
            (['--major-tautomer'], dict(action='store_true', default=False,
                                        help='If set, calculate the major tautomer')),
            (['--dearomatize'], dict(action='store_true', default=False,
                                     help='If set, dearomatize molecule')),
            (['--max-monomers'], dict(type=float, default=float('inf'),
                                      help='Maximum number of monomeric forms to build. Used for testing')),
            (['--alphabet'], dict(type=str, default=None, dest='alphabets', action='append',
                                  help='Id of alphabet to build. Defualt: build all alphabets')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        bpforms.util.build_alphabets(ph=args.ph, major_tautomer=args.major_tautomer, dearomatize=args.dearomatize,
                                     _max_monomers=args.max_monomers, alphabets=args.alphabets)
        print('Alphabets successfully built')


class VizAlphabetController(cement.Controller):
    """ Visualize an alphabet """

    class Meta:
        label = 'viz-alphabet'
        description = 'Visualize an alphabet'
        help = 'Visualize an alphabet'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['alphabet'], dict(type=str, help='Biopolymer alphabet')),
            (['path'], dict(type=str, help='Path to save visualization of alphabet')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        bpform_type = bpforms.util.get_form(args.alphabet)
        bpforms.util.gen_html_viz_alphabet(bpform_type, args.path)
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
