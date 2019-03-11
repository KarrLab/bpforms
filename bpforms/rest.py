""" REST JSON API

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_utils.util.chem import OpenBabelUtils
import bpforms
import bpforms.util
import flask
import flask_restplus
import flask_restplus.errors
import flask_restplus.fields
import math
import warnings

app = flask.Flask(__name__)


class PrefixMiddleware(object):
    def __init__(self, app, prefix=''):
        self.app = app
        self.prefix = prefix

    def __call__(self, environ, start_response):
        if environ['PATH_INFO'].startswith(self.prefix):
            environ['PATH_INFO'] = environ['PATH_INFO'][len(self.prefix):]
            environ['SCRIPT_NAME'] = self.prefix
            return self.app(environ, start_response)
        else:
            start_response('404', [('Content-Type', 'text/plain')])
            return ["This url does not belong to the app.".encode()]


app.wsgi_app = PrefixMiddleware(app.wsgi_app, prefix='/api')

api = flask_restplus.Api(app,
                         title='Bpforms JSON REST API',
                         description='JSON REST API for calculating properties of biopolymer forms',
                         contact='karr@mssm.edu',
                         version=bpforms.__version__,
                         license='MIT',
                         license_url='https://github.com/KarrLab/bpforms/blob/master/LICENSE',
                         doc='/')

bpform_ns = flask_restplus.Namespace('bpform', description='Calculate properties of biopolymer forms')
api.add_namespace(bpform_ns)

bpforms_model = bpform_ns.model('BpForm', {
    'alphabet': flask_restplus.fields.String(enum=list(bpforms.util.get_alphabets().keys()), required=True,
                                             title='Alphabet',
                                             description='Id of the alphabet of the biopolymer form'),
    'monomer_seq': flask_restplus.fields.String(required=True,
                                                title='Monomer sequence',
                                                description='Sequence of monomers of the biopolymer form',
                                                example='AA'),
    'circular': flask_restplus.fields.Boolean(default=False,
                                              required=False,
                                              title='Circularity',
                                              description='Circularity of the biopolymer form',
                                              example='False'),
    'ph': flask_restplus.fields.Float(default=float('NaN'), min=0., max=14., required=False,
                                      title='pH',
                                      description='pH at which to calculate the major microspecies of the biopolymer form',
                                      example=7.4),
    'major_tautomer': flask_restplus.fields.Boolean(default=False, required=False,
                                                    title='Calculate major tautomer',
                                                    description='If true, calculate the major tautomer',
                                                    example=True),
})


@bpform_ns.route("/")
class Bpform(flask_restplus.Resource):
    """ Optionally, calculate the major protonation and tautomerization form a biopolymer form and calculate its properties """

    @bpform_ns.doc('Optionally, calculate the major protonation and tautomerization form a biopolymer form and calculate its properties')
    @bpform_ns.expect(bpforms_model, validate=True)
    def post(self):
        """ Optionally, calculate the major protonation and tautomerization form a biopolymer form and calculate its properties """
        """
        Returns:
            :obj:`dict`
        """

        args = bpform_ns.payload
        alphabet = args['alphabet']
        monomer_seq = args['monomer_seq']
        circular = args.get('circular', False)
        ph = args.get('ph', float('NaN'))
        major_tautomer = args.get('major_tautomer', False)

        form_cls = bpforms.util.get_form(alphabet)

        try:
            form = form_cls(circular=circular).from_str(monomer_seq)
        except Exception as error:
            flask_restplus.abort(400, 'Unable to parse monomer sequence', errors={'monomer_seq': str(error)})

        smiles = None
        formula = None
        mol_wt = None
        charge = None
        with warnings.catch_warnings(record=True) as recorded_warnings:
            warnings.simplefilter('once', bpforms.BpFormsWarning)

            try:
                if math.isnan(ph):
                    formula = dict(form.get_formula())
                    mol_wt = form.get_mol_wt()
                    charge = form.get_charge()
                    structure = form.get_structure()
                else:
                    structure = form.get_major_micro_species(ph, major_tautomer=major_tautomer)
                    if structure is not None:
                        formula = OpenBabelUtils.get_formula(structure)
                        mol_wt = formula.get_molecular_weight()
                        formula = dict(formula)
                        charge = structure.GetTotalCharge()
                if structure is None:
                    smiles = None
                else:
                    smiles = OpenBabelUtils.export(structure, 'smiles')
            except Exception:
                pass

            if recorded_warnings:
                warning_message = ' '.join(str(recorded_warning.message) for recorded_warning in recorded_warnings)
            else:
                warning_message = None

        return {
            'alphabet': alphabet,
            'monomer_seq': str(form),
            'length': len(form),
            'structure': smiles,
            'formula': formula,
            'mol_wt': mol_wt,
            'charge': charge,
            'warnings': warning_message,
        }


alphabet_ns = flask_restplus.Namespace('alphabet', description='List alphabets and get their details')
api.add_namespace(alphabet_ns)


@alphabet_ns.route("/")
@alphabet_ns.doc(params={})
class AlphabetsResource(flask_restplus.Resource):
    """ Get list of alphabets """

    def get(self):
        """ Get a list of available alphabets """
        """
        Returns:
            :obj:`dict`: dictionary that maps that ids of available alphabets to dictionaries with
                properties of the alphabets
        """
        rv = {}
        for id, alphabet in bpforms.util.get_alphabets().items():
            rv[id] = {
                'id': alphabet.id,
                'name': alphabet.name,
                'description': alphabet.description,
            }
        return rv


@alphabet_ns.route("/<string:id>/")
@alphabet_ns.doc(params={
    'id': 'Id of the alphabet of the biopolymer form (e.g. "dna", "rna", or "protein")',
})
class AlpabetResource(flask_restplus.Resource):
    """ Get alphabets """

    def get(self, id):
        """ Get an alphabet """
        """
        Returns:
            :obj:`dict`: dictionary representation of an alphabet
        """
        try:
            alphabet_obj = bpforms.util.get_alphabet(id)
        except ValueError as error:
            flask_restplus.abort(400, 'Invalid alphabet "{}"'.format(id))

        return alphabet_obj.to_dict()
