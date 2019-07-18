""" REST JSON API

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

from .config import get_config
from wc_utils.util.chem import OpenBabelUtils
import bpforms
import bpforms.core
import bpforms.util
import flask
import flask_restplus
import flask_restplus.errors
import flask_restplus.fields
import math
import os
import pkg_resources
import urllib.parse
import warnings

config = get_config()['bpforms']['rest']

# setup app
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
    'seq': flask_restplus.fields.String(required=True,
                                        title='Sequence of monomeric forms',
                                        description='Sequence of monomeric forms of the biopolymer form',
                                        example='AA'),
    'circular': flask_restplus.fields.Boolean(default=False,
                                              required=False,
                                              title='Circularity',
                                              description='Circularity of the biopolymer form',
                                              example=False),
    'ph': flask_restplus.fields.Float(default=float('NaN'), min=0., max=14., required=False,
                                      title='pH',
                                      description='pH at which to calculate the major microspecies of the biopolymer form',
                                      example=7.4),
    'major_tautomer': flask_restplus.fields.Boolean(default=False, required=False,
                                                    title='Calculate major tautomer',
                                                    description='If true, calculate the major tautomer',
                                                    example=True),
    'dearomatize': flask_restplus.fields.Boolean(default=False, required=False,
                                                 title='Dearomatize the molecule',
                                                 description='If true, calculate dearomatize the molecule',
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
        seq = args['seq']
        circular = args.get('circular', False)
        ph = args.get('ph', float('NaN'))
        major_tautomer = args.get('major_tautomer', False)
        dearomatize = args.get('dearomatize', False)

        form_cls = bpforms.util.get_form(alphabet)

        try:
            form = form_cls().from_str(seq)
        except Exception as error:
            flask_restplus.abort(400, 'Unable to parse sequence of monomeric forms', errors={'seq': str(error)})
        form.circular = circular

        errors = form.validate()
        if errors:
            flask_restplus.abort(400, 'Form is invalid', errors={'seq': '. '.join(errors)})

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

                    if len(form.seq) <= config['max_len_get_structure']:
                        structure = form.get_structure()[0]
                    else:
                        structure = None
                        warnings.warn('Structure calculations are limited to forms with length <= {}'.format(
                            config['max_len_get_structure']), bpforms.BpFormsWarning)

                else:
                    if major_tautomer and len(form.seq) > config['max_len_get_major_micro_species_major_tautomer']:
                        warnings.warn('Major tautomer calculations are limited to forms with length <= {}'.format(
                            config['max_len_get_major_micro_species_major_tautomer']), bpforms.BpFormsWarning)
                        structure = None
                    elif len(form.seq) > config['max_len_get_major_micro_species']:
                        warnings.warn('Major microspecies calculations are limited to forms with length <= {}'.format(
                            config['max_len_get_major_micro_species']), bpforms.BpFormsWarning)
                        structure = None
                    else:
                        structure = form.get_major_micro_species(ph, major_tautomer=major_tautomer, dearomatize=dearomatize)

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
            'seq': str(form),
            'length': len(form),
            'structure': smiles,
            'formula': formula,
            'mol_wt': mol_wt,
            'charge': charge,
            'warnings': warning_message,
        }


alphabet_ns = flask_restplus.Namespace('alphabet', description='List alphabets and get information about alphabets')
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
        Args:
            id (:obj:`str`): id of alphabet

        Returns:
            :obj:`dict`: dictionary representation of an alphabet
        """
        return get_alphabet(id)


@bpforms.core.cache.memoize(typed=False, expire=30 * 24 * 60 * 60)
def get_alphabet(id):
    """ Get an alphabet

    Args:
        id (:obj:`str`): id of alphabet

    Returns:
        :obj:`dict`: dictionary representation of an alphabet
    """
    try:
        alphabet_obj = bpforms.util.get_alphabet(id)
    except ValueError as error:
        flask_restplus.abort(400, 'Invalid alphabet "{}"'.format(id))

    alphabet_dict = alphabet_obj.to_dict()

    for code, monomer in alphabet_obj.monomers.items():
        alphabet_dict['monomers'][code] = get_monomer_properties(id, code)

    return alphabet_dict


@alphabet_ns.route("/<string:alphabet>/<string:monomer>/", defaults={'format': 'json'})
@alphabet_ns.route("/<string:alphabet>/<string:monomer>/<string:format>/")
@alphabet_ns.doc(params={
    'alphabet': 'Id of the alphabet of the biopolymer form (e.g. "dna", "rna", or "protein")',
    'monomer': 'Code of the monomeric form (e.g. "C" or "m2C")',
    'format': 'Output format ("emf", "eps", "jpeg", "json", "msbmp", "pdf", "png" or "svg")',
})
class MonomerResource(flask_restplus.Resource):
    """ Get information about a monomer """

    def get(self, alphabet, monomer, format):
        """ Get a monomeric form """
        """
        Returns:
            :obj:`object`: dictionary representation of an monomer or SVG-encoded image of a monomer
        """
        return get_monomer(alphabet, monomer, format)


@bpforms.core.cache.memoize(typed=False, expire=30 * 24 * 60 * 60)
def get_monomer(alphabet, monomer, format):
    """ Get a monomeric form

    Args:
        alphabet (:obj:`str`): id of the alphabet
        monomer (:obj:`str`): code of a monomeric form
        format (:obj:`str`): output format ("emf", "eps", "jpeg", "json", "msbmp", "pdf", "png" or "svg")

    Returns:
        :obj:`object`: dictionary representation of an monomer or SVG-encoded image of a monomer
    """
    try:
        alphabet_obj = bpforms.util.get_alphabet(alphabet)
    except ValueError as error:
        flask_restplus.abort(400, 'Invalid alphabet "{}"'.format(alphabet))

    monomer = urllib.parse.unquote(monomer)
    monomer_obj = alphabet_obj.monomers.get(monomer, None)
    if monomer_obj is None:
        flask_restplus.abort(400, 'Monomer "{}" not in alphabet "{}"'.format(monomer, alphabet))

    if format == 'json':
        return get_monomer_properties(alphabet, monomer)

    elif format in ['emf', 'eps', 'jpeg', 'msbmp', 'png', 'pdf', 'svg']:
        if format == 'emf':
            mimetype = 'image/emf'
        elif format == 'eps':
            mimetype = 'application/postscript'
        elif format == 'jpeg':
            mimetype = 'image/jpeg'
        elif format == 'msbmp':
            mimetype = 'image/bmp'
        elif format == 'png':
            mimetype = 'image/png'
        elif format == 'pdf':
            mimetype = 'application/pdf'
        elif format == 'svg':
            mimetype = 'image/svg+xml'

        return flask.Response(monomer_obj.get_image(image_format=format, width=250, height=150),
                              mimetype=mimetype)

    else:
        flask_restplus.abort(400, 'Invalid format "{}"'.format(format))


@bpforms.core.cache.memoize(typed=False, expire=30 * 24 * 60 * 60)
def get_monomer_properties(alphabet, monomer):
    """ Get properties of a monomeric form

    Args:
        alphabet (:obj:`str`): id of an alphabet
        monomer (:obj:`str`): code of monomeric form

    Returns:
        :obj:`dict`: properties of monomeric form
    """
    alphabet_obj = bpforms.util.get_alphabet(alphabet)
    form_obj = bpforms.util.get_form(alphabet)()
    monomer_obj = alphabet_obj.monomers.get(monomer, None)
    monomer_dict = monomer_obj.to_dict(alphabet=alphabet_obj)
    monomer_dict['bonds_backbone'] = len(monomer_obj.backbone_bond_atoms) > 0
    monomer_dict['bonds_left'] = form_obj.can_monomer_bond_left(monomer_obj)
    monomer_dict['bonds_right'] = form_obj.can_monomer_bond_right(monomer_obj)
    monomer_dict['formula'] = dict(monomer_obj.get_formula())
    monomer_dict['mol_wt'] = monomer_obj.get_mol_wt()
    monomer_dict['charge'] = monomer_obj.get_charge()
    return monomer_dict
