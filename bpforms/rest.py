""" REST JSON API

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

import bpforms
import bpforms.util
import flask
import flask_restplus
import flask_restplus.errors
import flask_restplus.fields

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

bpform_parser = bpform_ns.parser()
bpform_parser.add_argument('alphabet', type=str,
                           location='json',
                           help='Id of the alphabet of the biopolymer form')
bpform_parser.add_argument('monomer_seq', type=str,
                           location='json',
                           help='Sequence of monomers of the biopolymer form')
bpform_parser.add_argument('ph', type=float, default=None,
                           location='json',
                           help='pH at which to calculate the major microspecies of the biopolymer form')
bpform_parser.add_argument('major_tautomer', type=bool, default=True,
                           location='json',
                           help='If true, calculate the major tautomer')


@bpform_ns.route("/")
class Bpform(flask_restplus.Resource):
    """ Optionally, calculate the major protonation and tautomerization form a biopolymer form and calculate its properties """

    @bpform_ns.expect(bpform_parser)
    def post(self):
        """ Optionally, calculate the major protonation and tautomerization form a biopolymer form and calculate its properties """
        """
        Returns:
            :obj:`dict`
        """

        args = bpform_parser.parse_args(strict=True)
        alphabet = args['alphabet']
        monomer_seq = args['monomer_seq']
        ph = args['ph']
        major_tautomer = args['major_tautomer']

        try:
            form_cls = bpforms.util.get_form(alphabet)
        except ValueError as error:
            flask_restplus.abort(400, 'Invalid alphabet "{}"'.format(alphabet))

        try:
            form = form_cls().from_str(monomer_seq)
        except Exception as error:
            flask_restplus.abort(400, 'Unable to parse monomer sequence', details=str(error))

        if ph is not None:
            form.protonate(ph, major_tautomer=major_tautomer)
        return {
            'alphabet': alphabet,
            'monomer_seq': str(form),
            'length': len(form),
            'formula': dict(form.get_formula()),
            'mol_wt': form.get_mol_wt(),
            'charge': form.get_charge(),
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
