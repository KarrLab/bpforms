""" REST JSON API

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-02-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

import bpforms
import bpforms.util
import flask_api
import flask_api.exceptions
import flask_api.response


class FlaskApi(flask_api.FlaskAPI):
    def __init__(self, name):
        """
        Args:
            name (:obj:`str`): name
        """
        super(FlaskApi, self).__init__(name)

    def handle_api_exception(self, exception):
        """ Handle API exception

        Args:
            exception (:obj:`flask_api.exceptions.APIException`): exception

        Returns:
            :obj:`flask_api.response.APIResponse`: response
        """
        if hasattr(exception, 'to_dict'):
            content = exception.to_dict()
        else:
            content = {'message': exception.detail}
        return flask_api.response.APIResponse(content, status=exception.status_code)


class ApiException(flask_api.exceptions.APIException):
    """ API exception to raise HTTP response with 400 error code

    Attributes:
        message (:obj:`str`): summary of the exception
        details (:obj:`dict`): details of the exception
    """
    status_code = 400
    detail = ''

    def __init__(self, message, **details):
        """
        Args:
            message (:obj:`str`): summary of the exception
            details (:obj:`dict`, optional): details of the exception
        """
        self.message = message
        self.details = details

    def to_dict(self):
        """ Generate a dictionary representation of the exception

        Returns:
            :obj:`dict`: dictionary representation of the exception
        """
        rv = dict(self.details)
        rv['message'] = self.message
        return rv

    def __str__(self):
        """ Generate a string representation of the exception

        Returns:
            :obj:`str`: summary of the exception
        """
        return self.message


app = FlaskApi(__name__)


@app.route("/api/")
def default():
    return {
        'description': ('The BpForms JSON REST API is a REST API for validating biopolymer forms '
                        'and calculating their lengths, formulae, molecular weights, and charges.'),
        'endpoints': [
            {
                '/api/alphabet': 'Get a list of available alphabets',
                '/api/alphabet/{alphabet: string}': 'Get the bases in an alphabet',
                '/api/bpform/properties/{alphabet: string}/{base sequence: string}(/{pH: float})?':
                ('Optionally protonates the biopolymer and then calculates the length, formula, molecular weight, '
                 'and charge of the biopolymer form specified by the alphabet (dna, rna, or protein) and base '
                 'sequence, and returns these properties as an associative array.',),
            },
        ],
        'version': bpforms.__version__,
    }


@app.route("/api/bpform/properties/<string:alphabet>/<string:base_seq>/", defaults={'ph': None})
@app.route("/api/bpform/properties/<string:alphabet>/<string:base_seq>/<string:ph>/")
def get_bpform_properties(alphabet, base_seq, ph):
    try:
        form_cls = bpforms.util.get_form(alphabet)
    except ValueError as error:
        raise ApiException('Invalid alphabet "{}"'.format(alphabet))

    try:
        form = form_cls().from_str(base_seq)
    except Exception as error:
        raise ApiException('Unable to parse base sequence', details=str(error))

    if ph is not None:
        try:
            ph = float(ph)
        except Exception:
            raise ApiException('pH must be a float')
        form.protonate(ph)
    return {
        'alphabet': alphabet,
        'base_seq': str(form),
        'length': len(form),
        'formula': dict(form.get_formula()),
        'mol_wt': form.get_mol_wt(),
        'charge': form.get_charge(),
    }


@app.route("/api/alphabet/")
def get_alphabets():
    return ['dna', 'rna', 'protein']


@app.route("/api/alphabet/<string:alphabet>/")
def get_alphabet(alphabet):
    try:
        alphabet_obj = bpforms.util.get_alphabet(alphabet)
    except ValueError as error:
        raise ApiException('Invalid alphabet "{}"'.format(alphabet))

    return alphabet_obj.to_dict()
