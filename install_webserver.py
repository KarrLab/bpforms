""" Install website

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-17
:Copyright: 2019, Karr Lab
:License: MIT
"""

import bcforms.core
import bcforms.rest
import bpforms.core
import bpforms.rest
import bpforms.util
import importlib
import os.path
import pathlib
import pkg_resources
import shutil
import subprocess
import sys

sys.path.append('docs')
sys.path.append('examples')
import build_examples
import modomics
import pro

def build(alphabet_ids=None, pro_max_num_proteins=None):
    """ Install website

    * Clear REST cache
    * Cache REST queries
    * Build images of monomers for HTML pages
    * Build images for examples
    * Restart server

    Args:
        alphabet_ids (:obj:`list` of :obj:`str`): list of ids of alphabets to cache; if :obj:`None`,
            cache all alphabets
        max_num_proteins (:obj:`int`, optional): maximum number of proteins in PRO to analyze
    """
    rest_client = bpforms.rest.app.test_client()

    if alphabet_ids is None:
        alphabets = bpforms.util.get_alphabets().values()
    else:
        alphabets = [bpforms.util.get_alphabet(alphabet_id) for alphabet_id in alphabet_ids]

    # pull from GitHub
    p = subprocess.Popen(['git', 'pull'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    err = p.communicate()[1].decode()
    if p.returncode != 0:
        raise Exception(err)
    if err:
        print(err, file=sys.stderr)

    # clear cache
    bpforms.core.cache.clear()
    bcforms.core.cache.clear()

    # cache alphabet REST queries and save JSON files for HTML pages
    data_dir = pkg_resources.resource_filename('bpforms', os.path.join('web', 'data'))
    alphabet_data_dir = pkg_resources.resource_filename('bpforms', os.path.join('web', 'data', 'alphabet'))
    if not os.path.isdir(data_dir):
        os.makedirs(data_dir)
    if not os.path.isdir(alphabet_data_dir):
        os.makedirs(alphabet_data_dir)

    rv = rest_client.get('/api/alphabet/')
    assert rv.status_code == 200
    with open(os.path.join(data_dir, 'alphabets.json'), 'wb') as file:
        file.write(rv.data)

    for alphabet in alphabets:
        rv = rest_client.get('/api/alphabet/' + alphabet.id + '/')
        assert rv.status_code == 200
        with open(os.path.join(alphabet_data_dir, alphabet.id + '.json'), 'wb') as file:
            file.write(rv.data)

    # cache crosslink REST queries and save json files for HTML pages
    bcforms_rest_client = bcforms.rest.app.test_client()
    rv = bcforms_rest_client.get('/api/crosslink/')
    assert rv.status_code == 200
    with open(os.path.join(data_dir, 'xlink' + '.json'), 'wb') as file:
        file.write(rv.data)

    # build images of monomers for alphabet web pages
    img_dir = pkg_resources.resource_filename('bpforms', os.path.join('web', 'img', 'alphabet'))

    for alphabet in alphabets:
        alphabet_img_dir = os.path.join(img_dir, alphabet.id)
        if not os.path.isdir(alphabet_img_dir):
            os.makedirs(alphabet_img_dir)

        for code, monomer in alphabet.monomers.items():
            with open(os.path.join(alphabet_img_dir, code + '.png'), 'wb') as file:
                file.write(monomer.get_image(image_format='png', width=250, height=150))

    # build images of crosslinks for webpage
    xlink_img_dir = pkg_resources.resource_filename('bpforms', os.path.join('web', 'img', 'crosslink'))
    if not os.path.isdir(xlink_img_dir):
        os.makedirs(xlink_img_dir)

    xlinks = list(bcforms.core.parse_yaml(bcforms.core._xlink_filename).keys())
    for xlink in xlinks:
        with open(os.path.join(xlink_img_dir, xlink + '.png'), 'wb') as file:
            file.write(bcforms.core.draw_xlink(xlink))

    # build examples
    build_examples.build()
    modomics.run()
    if pro_max_num_proteins is not None and pro_max_num_proteins:
        pro.run(max_num_proteins=pro_max_num_proteins)

    # restart server
    restart_filename = os.path.join('..', 'tmp', 'restart')
    if os.path.isfile(restart_filename):
        pathlib.Path(restart_filename).touch()


if __name__ == "__main__":
    build()
