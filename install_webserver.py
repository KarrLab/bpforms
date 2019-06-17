""" Install website

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-17
:Copyright: 2019, Karr Lab
:License: MIT
"""

import bpforms.rest
import bpforms.util
import os.path
import pathlib
import pkg_resources
import shutil
import subprocess
import sys

sys.path.append('docs')
import build_examples


def build(alphabet_ids=None):
    """ Install website

    * Clear REST cache
    * Cache REST queries
    * Build images of monomers for HTML pages
    * Build images for examples
    * Restart server

    Args:
        alphabet_ids (:obj:`list` of :obj:`str`): list of ids of alphabets to cache; if :obj:`None`, 
            cache all alphabets
    """

    # pull from GitHub
    p = subprocess.Popen(['git', 'pull'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    err = p.communicate()[1].decode()
    if p.returncode != 0:
        raise Exception(err)
    if err:
        print(err, file=sys.stderr)

    # clear cache
    cache_dirname = os.path.join(os.path.expanduser('~'), '.cache', 'bpforms')
    if os.path.isdir(cache_dirname):
        shutil.rmtree(cache_dirname)

    # build images of monomers for alphabet web pages
    # cache alphabet REST queries
    img_dir = pkg_resources.resource_filename('bpforms', os.path.join('web', 'img', 'alphabet'))

    if alphabet_ids is None:
        alphabets = bpforms.util.get_alphabets().values()
    else:
        alphabets = [bpforms.util.get_alphabet(alphabet_id) for alphabet_id in alphabet_ids]

    for alphabet in alphabets:
        alphabet_img_dir = os.path.join(img_dir, alphabet.id)
        if not os.path.isdir(alphabet_img_dir):
            os.makedirs(alphabet_img_dir)

        for code, monomer in alphabet.monomers.items():
            with open(os.path.join(alphabet_img_dir, code + '.png'), 'wb') as file:
                file.write(monomer.get_image(image_format='png', width=250, height=150))

        bpforms.rest.get_alphabet(alphabet.id)

    # build examples
    build_examples.build()

    # restart server
    restart_filename = os.path.join('..', 'tmp', 'restart')
    if os.path.isfile(restart_filename):
        pathlib.Path(restart_filename).touch()


if __name__ == "__main__":
    build()
