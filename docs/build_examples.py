""" Generate images for examples

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-14
:Copyright: 2019, Karr Lab
:License: MIT
"""

import bpforms
import os
import pkg_resources


default_dirname = pkg_resources.resource_filename('bpforms', os.path.join('web', 'img', 'example'))


def build(dirname=default_dirname):
    """ Build artifacts (e.g., images) for examples

    Args:
        dirname (:obj:`str`, optional): directory to save examples
    """
    # make directory for save images
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    # save images
    draw_polymer(bpforms.DnaForm, '{dI}ACGC', dirname, 'form-DNA',
                 show_atom_nums=False, width=203)
    draw_polymer(bpforms.RnaForm, 'AC{9A}GC', dirname, 'form-RNA',
                 show_atom_nums=False, width=203)
    draw_polymer(bpforms.ProteinForm, 'ACUC', dirname, 'form-Protein',
                 show_atom_nums=False, width=203)

    draw_monomer(bpforms.DnaForm, '''
        [id: "dI" 
            | name: "hypoxanthine"
            | structure: "O=C1NC=NC2=C1N=CN2"
            ]
        ''', dirname, 'structure-dI', show_atom_nums=True)

    draw_polymer(bpforms.DnaForm, '''
        [id: "dI" 
            | name: "hypoxanthine"
            | structure: "O=C1NC=NC2=C1N=CN2"
            | backbone-bond-atom: N10
            | backbone-displaced-atom: H10
            ]
        ''', dirname, 'monomer-backbone-bonds-dI', show_atom_nums=True)

    draw_monomer(bpforms.ProteinForm, '''
        [id: "AA0305" 
            | name: "N5-methyl-L-arginine"
            | structure: "O=C[C@H](CCCN(C(=[NH2])N)C)[NH3+]"
            | backbone-bond-atom: C2
            | backbone-displaced-atom: H2
            | right-bond-atom: C2
            | left-bond-atom: N15-1
            | left-displaced-atom: H15+1
            | left-displaced-atom: H15+1
            ]
        ''', dirname, 'left-right-bonds-AA0305', show_atom_nums=True)

    draw_polymer(bpforms.ProteinForm, '''
        A[id: "AA0305" 
            | name: "N5-methyl-L-arginine"
            | structure: "O=C[C@H](CCCN(C(=[NH2])N)C)[NH3+]"
            | backbone-bond-atom: C2
            | backbone-displaced-atom: H2
            | right-bond-atom: C2
            | left-bond-atom: N15-1
            | left-displaced-atom: H15+1
            | left-displaced-atom: H15+1            
            ]A
        ''', dirname, 'left-right-bonds-A-AA0305-A', show_atom_nums=False)

    draw_polymer(bpforms.DnaForm, '''
        AC | circular
        ''', dirname, 'circular-DNA-AC', show_atom_nums=False)

    draw_polymer(bpforms.ProteinForm, '''
        CAC | crosslink: [
            left-bond-atom: 1S1 |
            left-displaced-atom: 1H1 |
            right-bond-atom: 3S1 |
            right-displaced-atom: 3H1
        ]''', dirname, 'crosslink-protein-sulfide-bond', show_atom_nums=False)


def draw_monomer(Form, monomer, dirname, filename, show_atom_nums=False,
                 width=250, format='svg'):
    """ Generate and save an image of a monomer for an example

    Args:
        Form (:obj:`cls`): type of form (e.g., :obj:`bpforms.DnaForm`)
        monomer (:obj:`str`): string representation of the monomeric form
        dirname (:obj:`str`): directory to save image
        filename (:obj:`str`): filename to save image
        show_atom_nums (:obj:`bool`, optional): if :obj:`True`, show the numbers of the atoms
        width (:obj:`int`, optional): width of image
        format (:obj:`str`, optional): format for image
    """
    monomer = monomer.replace('\n', '').strip()
    form = Form().from_str(monomer)
    assert form.validate() == []
    if format == 'svg':
        mode = 'w'
    else:
        mode = 'wb'
    with open(os.path.join(dirname, filename + '.' + format), mode) as file:
        img = form.seq[0].get_image(width=width, height=150, image_format=format, show_atom_nums=show_atom_nums)
        file.write(img)


def draw_polymer(Form, polymer, dirname, filename, show_atom_nums=False,
                 width=250, format='svg'):
    """ Generate and save an image of a polymer for an example

    Args:
        Form (:obj:`cls`): type of form (e.g., :obj:`bpforms.DnaForm`)
        polymer (:obj:`str`): string representation of the biopolymer form
        dirname (:obj:`str`): directory to save image
        filename (:obj:`str`): filename to save image
        show_atom_nums (:obj:`bool`, optional): if :obj:`True`, show the numbers of the atoms
        width (:obj:`int`, optional): width of image
        format (:obj:`str`, optional): format for image
    """
    polymer = polymer.replace('\n', '').strip()
    form = Form().from_str(polymer)
    assert form.validate() == []
    if format == 'svg':
        mode = 'w'
    else:
        mode = 'wb'
    with open(os.path.join(dirname, filename + '.' + format), mode) as file:
        img = form.get_image(width=width, height=150, image_format=format, show_atom_nums=show_atom_nums)
        file.write(img)
