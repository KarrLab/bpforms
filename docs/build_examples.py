""" Generate images for examples

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-14
:Copyright: 2019, Karr Lab
:License: MIT
"""

import bpforms
import os
import pkg_resources
import wc_utils.util.chem


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
    draw_polymer(bpforms.ProteinForm, 'AC{U}C', dirname, 'form-Protein',
                 show_atom_nums=False, width=203)

    draw_polymer(bpforms.DnaForm, 'C{m2A}G{m2C}', dirname, 'alphabet-DNA',
                 show_atom_nums=False, width=203)
    draw_polymer(bpforms.RnaForm, 'A{21C}GC', dirname, 'alphabet-RNA',
                 show_atom_nums=False, width=203)
    draw_polymer(bpforms.ProteinForm, '{AA0037}E{AA0038}', dirname, 'alphabet-Protein',
                 show_atom_nums=False, width=203)

    draw_monomer(bpforms.DnaForm, '{m2A}', dirname, 'dna-m2A',
                 show_atom_nums=False, width=203)
    draw_monomer(bpforms.RnaForm, '{21C}', dirname, 'rna-21C',
                 show_atom_nums=False, width=203)
    draw_monomer(bpforms.ProteinForm, '{AA0037}', dirname, 'protein-AA0037',
                 show_atom_nums=False, width=203)

    draw_monomer(bpforms.DnaForm, '''
        [id: "m2C"
            | name: "2-O-methylcytosine"
            | structure: "COC1=NC(=CCN1C1CC(C(O1)COP(=O)([O-])[O-])O)N"
            ]
        ''', dirname, 'structure-m2C', show_atom_nums=True)

    draw_monomer(bpforms.ProteinForm, '''
        [id: "AA0305"
            | name: "N5-methyl-L-arginine"
            | structure: "OC(=O)[C@H](CCCN(C(=[NH2])N)C)[NH3+]"
            | r-bond-atom: C2
            | r-displaced-atom: O1
            | r-displaced-atom: H1
            | l-bond-atom: N16-1
            | l-displaced-atom: H16+1
            | l-displaced-atom: H16
            ]
        ''', dirname, 'left-right-bonds-AA0305', show_atom_nums=True)

    draw_polymer(bpforms.ProteinForm, '''
        A[id: "AA0305"
            | name: "N5-methyl-L-arginine"
            | structure: "OC(=O)[C@H](CCCN(C(=[NH2])N)C)[NH3+]"
            | r-bond-atom: C2
            | r-displaced-atom: O1
            | r-displaced-atom: H1
            | l-bond-atom: N16-1
            | l-displaced-atom: H16+1
            | l-displaced-atom: H16
            ]A
        ''', dirname, 'left-right-bonds-A-AA0305-A', show_atom_nums=False)

    draw_polymer(bpforms.DnaForm, '''
        AC | circular
        ''', dirname, 'circular-DNA-AC', show_atom_nums=False)

    draw_polymer(bpforms.ProteinForm, '''
        CAC | x-link: [
            l-bond-atom: 1S11 |
            l-displaced-atom: 1H11 |
            r-bond-atom: 3S11 |
            r-displaced-atom: 3H11
        ]''', dirname, 'crosslink-protein-sulfide-bond', show_atom_nums=False)

    draw_polymer(bpforms.ProteinForm, '''
        C:AC | x-link: [
            l-bond-atom: 1S11 |
            l-displaced-atom: 1H11 |
            r-bond-atom: 3S11 |
            r-displaced-atom: 3H11
        ]''', dirname, 'nick-protein', show_atom_nums=False)


def draw_monomer(Form, monomer, dirname, filename, show_atom_nums=False,
                 width=250, height=150, format='svg'):
    """ Generate and save an image of a monomer for an example

    Args:
        Form (:obj:`cls`): type of form (e.g., :obj:`bpforms.DnaForm`)
        monomer (:obj:`str`): string representation of the monomeric form
        dirname (:obj:`str`): directory to save image
        filename (:obj:`str`): filename to save image
        show_atom_nums (:obj:`bool`, optional): if :obj:`True`, show the numbers of the atoms
        width (:obj:`int`, optional): width of image
        height (:obj:`int`, optional): height of image
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
        img = form.seq[0].get_image(width=width, height=height, image_format=format, show_atom_nums=show_atom_nums)
        file.write(img)


def draw_polymer(Form, polymer, dirname, filename, show_atom_nums=False,
                 width=250, height=150, format='svg'):
    """ Generate and save an image of a polymer for an example

    Args:
        Form (:obj:`cls`): type of form (e.g., :obj:`bpforms.DnaForm`)
        polymer (:obj:`str`): string representation of the biopolymer form
        dirname (:obj:`str`): directory to save image
        filename (:obj:`str`): filename to save image
        show_atom_nums (:obj:`bool`, optional): if :obj:`True`, show the numbers of the atoms
        width (:obj:`int`, optional): width of image
        height (:obj:`int`, optional): height of image
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
        img = form.get_image(width=width, height=height, image_format=format, show_atom_nums=show_atom_nums)
        file.write(img)

    # print properties
    print(str(form))
    print('''<p class="form-properties">
        Length: {}<br/>
        Formula: {}<br/>
        Molecular weight: {:.1f}<br/>
        Charge: {}
        </p>'''.format(
        len(form.seq),
        ''.join('{}<sub>{}</sub>'.format(el, int(count)) for el, count in form.get_formula().items()),
        form.get_mol_wt(),
        form.get_charge()))
