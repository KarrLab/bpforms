""" Code to help build alphabets

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-08-14
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bpforms.core import (cache, Monomer, Identifier, IdentifierSet, SynonymSet)
from xml.etree.ElementTree import ElementTree
import io
import mendeleev
import openbabel
import os
import pkg_resources
import requests
import tarfile


def download_pdb_ccd():
    """ Download PDB CCD

    Returns:
        :obj:`str`: path to tar.gz file for the PDB CCD
    """
    dirname = pkg_resources.resource_filename('bpforms', os.path.join('alphabet', 'pdb'))
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    filename = os.path.join(dirname, 'components-pub-xml.tar.gz')
    if not os.path.isfile(filename):
        response = requests.get('http://ligand-expo.rcsb.org/dictionaries/components-pub-xml.tar.gz')
        response.raise_for_status()
        with open(filename, 'wb') as file:
            file.write(response.content)

    return filename


@cache.memoize(typed=False, expire=30 * 24 * 60 * 60)  # , filename_args=[0])
def parse_pdb_ccd(filename, valid_types, max_monomers):
    """ Parse entries out of the PDB CCD

    Args:
        filename (:obj:`str`): path to tar.gz file for PDB CCD
        valid_types (:obj:`tuple` of :obj:`str`): list
            of types of entries to retrieve
        max_monomers (:obj:`float`): maximum number of
            entries to process

    Returns:
        :obj:`list` of :obj:`tuple`:  list of metadata and
            structures of the entries
    """
    entries = []
    with tarfile.open(filename, 'r:gz') as tar_file:
        i_file = 0
        n_files = len(tar_file.getmembers())
        for file_info in tar_file:
            i_file += 1
            if i_file % 1000 == 1:
                print('Processing file {} of {}'.format(i_file, n_files))

            if os.path.splitext(file_info.name)[-1] != '.xml':
                continue

            xml_file = tar_file.extractfile(file_info)
            entry = parse_pdb_ccd_entry(
                xml_file, valid_types)
            if entry is not None:
                entries.append(entry)
            if max_monomers is not None and len(entries) >= max_monomers:
                break
    return entries


def parse_pdb_ccd_entry(xml_file, valid_types):
    """ Parse an entry of the PDB CCD

    Args:
        xml_file (:obj:`io.BufferedReader`): XML file
            that defines an entry of the PDB CCD
        valid_types (:obj:`list` of :obj:`str`): list
            of types of entries to retrieve

    Returns:
        :obj:`tuple`:

            * :obj:`Monomer`: metadata about the entry
            * :obj:`str`: id of base monomer
            * :obj:`str`: SMILES-encoded structure of the entry
            * :obj:`dict`: structure of the entry
            * :obj:`dict`: dictionary that maps atom ids to their
                coordinates
    """
    ns = '{http://pdbml.pdb.org/schema/pdbx-v40.xsd}'

    xml_root = ElementTree().parse(xml_file)
    xml_group = xml_root.find(ns + 'chem_compCategory')
    if xml_group is None:
        return  # pragma: no cover # element is always present

    # get id
    xml_comp = xml_group.find(ns + 'chem_comp')
    if xml_comp is None:
        return  # pragma: no cover # element is always present
    id = xml_comp.get('id')
    identifiers = IdentifierSet([Identifier('pdb-ccd', id)])

    # check that compound has been released, is an amino acid, and is no ambiguous
    xml_el = xml_comp.find(ns + 'pdbx_release_status')
    if xml_el is None or xml_el.text != 'REL':
        return

    xml_el = xml_comp.find(ns + 'type')
    if xml_el is None or xml_el.text not in valid_types:
        return

    xml_el = xml_comp.find(ns + 'pdbx_ambiguous_flag')
    if xml_el is None or xml_el.text != 'N':
        return  # pragma: no cover # element is always present

    # get name
    xml_el = xml_comp.find(ns + 'name')
    if xml_el is None:
        name = None  # pragma: no cover # element is always present
    else:
        name = xml_el.text.lower()

    # retrieve synonyms
    synonyms = SynonymSet()
    xml_el = xml_comp.find(ns + 'one_letter_code')
    if xml_el is not None:
        synonyms.add(xml_el.text)

    for xml_group in xml_root.findall(ns + 'pdbx_chem_comp_identifierCategory'):
        for xml_subgroup in xml_group.findall(ns + 'pdbx_chem_comp_identifier'):
            if xml_subgroup.get('type') == 'SYSTEMATIC NAME':
                synonyms.add(xml_subgroup.find(ns + 'identifier').text)

    xml_el = xml_comp.find(ns + 'mon_nstd_parent_comp_id')
    if xml_el is not None:
        base_monomer = xml_el.text
    else:
        base_monomer = None

    # retrieve structure
    smiles = None
    for xml_group in xml_root.findall(ns + 'pdbx_chem_comp_descriptorCategory'):
        for xml_subgroup in xml_group.findall(ns + 'pdbx_chem_comp_descriptor'):
            if xml_subgroup.get('type') == 'SMILES_CANONICAL' and \
                    xml_subgroup.get('program') == 'OpenEye OEToolkits':
                smiles = xml_subgroup.find(ns + 'descriptor').text
    if smiles is None:
        return  # pragma: no cover # element is always present

    # discard entries with coordinating metals
    mol = openbabel.OBMol()
    inchi_conv = openbabel.OBConversion()
    assert inchi_conv.SetInFormat('smi'), 'Unable to set format to SMILES'
    assert inchi_conv.SetOutFormat('inchi'), 'Unable to set format to InChI'
    inchi_conv.ReadString(mol, smiles)
    inchi = inchi_conv.WriteString(mol)
    formula = inchi.split('/')[1]
    if '.' in formula:
        return

    mol = {
        'complete': True,
        'atoms': [],
        'bonds': [],
    }
    atoms = {}

    xml_group = xml_root.find(ns + 'chem_comp_atomCategory')
    xml_atoms = xml_group.findall(ns + 'chem_comp_atom')
    for xml_atom in xml_atoms:
        if xml_atom.get('comp_id') != id:
            continue

        atom_id = xml_atom.get('atom_id')
        charge = int(float(xml_atom.find(ns + 'charge').text))
        element = xml_atom.find(ns + 'type_symbol').text
        element = element[0] + element[1:].lower()

        mol['atoms'].append({
            'atomic_num': getattr(mendeleev, element).atomic_number,
            'charge': charge,
        })

        atoms[atom_id] = {
            'position': int(float(xml_atom.find(ns + 'pdbx_ordinal').text)),
            'element': element,
            'charge': charge,
        }

    xml_group = xml_root.find(ns + 'chem_comp_bondCategory')
    xml_bonds = xml_group.findall(ns + 'chem_comp_bond')
    for xml_bond in xml_bonds:
        if xml_bond.get('comp_id') != id:
            continue

        i_atom_1 = atoms[xml_bond.get('atom_id_1')]['position']
        i_atom_2 = atoms[xml_bond.get('atom_id_2')]['position']

        order_str = xml_bond.find(ns + 'value_order').text
        if order_str == 'sing':
            order = 1
        elif order_str == 'doub':
            order = 2
        elif order_str == 'trip':
            order = 3
        elif order_str == 'quad':
            order = 4
        elif order_str == 'arom':
            order = 5
        else:  # ['delo', 'pi', 'poly']
            mol['complete'] = False

        mol['bonds'].append({
            'begin': i_atom_1,
            'end': i_atom_2,
            'order': order,
        })

    return (Monomer(id=id, name=name, synonyms=synonyms, identifiers=identifiers),
            base_monomer, smiles, mol, atoms)


def get_pdb_ccd_open_babel_mol(pdb_mol):
    """ Generate an Open Babel representation of a PDB CCD entry

    Args:
        pdb_mol (:obj:`dict`): structure of a entry

    Returns:
        :obj:`openbabel.OBMol`: structure of a entry
    """
    if not pdb_mol['complete']:
        return None

    mol = openbabel.OBMol()

    for pdb_atom in pdb_mol['atoms']:
        atom = openbabel.OBAtom()
        atom.SetAtomicNum(pdb_atom['atomic_num'])
        atom.SetFormalCharge(pdb_atom['charge'])
        mol.AddAtom(atom)

    for pdb_bond in pdb_mol['bonds']:
        bond = openbabel.OBBond()
        bond.SetBegin(mol.GetAtom(pdb_bond['begin']))
        bond.SetEnd(mol.GetAtom(pdb_bond['end']))
        bond.SetBondOrder(pdb_bond['order'])
        assert mol.AddBond(bond)

    return mol


def get_can_smiles(mol):
    """ Get the canonical SMILES representation of a molecule without its stereochemistry

    Args:
        mol (:obj:`openbabel.OBMol`): molecule

    Returns:
        :obj:`str`: SMILES representation of a molecule without its stereochemistry
    """
    conv = openbabel.OBConversion()
    assert conv.SetInFormat('smi')
    assert conv.SetOutFormat('smi')
    conv.SetOptions('c', conv.OUTOPTIONS)
    smiles = conv.WriteString(mol, True)
    mol = openbabel.OBMol()
    conv.ReadString(mol, smiles)

    for atom in openbabel.OBMolAtomIter(mol):
        atom.UnsetStereo()

    for bond in openbabel.OBMolBondIter(mol):
        bond.UnsetHash()
        bond.UnsetWedge()
        bond.UnsetUp()
        bond.UnsetDown()

    mol.DeleteData(openbabel.StereoData)

    conv = openbabel.OBConversion()
    assert conv.SetInFormat('smi')
    assert conv.SetOutFormat('smi')
    conv.SetOptions('c', conv.OUTOPTIONS)
    smiles = conv.WriteString(mol, True)

    conv = openbabel.OBConversion()
    assert conv.SetInFormat('smi')
    assert conv.SetOutFormat('smi')
    conv.SetOptions('c', conv.OUTOPTIONS)
    mol = openbabel.OBMol()
    conv.ReadString(mol, smiles)
    smiles = conv.WriteString(mol, True)

    return smiles
