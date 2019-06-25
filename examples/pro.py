""" Generate BpForms for all of the proteins in PRO, verify
them, and calculate their properties

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-24
:Copyright: 2019, Karr Lab
:License: MIT
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import bpforms
import csv
import os
import pickle
import re
import requests
import requests_cache

IN_URL = 'https://proconsortium.org/download/current/pro_nonreasoned.obo'
OBO_FILENAME = os.path.join('examples', 'pro_nonreasoned.obo')
PKL_FILENAME = os.path.join('examples', 'pro_nonreasoned.pkl')
PKL_PARSED_FILENAME = os.path.join('examples', 'pro_nonreasoned.parsed.pkl')
UNIPROT_ENDPOINT = 'https://www.uniprot.org/uniprot/{}.fasta'

OUT_PICKLE_FILENAME = os.path.join('examples', 'pro_nonreasoned.pkl')
OUT_TSV_FILENAME = os.path.join('examples', 'pro_nonreasoned.tsv')
OUT_FASTA_FILENAME = os.path.join('examples', 'pro_nonreasoned.fasta')
IN_MONOMERS_FILENAME = os.path.join('examples', 'pro_nonreasoned.monomers.csv')

cache_name = os.path.join('examples', 'pro')
session = requests_cache.core.CachedSession(cache_name, backend='sqlite', expire_after=None)
session.mount('https://www.uniprot.org/', requests.adapters.HTTPAdapter(max_retries=5))

AA_CHARS_TO_CODES = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Cys': 'C',
    'Glu': 'E',
    'Gln': 'Q',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V',
}


def run(in_monomers_filename=IN_MONOMERS_FILENAME,
        out_pickle_filename=OUT_PICKLE_FILENAME, out_tsv_filename=OUT_TSV_FILENAME, out_fasta_filename=OUT_FASTA_FILENAME):
    """ Download PRO ontology, generate proteoforms, and encode with BpForms

    Args:
        in_monomers_filename (:obj:`str`, optional): path to list of ids of monomeric forms used 
            by PRO and their alphabet code in tab-separated format
        out_pickle_filename (:obj:`str`, optional): path to save results in pickle format
        out_tsv_filename (:obj:`str`, optional): path to save results in tab-separated format
        out_fasta_filename (:obj:`str`, optional): path to save results in FASTA format        

    Returns:
        :obj:`list` of :obj:`dict`: proteoforms encoded with BpForms
    """
    # get the PRO ontology and extract the modified proteins from the ontology
    proteins = get_pro()

    # parse the modified proteins and retrieve their sequences
    parsed_proteins = []
    for i_protein, protein in enumerate(proteins):
        parsed_proteins.append(parse_protein(protein))

    # save the parsed proteins in pickle format
    with open(out_pickle_filename, 'wb') as file:
        pickle.dump(parsed_proteins, file)

    # read list of monomers
    monomers = {}
    with open(in_monomers_filename, 'r') as file:
        reader = csv.DictReader(file, dialect='excel')
        for row in reader:
            monomers[row['PRO id']] = bpforms.protein_alphabet.monomers.get(row['BpForms code'], None)

    # generate list of modified monomeric forms
    for protein in parsed_proteins:
        for modification in protein['modifications']:
            if modification['monomer'] not in monomers:
                monomers[modification['monomer']] = None

    # print list of unmapped monomers
    unmapped_monomers = []
    for monomer, code in monomers.items():
        if not code:
            unmapped_monomers.append(monomer)
    unmapped_monomers.sort()
    if unmapped_monomers:
        print('Several PRO monomeric forms have not been mapped to BpForms monomeric forms:\n  {}'.format(
            '\n  '.join(unmapped_monomers)))

    # generate BpForms for each protein
    for i_protein, protein in enumerate(parsed_proteins):
        if i_protein % 10 == 0:
            print('Generating BpForms {} of {}'.format(i_protein + 1, len(parsed_proteins)))

        protein['modified_seq'] = None
        if not protein['uniprot_id']:
            continue
        if not protein['seq']:
            continue
        if not protein['modifications']:
            continue
        if protein['errors']:
            continue
        form = gen_bpform(protein, monomers)
        protein['modified_seq'] = str(form)
        if not form.validate():
            protein['formula'] = form.get_formula()
            protein['mol_wt'] = form.get_mol_wt()
            protein['charge'] = form.get_charge()

    # save the proteoforms in TSV format
    with open(out_tsv_filename, 'w') as file:
        writer = csv.writer(file, dialect='excel-tab')
        writer.writerow(['PRO id', 'UniProt id',
                         'Unmodified sequence (IUBMB)', 'Processing', 'Modifications',
                         'Proteoform sequence (BpForms)', 'Concrete', 'Formula', 'Molecular weight', 'Charge',
                         'Issues'])

        for parsed_protein in parsed_proteins:
            if parsed_protein['errors']:
                errors = '. '.join(parsed_protein['errors']) + '.'
            else:
                errors = ''
            writer.writerow([
                parsed_protein['id'],
                parsed_protein['uniprot_id'],
                parsed_protein['seq'],
                ', '.join('{}-{}'.format(p['start'], p['end']) for p in parsed_protein['processing']),
                ', '.join('{} --> {} ({})'.format(m['residue'] or '?', m['monomer'], ', '.join(str(p) for p in m['positions']))
                          for m in parsed_protein['modifications']),
                parsed_protein.get('modified_seq', None),
                parsed_protein.get('concrete', False),
                parsed_protein.get('formula', None),
                parsed_protein.get('mol_wt', None),
                parsed_protein.get('charge', None),
                errors,
            ])

    # save the proteoforms in FASTA format
    seqs = (SeqRecord(id='{} | {}'.format(protein['id'], protein['uniprot_id']),
                      seq=Seq(protein['modified_seq']),
                      description='')
            for protein in parsed_proteins
            if protein['modified_seq'])
    SeqIO.write(seqs, out_fasta_filename, "fasta")

    # return proteins
    return parsed_proteins


def get_pro(obo_filename=OBO_FILENAME, pkl_filename=PKL_FILENAME):
    """ Get the PRO ontology and extract the modified proteins from the ontology

    Args:
        obo_filename (:obj:`str`, optional): filename to save PRO in OBO format
        pkl_filename (:obj:`str`, optional): filename to save/read PRO from pickled file

    Returns:
        :obj:`list` of :obj:`dict`: list of PRO ontology terms for modified proteins
    """
    # download or read PRO
    if not os.path.isfile(obo_filename):
        # download PRO
        response = requests.get(IN_URL)
        response.raise_for_status()
        with open(obo_filename, 'wb') as file:
            file.write(response.content)

        # parse PRO
        proteins = []
        protein = None
        with open(obo_filename, 'r') as file:
            for line in file:
                line = line.rstrip('\n')
                if line.startswith('['):
                    if line.startswith('[Term]'):
                        protein = {}
                    else:
                        protein = None
                elif line and protein is not None:
                    key, _, value = line.partition(': ')
                    if key not in protein:
                        protein[key] = []
                    protein[key].append(value)

                    if key == 'comment' and value.startswith('Category=organism-modification.'):
                        proteins.append(protein)

        # save PRO in pickle format
        with open(pkl_filename, 'wb') as file:
            pickle.dump(proteins, file)
    else:
        # load PRO from pickle format
        with open(pkl_filename, 'rb') as file:
            proteins = pickle.load(file)

    # return PRO
    return proteins


def parse_protein(protein):
    """ Parse the modification information from a term for a modified protein

    Args:
        protein (:obj:`dict`): term for a modified protein

    Returns:
        :obj:`dict` with PRO id, UniProt id, processing start position, processing end position, unmodified sequence, and modifications
    """
    assert len(protein['id']) == 1
    id = protein['id'][0]
    errors = []

    seq_synonyms = []
    for synonym in protein.get('synonym', []):
        if synonym.startswith('"UniProtKB:') and '" EXACT PRO-proteoform-std ' in synonym:
            seq_synonyms.append(synonym)
    if not seq_synonyms:
        errors.append('No synonym which defines a modified sequence')
        return {
            'id': id,
            'uniprot_id': None,
            'processing': [],
            'modifications': [],
            'seq': None,
            'errors': errors,
        }
    elif len(seq_synonyms) > 1:
        errors.append('Multiple synonyms which define modified sequences')
    synonym = seq_synonyms[0]

    uniprot_id, _, processing_modifications_type = synonym.partition(', ')
    uniprot_id = uniprot_id.partition(':')[2]

    response = session.get(UNIPROT_ENDPOINT.format(uniprot_id))
    response.raise_for_status()
    seq = response.content.decode('utf-8').partition('\n')[2].replace('\n', '')
    if not seq:
        errors.append('No sequence for UniProt entry; entry may be deprecated')

    processing_modifications = processing_modifications_type.partition('"')[0]
    processing = []
    while True:
        match = re.match(r'^(\?|\d+)\-(\?|\d+)(, |$)', processing_modifications)
        if match:
            if match.group(1) == '?':
                start = None
                errors.append('Unknown processing start position')
            else:
                start = int(float(match.group(1)))
                if start <= 0 or start > len(seq):
                    errors.append('Start position must be within sequence')
            if match.group(2) == '?':
                end = None
                errors.append('Unknown processing end position')
            else:
                end = int(float(match.group(2)))
                if end <= 0 or end > len(seq):
                    errors.append('End position must be within sequence')

            if start and end and start > end:
                errors.append('End position must be after start position')

            processing.append({
                'start': start,
                'end': end,
            })
            processing_modifications = processing_modifications[len(match.group(0)):]
        else:
            break
    if processing_modifications:
        modifications_str = processing_modifications.split('|')
    else:
        modifications_str = []

    modifications = []
    for modification in modifications_str:
        if modification and modification[0] == '(' and modification[-1] == ')':
            modification = modification[1:-1]

        if ' or ' in modification or ' and/or ' in modification:
            errors.append('Boolean logic not supported')
            continue

        if ', ' in modification:
            residue_positions, _, monomer = modification.partition(', ')
            if monomer == 'PR:000026291':  # unmodified monomer
                continue

            residue_codes = set()
            positions = []
            for residue_position in residue_positions.split('/'):
                residue_chars, _, position = residue_position.partition('-')
                residue_code = AA_CHARS_TO_CODES[residue_chars]
                position = int(float(position))
                if seq[position - 1] != residue_code:
                    errors.append('Position {} != {}'.format(position, residue_code))
                residue_codes.add(residue_code)
                positions.append(position)
            if len(residue_codes) != 1:
                residue_code = None
                errors.append('Residues {{{}}} annotated with the same modification {}'.format(
                    ', '.join(residue_codes), monomer))
        else:
            residue_code = None
            positions = []
            monomer = modification

        modifications.append({
            'residue': residue_code,
            'positions': positions,
            'monomer': monomer
        })

    return {
        'id': id,
        'uniprot_id': uniprot_id,
        'processing': processing,
        'modifications': modifications,
        'seq': seq,
        'errors': errors,
    }


def gen_bpform(protein, pro_ids_to_bpform_monomers):
    """ Generate BpForm for a modified protein in PRO

    Args:
        protein (:obj:`dict`): term for modified protein
        pro_ids_to_bpform_monomers (:obj:`dict`): dictionary which maps ids of monomeric forms
            used by PRO to monomeric forms in the BpForms protein alphabet

    Returns:
        :obj:`bpforms.ProteinForm`: BpForm for a term in PRO
    """
    form = bpforms.ProteinForm()
    monomers = bpforms.protein_alphabet.monomers

    # generate BpForm for unmodified sequence
    for base in protein['seq']:
        form.seq.append(monomers[base])

    # apply modifications
    concrete = True
    for modification in protein['modifications']:
        monomer = pro_ids_to_bpform_monomers[modification['monomer']]
        if monomer is None:
            concrete = False
            monomer = bpforms.Monomer(
                id=modification['monomer'])
            if modification['residue']:
                monomer.base_monomers.add(monomers[modification['residue']])

        if modification['positions']:
            for position in modification['positions']:
                if form.seq[position - 1] == monomers[protein['seq'][position - 1]]:
                    form.seq[position - 1] = monomer
                else:
                    protein['errors'].append('Unable to set monomeric form at position {}'.format(
                        position))

        elif modification['residue']:
            concrete = False
            monomer.start_position = protein['seq'].find(modification['residue']) + 1
            monomer.end_position = protein['seq'].rfind(modification['residue']) + 1
            form.seq[monomer.start_position - 1] = monomer
            set_monomer = False
            for i_monomer in range(monomer.start_position, monomer.end_position + 1):
                if form.seq[i_monomer - 1] == monomers[protein['seq'][i_monomer - 1]]:
                    set_monomer = True
                    form.seq[i_monomer] = monomer
                    break
            if not set_monomer:
                protein['errors'].append('Unable to set monomeric form')
        else:
            concrete = False
            monomer.start_position = 1
            monomer.end_position = len(protein['seq'])
            set_monomer = False
            for i_monomer in range(monomer.start_position, monomer.end_position + 1):
                if form.seq[i_monomer - 1] == monomers[protein['seq'][i_monomer - 1]]:
                    form.seq[i_monomer] = monomer
                    set_monomer = True
                    break
            if not set_monomer:
                protein['errors'].append('Unable to set monomeric form')

    # apply processing
    if protein['processing']:
        procesed_seq = []
        for processing in protein['processing']:
            procesed_seq.extend(form.seq[processing['start']-1:processing['end']])
        form.seq = procesed_seq

    # validate
    protein['concrete'] = concrete
    protein['errors'].extend(form.validate())    

    # return proteoform represented with BpForms
    return form
