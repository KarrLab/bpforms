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
from matplotlib import pyplot
from xml.etree import ElementTree
import bpforms
import copy
import csv
import matplotlib
import numpy
import os
import pickle
import re
import requests
import requests_cache

IN_URL = 'https://proconsortium.org/download/current/pro_nonreasoned.obo'
IN_OBO_FILENAME = os.path.join('examples', 'pro_nonreasoned.obo')
IN_PKL_FILENAME = os.path.join('examples', 'pro_nonreasoned.pkl')
IN_TSV_FILELANE = os.path.join('examples', 'pro_input.in.tsv') # from Darren Natale
IN_MONOMERS_FILENAME = os.path.join('examples', 'pro.monomers.csv')
UNIPROT_SEQ_ENDPOINT = 'https://www.uniprot.org/uniprot/{}.fasta'
UNIPROT_XML_ENDPOINT = 'https://www.uniprot.org/uniprot/{}.xml'

OUT_PICKLE_FILENAME = os.path.join('examples', 'pro_input.out.pkl')
OUT_PICKLE_FILENAME_2 = os.path.join('examples', 'pro_input.out.2.pkl')
OUT_TSV_FILENAME = os.path.join('examples', 'pro_input.out.tsv')
OUT_FASTA_FILENAME = os.path.join('examples', 'pro_input.fasta')
OUT_FIG_FILENAME = os.path.join('examples', 'pro_input.svg')

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


def run(in_obo_filename=IN_OBO_FILENAME, in_pkl_filename=IN_PKL_FILENAME, in_tsv_filename=IN_TSV_FILELANE,
        in_monomers_filename=IN_MONOMERS_FILENAME,
        max_num_proteins=None,
        out_pickle_filename=OUT_PICKLE_FILENAME, out_pickle_filename_2=OUT_PICKLE_FILENAME_2,
        out_tsv_filename=OUT_TSV_FILENAME, out_fasta_filename=OUT_FASTA_FILENAME,
        out_fig_filename=OUT_FIG_FILENAME):
    """ Download PRO ontology, generate proteoforms, and encode with BpForms

    Args:
        in_obo_filename (:obj:`str`, optional): path to save/read PRO ontology in OBO format
        in_pkl_filename (:obj:`str`, optional): path to save/read parsed content of PRO ontology
        in_tsv_filename (:obj:`str`, optional): path to read PRO entries in TSV format
        in_monomers_filename (:obj:`str`, optional): path to list of ids of monomeric forms used 
            by PRO and their alphabet code in tab-separated format
        max_num_proteins (:obj:`int`, optional): maximum number of proteins to analyze
        out_pickle_filename (:obj:`str`, optional): path to save results in pickle format
        out_pickle_filename_2 (:obj:`str`, optional): path to save results in pickle format
        out_tsv_filename (:obj:`str`, optional): path to save results in tab-separated format
        out_fasta_filename (:obj:`str`, optional): path to save results in FASTA format
        out_fig_filename (:obj:`str`, optional): path to save plot of results

    Returns:
        :obj:`list` of :obj:`dict`: proteoforms encoded with BpForms
    """
    # get the PRO ontology and extract the modified proteins from the ontology
    # proteins = get_pro_from_obo(obo_filename=in_obo_filename, pkl_filename=in_pkl_filename, max_num_proteins=max_num_proteins)
    proteins = get_pro_from_tsv(in_tsv_filename, max_num_proteins=max_num_proteins)

    # parse the modified proteins and retrieve their sequences
    if not os.path.isfile(out_pickle_filename):
        # parse the modified proteins and retrieve their sequences
        parsed_proteins = []
        for i_protein, protein in enumerate(proteins):
            if i_protein % 100 == 0:
                print('Parsing protein {} of {}'.format(i_protein + 1, len(proteins)))
            parsed_proteins.append(parse_protein(protein))

        # save the parsed proteins in pickle format
        with open(out_pickle_filename, 'wb') as file:
            pickle.dump(parsed_proteins, file)
    else:
        # load saved parsed proteins in pickle format
        with open(out_pickle_filename, 'rb') as file:
            parsed_proteins = pickle.load(file)

    # read list of monomers
    monomers = {}
    with open(in_monomers_filename, 'r') as file:
        reader = csv.DictReader(file, dialect='excel')
        for row in reader:
            monomers[row['PRO id']] = {
                'mod': bpforms.protein_alphabet.monomers.get(row['BpForms code'], None),
                'origin': [],
                }
            if row['Base monomer']:
                monomers[row['PRO id']]['origin'] = row['Base monomer'].split(', ')

    # generate list of modified monomeric forms
    for protein in parsed_proteins:
        for modification in protein['modifications']:
            if modification['monomer'] not in monomers:
                monomers[modification['monomer']] = {
                    'mod': None,
                    'origin': [],
                }

    # print list of unmapped monomers
    unmapped_monomers = []
    for monomer, code in monomers.items():
        if not code['mod']:
            unmapped_monomers.append(monomer)
    unmapped_monomers.sort()
    if unmapped_monomers:
        print('Several PRO monomeric forms have not been mapped to BpForms monomeric forms:\n  {}'.format(
            '\n  '.join(unmapped_monomers)))

    # check for inconsistencies between residue and modified monomeric form
    monomer_codes = {}
    for code, monomer in bpforms.protein_alphabet.monomers.items():
        monomer_codes[monomer] = code

    for protein in parsed_proteins:
        for modification in protein.get('modifications', []):
            if modification['residue'] and modification['monomer']:
                monomer = monomers.get(modification['monomer'], None)
                if (monomer['mod'] and monomer['mod'].get_canonical_code(monomer_codes) != modification['residue']) \
                    or (monomer['origin'] and modification['residue'] not in monomer['origin']):
                    codes = set(monomer['origin'])
                    if monomer['mod']:
                        codes.add(monomer['mod'].get_canonical_code(monomer_codes))                    
                    msg = 'Modified monomeric form {} potentially inconsistent with residue {} != {}'.format(
                        modification['monomer'], modification['residue'], 
                        ', '.join(codes))
                    print(protein['id'] + ': ' + msg)

    # generate BpForms for each protein
    if not os.path.isfile(out_pickle_filename_2):
        for i_protein, protein in enumerate(parsed_proteins):
            if i_protein % 100 == 0:
                print('Generating BpForms {} of {}'.format(i_protein + 1, len(parsed_proteins)))

            protein['modified_seq'] = None
            if not protein['uniprot_id']:
                continue
            if not protein['seq']:
                continue
            if protein['pro_errors']:
                continue

            processed_form = gen_bpform(protein, monomers, monomer_codes, apply_modifications=False)
            protein['processed_seq'] = str(processed_form)
            if not processed_form.validate():
                processed_formula = processed_form.get_formula()
                protein['processed_formula'] = str(processed_formula)
                protein['processed_mol_wt'] = processed_form.get_mol_wt()
                protein['processed_charge'] = processed_form.get_charge()

            if not protein['modifications']:
                continue

            modified_form = gen_bpform(protein, monomers, monomer_codes, include_annotations=False)
            protein['modified_seq'] = str(modified_form)

            modified_form = gen_bpform(protein, monomers, monomer_codes)
            if not modified_form.validate():
                modified_formula = modified_form.get_formula()
                protein['modified_formula'] = str(modified_formula)
                protein['modified_mol_wt'] = modified_form.get_mol_wt()
                protein['modified_charge'] = modified_form.get_charge()

                protein['modifications_formula'] = str(modified_formula - processed_formula)
                protein['modifications_mol_wt'] = protein['modified_mol_wt'] - protein['processed_mol_wt']
                protein['modifications_charge'] = protein['modified_charge'] - protein['processed_charge']

        # save the parsed proteins in pickle format
        with open(out_pickle_filename_2, 'wb') as file:
            pickle.dump(parsed_proteins, file)
    else:
        with open(out_pickle_filename_2, 'rb') as file:
            parsed_proteins = pickle.load(file)

    # save the proteoforms in TSV format
    with open(out_tsv_filename, 'w') as file:
        writer = csv.writer(file, dialect='excel-tab')
        writer.writerow(['PRO id', 'UniProt id', 'Organism',
                         'Unmodified sequence (IUBMB)',
                         'Processing', 'Processsed sequence (IUBMB)', 'Processsed formula', 'Processsed molecular weight', 'Processsed charge',
                         'Modifications', 'Modified sequence (BpForms)', 'Is modified sequence concrete', 'Modified formula', 'Modified molecular weight', 'Modified charge',
                         'Modifications formula', 'Modifications molecular weight', 'Modifications charge',
                         'PRO issues', 'Monomeric form issues'])

        for parsed_protein in parsed_proteins:
            if parsed_protein.get('pro_errors', None):
                pro_errors = '. '.join(parsed_protein['pro_errors']) + '.'
            else:
                pro_errors = None

            if parsed_protein.get('modified_errors', None):
                modified_errors = '. '.join(parsed_protein['modified_errors']) + '.'
            else:
                modified_errors = None

            writer.writerow([
                parsed_protein['id'],
                parsed_protein.get('uniprot_id', None),
                parsed_protein.get('organism', None),
                parsed_protein.get('seq', None),
                ', '.join('{}-{}'.format(p['start'], p['end']) for p in parsed_protein['processing']),
                parsed_protein.get('processed_seq', None),
                parsed_protein.get('processed_formula', None),
                parsed_protein.get('processed_mol_wt', None),
                parsed_protein.get('processed_charge', None),
                ', '.join('{} --> {} ({})'.format(m['residue'] or '?', m['monomer'], ', '.join(str(p) for p in m['positions']))
                          for m in parsed_protein['modifications']),
                parsed_protein.get('modified_seq', None),
                parsed_protein.get('modified_concrete', False),
                parsed_protein.get('modified_formula', None),
                parsed_protein.get('modified_mol_wt', None),
                parsed_protein.get('modified_charge', None),
                parsed_protein.get('modifications_formula', None),
                parsed_protein.get('modifications_mol_wt', None),
                parsed_protein.get('modifications_charge', None),
                pro_errors,
                modified_errors,
            ])

    # save the proteoforms in FASTA format
    seqs = (SeqRecord(id='{} | {}'.format(protein['id'], protein['uniprot_id']),
                      seq=Seq(protein['modified_seq']),
                      description='')
            for protein in parsed_proteins
            if protein['modified_seq'])
    SeqIO.write(seqs, out_fasta_filename, "fasta")

    # analyze frequency of modifications
    plot_modifications(parsed_proteins, fig_filename=out_fig_filename)

    # return proteins
    return proteins, parsed_proteins


def get_pro_from_obo(obo_filename=IN_OBO_FILENAME, pkl_filename=IN_PKL_FILENAME, max_num_proteins=None):
    """ Get the PRO ontology and extract the modified proteins from the ontology

    Args:
        obo_filename (:obj:`str`, optional): filename to save PRO in OBO format
        pkl_filename (:obj:`str`, optional): filename to save/read PRO from pickled file
        max_num_proteins (:obj:`int`, optional): maximum number of proteins to analyze

    Returns:
        :obj:`list` of :obj:`dict`: list of PRO ontology terms for modified proteins
    """
    # download PRO
    if not os.path.isfile(obo_filename):
        response = requests.get(IN_URL)
        response.raise_for_status()
        with open(obo_filename, 'wb') as file:
            file.write(response.content)

    # parse PRO or read from cache
    if not os.path.isfile(pkl_filename):
        # parse PRO
        proteins = []
        protein = None
        with open(obo_filename, 'r') as file:
            for line in file:
                line = line.rstrip('\n')
                if line.startswith('['):
                    if line.startswith('[Term]'):
                        if max_num_proteins is not None and len(proteins) >= max_num_proteins:
                            break
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

        if max_num_proteins is not None and max_num_proteins < len(proteins):
            proteins = proteins[0:max_num_proteins]

    # return PRO
    return proteins


def get_pro_from_tsv(filename, max_num_proteins=None):
    """ Extract PRO entries from TSV file

    Args:
        obo_filename (:obj:`str`, optional): filename to save PRO in OBO format
        max_num_proteins (:obj:`int`, optional): maximum number of proteins to analyze

    Returns:
        :obj:`list` of :obj:`dict`: list of PRO ontology terms for modified proteins
    """
    proteins = []
    with open(filename, 'r') as file:
        reader = csv.DictReader(file, fieldnames=('id', 'category', 'synonym_type', 'seq'), dialect='excel-tab')
        for row in reader:
            proteins.append({
                'id': [row['id']],
                'category': [row['category']],
                'synonym': ['"{}" {} PRO-proteoform-std'.format(row['seq'], row['synonym_type'])],
                })

            if max_num_proteins is not None and len(proteins) >= max_num_proteins:
                break

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
        if synonym.startswith('"UniProtKB:') and ' PRO-proteoform-std' in synonym:
            seq_synonyms.append(synonym)
    if not seq_synonyms:
        errors.append('No synonym which defines a modified sequence')
        return {
            'id': id,
            'uniprot_id': None,
            'processing': [],
            'modifications': [],
            'seq': None,
            'pro_errors': errors,
        }
    elif len(seq_synonyms) > 1:
        errors.append('Multiple synonyms which define modified sequences')
    synonym = seq_synonyms[0]

    uniprot_id, _, processing_modifications_type = synonym.partition(', ')
    uniprot_id = uniprot_id.partition(':')[2]

    organism_name = None
    response = session.get(UNIPROT_XML_ENDPOINT.format(uniprot_id))
    response.raise_for_status()
    if response.content:
        xml_root = ElementTree.fromstring(response.content)
        entry = xml_root.find('{http://uniprot.org/uniprot}entry')
        organism = entry.find('{http://uniprot.org/uniprot}organism')
        names = organism.findall('{http://uniprot.org/uniprot}name')
        for name in names:
            if name.get('type') == 'scientific':
                organism_name = name.text
                break

    response = session.get(UNIPROT_SEQ_ENDPOINT.format(uniprot_id))
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
    if processing_modifications.startswith('which') \
        or processing_modifications.startswith('with') \
        or 'MOD:00046 OR Thr-163, MOD:00047' in processing_modifications:
        modifications_str = []
        errors.append('Unable to parse sequence')
    elif processing_modifications:
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

            residue_codes = set()
            positions = []
            for residue_position in residue_positions.split('/'):
                residue_chars, _, position = residue_position.partition('-')
                residue_code = AA_CHARS_TO_CODES[residue_chars]
                position = int(float(position))
                if position > len(seq):
                    errors.append('Position {} is greater than the sequence length {}'.format(position, len(seq))) 
                elif seq[position - 1] != residue_code:
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
        'organism': organism_name,
        'processing': processing,
        'modifications': modifications,
        'seq': seq,
        'pro_errors': errors,
    }


def gen_bpform(protein, pro_ids_to_bpform_monomers, monomer_codes, apply_modifications=True, include_annotations=True):
    """ Generate BpForm for a modified protein in PRO

    Args:
        protein (:obj:`dict`): term for modified protein
        pro_ids_to_bpform_monomers (:obj:`dict`): dictionary which maps ids of monomeric forms
            used by PRO to monomeric forms in the BpForms protein alphabet
        monomer_codes (:obj:`dict`): dictionary that maps monomers to their codes in the alphabet
        apply_modifications (:obj:`bool`, optional): if :obj:`True`, include modifications in proteoform
        include_annotations (:obj:`bool`, optional): if :obj:`True`, include metadata about modified monomers

    Returns:
        :obj:`bpforms.ProteinForm`: BpForm for a term in PRO
    """
    form = bpforms.ProteinForm()
    monomers = bpforms.protein_alphabet.monomers

    # generate BpForm for unmodified sequence
    for base in protein['seq']:
        form.seq.append(monomers[base])

    # apply processing
    modifications = copy.deepcopy(protein['modifications'])
    seq = protein['seq']
    if protein['processing']:
        procesed_seq = []
        seq = ''
        for processing in protein['processing']:
            procesed_seq.extend(form.seq[processing['start']-1:processing['end']])
            seq += protein['seq'][processing['start']-1:processing['end']]
        form.seq = procesed_seq

        for modification in modifications:
            modification['processed_positions'] = []
            for position in modification['positions']:
                seq_len = 0
                processed_position = None
                for processing in protein['processing']:
                    if position >= processing['start'] and position <= processing['end']:
                        processed_position = seq_len + position - processing['start'] + 1
                        break
                    seq_len += processing['end'] - processing['start'] + 1
                if processed_position is not None:
                    modification['processed_positions'].append(processed_position)
    else:
        for modification in modifications:
            modification['processed_positions'] = modification['positions']

    # apply modifications
    if apply_modifications:
        concrete = True
        protein['modified_errors'] = []

        for modification in modifications:
            monomer = pro_ids_to_bpform_monomers[modification['monomer']]['mod']
            origin = pro_ids_to_bpform_monomers[modification['monomer']]['origin']
            if modification['monomer'] == 'PR:000026291':
                if include_annotations:
                    monomer = bpforms.Monomer().from_dict(
                        monomers[modification['residue']].to_dict(
                        alphabet=bpforms.protein_alphabet),
                        alphabet=bpforms.protein_alphabet)
                    monomer.base_monomers.add(monomers[modification['residue']])
                else:
                    monomer = bpforms.Monomer(id=modification['residue'])
                monomer.name = None
                monomer.synonyms = []
                monomer.identifiers = [bpforms.Identifier('pro', modification['monomer'])]
                monomer.comments = None

            elif monomer is None:
                concrete = False
                monomer = bpforms.Monomer(id=modification['monomer'])
                if modification['residue']:
                    monomer.base_monomers.add(monomers[modification['residue']])
                else:
                    monomer.base_monomers.update(monomers[base] for base in origin)

            if modification['positions']:
                for position in modification['processed_positions']:
                    if form.seq[position - 1] == monomers[seq[position - 1]]:
                        form.seq[position - 1] = monomer
                    else:
                        protein['modified_errors'].append('Unable to set monomeric form at position {}'.format(
                            position))

            elif modification['residue']:
                concrete = False

                if pro_ids_to_bpform_monomers[modification['monomer']]['mod'] is None:
                    base_monomers = monomer.base_monomers
                else:
                    base_monomers = set()
                monomer = bpforms.Monomer(id=monomer.id, base_monomers=base_monomers)

                monomer.start_position = seq.find(modification['residue']) + 1
                monomer.end_position = seq.rfind(modification['residue']) + 1
                set_monomer = False
                for i_monomer in range(monomer.start_position, monomer.end_position + 1):
                    if form.seq[i_monomer - 1] == monomers[seq[i_monomer - 1]]:
                        set_monomer = True
                        form.seq[i_monomer - 1] = monomer
                        break
                if not set_monomer:
                    protein['modified_errors'].append('Unable to set monomeric form')
            else:
                concrete = False

                canonical_code = monomer.get_canonical_code(monomer_codes)
                if pro_ids_to_bpform_monomers[modification['monomer']]['mod'] is None:
                    base_monomers = monomer.base_monomers
                else:
                    base_monomers = set()
                monomer = bpforms.Monomer(id=monomer.id, base_monomers=base_monomers)

                if canonical_code and canonical_code != '?':
                    start_position = seq.find(canonical_code) + 1
                    end_position = seq.rfind(canonical_code) + 1
                    if start_position == 0:
                        protein['modified_errors'].append('Sequence does not contain residue {} for modification {}'.format(
                            canonical_code, monomer.id))
                    else:
                        monomer.start_position = start_position
                        monomer.end_position = end_position

                elif origin:
                    start_position = float('inf')
                    end_position = -float('inf')
                    for base in origin:
                        start_pos = seq.find(base) + 1
                        if start_pos > 0:
                              start_position = min(start_position, start_pos)

                        end_pos = seq.rfind(base) + 1
                        if end_pos > 0:
                            end_position = max(end_position, end_pos)

                    if numpy.isinf(start_position):
                        protein['modified_errors'].append('Sequence does not contain residues {} for modification {}'.format(
                            ', '.join(origin), monomer.id))
                    else:
                        monomer.start_position = start_position
                        monomer.end_position = end_position

                else:
                    monomer.start_position = 1
                    monomer.end_position = len(seq)

                if monomer.start_position:
                    set_monomer = False
                    for i_monomer in range(monomer.start_position, monomer.end_position + 1):
                        if form.seq[i_monomer - 1] == monomers[seq[i_monomer - 1]]:
                            form.seq[i_monomer - 1] = monomer
                            set_monomer = True
                            break
                    if not set_monomer:
                        protein['modified_errors'].append('Unable to set monomeric form')

    # validate
    if apply_modifications:
        protein['modified_concrete'] = concrete
        protein['modified_errors'].extend(form.validate())

    # return proteoform represented with BpForms
    return form


def plot_modifications(proteins, organism='Homo sapiens', fig_filename=OUT_FIG_FILENAME):
    """ Plot a summary of the modifications in PRO

    Args:
        proteins (:obj:`list` of :obj:`dict`): entries in PRO ontology
        organism (:obj:`str`, optional): organism to analyze
        fig_filename (:obj:`str`, optional): path to save analysis
    """
    code_freq = {}
    canonical_code_freq = {}
    for protein in proteins:
        if (organism is None or protein.get('organism', None) == organism) and protein.get('modified_seq', None):
            for modification in protein['modifications']:
                if modification['residue'] and modification['monomer']:
                    n_mods = max(1, len(modification['positions']))
                    if modification['residue'] not in canonical_code_freq:
                        canonical_code_freq[modification['residue']] = 0
                    if modification['monomer'] not in code_freq:
                        code_freq[modification['monomer']] = 0
                    canonical_code_freq[modification['residue']] += n_mods
                    code_freq[modification['monomer']] += n_mods

    pyplot.style.use('ggplot')
    fig, axes = pyplot.subplots(nrows=1, ncols=2, gridspec_kw={'width_ratios': [1, 4]})
    fig.set_size_inches(9.3, 1.5)
    plot_codes(canonical_code_freq,
               'Frequency of modifications',
               axes[0], ignore_canonical=False)
    plot_codes(code_freq,
               'Frequency of modified monomeric forms',
               axes[1], ignore_canonical=True)
    fig.savefig(fig_filename, transparent=True,
                bbox_inches=matplotlib.transforms.Bbox([[0.69, -0.5], [8.35, 1.5]]))
    pyplot.close(fig)


def plot_codes(code_freq, title, axis, ignore_canonical=False):
    id_freqs = []
    for code, count in code_freq.items():
        if ignore_canonical and code in ['A', 'C', 'G', 'U']:
            continue
        id_freqs.append((code, count))
    id_freqs.sort()

    y_pos = numpy.arange(len(id_freqs))
    freq = numpy.array([id_freq[-1] for id_freq in id_freqs])
    freq = freq / numpy.sum(freq) * 100.
    x_tick_labels = {id: y_pos for y_pos, (id, _) in enumerate(id_freqs)}

    axis.bar(y_pos, freq, align='center', alpha=0.5)
    axis.set_xticks(y_pos)
    axis.set_xticklabels(x_tick_labels, rotation=270, fontsize=6, fontfamily='Raleway')
    axis.set_ylabel('Frequency (%)', fontdict={
        'fontsize': 10,
        'fontweight': 'regular',
        'fontfamily': 'Raleway',
    })
    axis.set_title(title, fontdict={
        'fontsize': 10,
        'fontweight': 'regular',
        'fontfamily': 'Raleway',
    })
    axis.set_xlim((-0.75, len(id_freqs) - 0.25))
