""" Generate BpForms for all of the tRNA in MODOMICS, verify
them, and calculate their properties

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-20
:Copyright: 2019, Karr Lab
:License: MIT
"""

from matplotlib import pyplot
import bpforms
import bs4
import csv
import matplotlib
import numpy
import os
import requests
import requests_cache

URL = 'http://modomics.genesilico.pl/sequences/list/'


def run():
    # create dict of MODOMICS single character monomer codes
    monomer_codes = {}
    for monomer in bpforms.rna_alphabet.monomers.values():
        for synonym in monomer.synonyms:
            if len(synonym) == 1:
                monomer_codes[synonym.upper()] = monomer
                break

    # create cache for web queries
    cache_name = os.path.join('examples', 'modomics')
    session = requests_cache.core.CachedSession(cache_name, backend='sqlite', expire_after=None)
    session.mount('http://modomics.genesilico.pl', requests.adapters.HTTPAdapter(max_retries=5))

    # parse rRNA and tRNA data
    ssu_rna_forms = run_rrna(session, monomer_codes, 'SSU', os.path.join('examples', 'modomics.ssu-rrna.tsv'))
    lsu_rna_forms = run_rrna(session, monomer_codes, 'LSU', os.path.join('examples', 'modomics.lsu-rrna.tsv'))
    trna_forms, trna_canonical_code_freq, trna_code_freq = run_trna(session, monomer_codes, os.path.join('examples', 'modomics.trna.tsv'))

    # plot distribution of tRNA monomeric forms
    plot_trna_code_freq(monomer_codes, trna_code_freq)

    # return results
    return ssu_rna_forms, lsu_rna_forms, trna_forms, trna_canonical_code_freq, trna_code_freq


def run_rrna(session, monomer_codes, type, out_filename):
    response = session.get(URL + type + '/')
    response.raise_for_status()

    doc = bs4.BeautifulSoup(response.text, 'lxml')
    table = doc.find('table', {'class': 'tabseq'})
    tbody = table.find('tbody')
    rows = tbody.find_all('tr')
    rna_forms = []
    for row in rows:
        if not isinstance(row, bs4.element.Tag):
            continue

        cells = row.find_all('td')

        rna_form = bpforms.RnaForm()
        unsupported_codes = set()
        for code in cells[5].text.strip().replace('-', '').replace('_', '').upper():
            monomer = monomer_codes.get(code, None)
            if monomer is None:
                unsupported_codes.add(code)
                monomer = bpforms.Monomer(id=code)
            rna_form.seq.append(monomer)

        rna_forms.append({
            'GenBank': cells[2].text,
            'Organism': cells[3].text,
            'Organellum': cells[4].text,
            'Sequence (MODOMICS)': cells[5].text.strip().replace('-', '').replace('_', ''),
        })
        analyze_form(rna_form, unsupported_codes, rna_forms[-1])

    # save results to tab-separated file
    save_results(rna_forms, ['GenBank'], out_filename)

    return rna_forms


def run_trna(session, monomer_codes, out_filename):
    response = session.get(URL + 'tRNA/')
    response.raise_for_status()

    text = response.text \
        .replace('<td><AU</td>', '<td>&lt;AU</td>') \
        .replace('<td><AA</td>', '<td>&lt;AA</td>')
    doc = bs4.BeautifulSoup(text, 'lxml')
    table = doc.find('table', {'class': 'tabseq'})
    tbody = table.find('tbody')
    rows = tbody.find_all('tr')
    rna_forms = []

    code_freq = {}
    canonical_code_freq = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
    for row in rows:
        if not isinstance(row, bs4.element.Tag):
            continue

        cells = row.find_all('td')

        rna_form = bpforms.RnaForm()
        unsupported_codes = set()
        for code in cells[5].text.strip().replace('-', '').replace('_', '').upper():
            monomer = monomer_codes.get(code, None)
            if monomer is None:
                unsupported_codes.add(code)
                monomer = bpforms.Monomer(id=code)
            else:
                if code not in code_freq:
                    code_freq[code] = 0
                code_freq[code] += 1
            rna_form.seq.append(monomer)

        rna_forms.append({
            'Amino acid type': cells[1].text,
            'Anticodon': cells[2].text,
            'Organism': cells[3].text,
            'Organellum': cells[4].text,
            'Sequence (MODOMICS)': cells[5].text.strip().replace('-', '').replace('_', ''),
        })
        analyze_form(rna_form, unsupported_codes, rna_forms[-1])

        canonical_code_freq['A'] += \
            rna_forms[-1]['Sequence (IUPAC)'].count('A') \
            - rna_form.seq.count(bpforms.rna_alphabet.monomers.A)
        canonical_code_freq['C'] += \
            rna_forms[-1]['Sequence (IUPAC)'].count('C') \
            - rna_form.seq.count(bpforms.rna_alphabet.monomers.C)
        canonical_code_freq['G'] += \
            rna_forms[-1]['Sequence (IUPAC)'].count('G') \
            - rna_form.seq.count(bpforms.rna_alphabet.monomers.G)
        canonical_code_freq['U'] += \
            rna_forms[-1]['Sequence (IUPAC)'].count('U') \
            - rna_form.seq.count(bpforms.rna_alphabet.monomers.U)

    # save results to tab-separated file
    save_results(rna_forms, ['Amino acid type', 'Anticodon'], out_filename)

    with open(os.path.join('examples', 'modomics.trna.canonical-code-freq.tsv'), 'w') as file:
        writer = csv.DictWriter(file, fieldnames=['Code', 'Frequency'], dialect='excel-tab')
        writer.writeheader()
        for code, freq in canonical_code_freq.items():
            writer.writerow({'Code': code, 'Frequency': freq})

    with open(os.path.join('examples', 'modomics.trna.code-freq.tsv'), 'w') as file:
        writer = csv.DictWriter(file, fieldnames=['Code', 'Frequency'], dialect='excel-tab')
        writer.writeheader()
        for code, freq in code_freq.items():
            writer.writerow({'Code': code, 'Frequency': freq})

    return rna_forms, canonical_code_freq, code_freq


def analyze_form(rna_form, unsupported_codes, results_dict):
    results_dict['Sequence (BpForms)'] = str(rna_form)
    results_dict['Sequence (IUPAC)'] = canonical_seq = rna_form.get_canonical_seq()
    results_dict['Length'] = len(rna_form.seq)

    results_dict['Number of modifications'] = len(rna_form.seq) \
        - rna_form.seq.count(bpforms.rna_alphabet.monomers.A) \
        - rna_form.seq.count(bpforms.rna_alphabet.monomers.C) \
        - rna_form.seq.count(bpforms.rna_alphabet.monomers.G) \
        - rna_form.seq.count(bpforms.rna_alphabet.monomers.U)
    results_dict['Number of modified A'] = canonical_seq.count('A') - rna_form.seq.count(bpforms.rna_alphabet.monomers.A)
    results_dict['Number of modified C'] = canonical_seq.count('C') - rna_form.seq.count(bpforms.rna_alphabet.monomers.C)
    results_dict['Number of modified G'] = canonical_seq.count('G') - rna_form.seq.count(bpforms.rna_alphabet.monomers.G)
    results_dict['Number of modified U'] = canonical_seq.count('U') - rna_form.seq.count(bpforms.rna_alphabet.monomers.U)

    if unsupported_codes:
        results_dict['BpForms errors'] = 'MODOMICS sequence uses monomeric forms {}'.format(
            ', '.join(unsupported_codes))
    else:
        results_dict['Formula'] = str(rna_form.get_formula())
        results_dict['Molecular weight'] = rna_form.get_mol_wt()
        results_dict['Charge'] = rna_form.get_charge()

        canonical_form = bpforms.RnaForm().from_str(canonical_seq)
        results_dict['Canonical formula'] = str(canonical_form.get_formula())
        results_dict['Canonical molecular weight'] = canonical_form.get_mol_wt()
        results_dict['Canonical charge'] = canonical_form.get_charge()

        results_dict['Extra formula'] = str(rna_form.get_formula() - canonical_form.get_formula())
        results_dict['Extra molecular weight'] = rna_form.get_mol_wt() - canonical_form.get_mol_wt()
        results_dict['Extra charge'] = rna_form.get_charge() - canonical_form.get_charge()

        results_dict['BpForms errors'] = ' '.join(rna_form.validate())


def save_results(rna_forms, variable_field_names, out_filename):
    with open(out_filename, 'w') as file:
        writer = csv.DictWriter(file,
                                fieldnames=variable_field_names + [
                                    'Organism', 'Organellum',
                                    'Sequence (MODOMICS)',
                                    'Sequence (BpForms)',
                                    'Sequence (IUPAC)',
                                    'Length',
                                    'Number of modifications',
                                    'Number of modified A',
                                    'Number of modified C',
                                    'Number of modified G',
                                    'Number of modified U',
                                    'Formula', 'Molecular weight', 'Charge',
                                    'Canonical formula', 'Canonical molecular weight', 'Canonical charge',
                                    'Extra formula', 'Extra molecular weight', 'Extra charge',
                                    'BpForms errors'],
                                dialect='excel-tab')
        writer.writeheader()

        for rna_form in rna_forms:
            writer.writerow(rna_form)


def plot_trna_code_freq(monomer_codes, trna_code_freq):
    monomer_ids = {}
    for code, monomer in bpforms.rna_alphabet.monomers.items():
        monomer_ids[monomer] = code

    trna_canonical_code_freq = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
    for code, count in trna_code_freq.items():
        if code in ['A', 'C', 'G', 'U']:
            continue
        canonical_code = monomer_codes[code].get_canonical_code(monomer_ids)
        trna_canonical_code_freq[canonical_code] += count

    pyplot.style.use('ggplot')
    fig, axes = pyplot.subplots(nrows=1, ncols=2, gridspec_kw={'width_ratios': [1, 16]})
    fig.set_size_inches(9.3, 1.5)
    plot_codes(trna_canonical_code_freq, monomer_codes,
               'Frequency of modifications',
               axes[0], ignore_canonical=False)
    plot_codes(trna_code_freq, monomer_codes,
               'Frequency of modified monomeric forms',
               axes[1], ignore_canonical=True)
    fig.savefig(os.path.join('bpforms', 'web', 'img', 'example', 'modomics.trna.code-freq.svg'),
                transparent=True,
                bbox_inches=matplotlib.transforms.Bbox([[0.4, -0.1], [8.35, 1.5]]))
    pyplot.close(fig)


def plot_codes(code_freq, monomer_codes, title, axis, ignore_canonical=False):
    monomer_ids = {}
    for code, monomer in bpforms.rna_alphabet.monomers.items():
        monomer_ids[monomer] = code

    id_freqs = []
    for code, count in code_freq.items():
        if ignore_canonical and code in ['A', 'C', 'G', 'U']:
            continue
        id_freqs.append((monomer_codes[code].get_canonical_code(monomer_ids),
                         monomer_codes[code].id,
                         count))
    id_freqs.sort()

    y_pos = numpy.arange(len(id_freqs))
    freq = numpy.array([id_freq[2] for id_freq in id_freqs])
    freq = freq / numpy.sum(freq) * 100.
    x_tick_labels = {id: y_pos for y_pos, (_, id, _) in enumerate(id_freqs)}

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
