""" Calculate properties for BpForms modeled in Bouhaddou et al., PLoS Comput Biol, 2018.

This example uses both the Python API and the JSON REST API

Bouhaddou M, Barrette AM, Stern AD, Koch RJ, DiStefano MS, Riesel EA, Santos LC, Tan AL, Mertz AE & Birtwistle MR.
A mechanistic pan-cancer pathway model informed by multi-omics data interprets stochastic cell fate responses to drugs and mitogens.
PLoS Comput Biol 2018, 14(3): e1005985. doi: `10.1371/journal.pcbi.1005985 <http://dx.plos.org/10.1371/journal.pcbi.1005985>`_.

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Marc Birtwistle <mbirtwi@clemson.edu>
:Author: Cemal Erdem <cemale@clemson.edu>
:Date: 2019-06-18
:Copyright: 2019, Karr Lab
:License: MIT
"""

from Bio import SeqIO
import bpforms
import csv
import os.path
import requests

IN_FILENAME = os.path.join('examples', 'bouhaddou_et_al_plos_comput_biol_2018.fasta')
OUT_FILENAME = os.path.join('examples', 'bouhaddou_et_al_plos_comput_biol_2018.tsv')

ENDPOINT = 'https://www.bpforms.org'


def run():
    # read BpForms from FASTA file
    seqs = []
    with open(IN_FILENAME, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            seqs.append({'id': record.id, 'seq': str(record.seq)})

    seq_props = calc_bpforms_props_with_python_api(seqs)
    seq_props = calc_bpforms_props_with_rest_api(seqs)

    # save computed properties to .tsv file
    with open(OUT_FILENAME, 'w') as file:
        writer = csv.DictWriter(file, fieldnames=['Species', 'Formula', 'Molecular weight', 'Charge', 'Length'], dialect='excel-tab')
        writer.writeheader()
        for seq_prop in seq_props:
            writer.writerow(seq_prop)


def calc_bpforms_props_with_python_api(seqs):
    # calculate properties
    seq_props = []
    for seq in seqs:
        form = bpforms.ProteinForm().from_str(seq['seq'])
        seq_props.append({
            'Species': seq['id'],
            'Formula': form.get_formula(),
            'Molecular weight': form.get_mol_wt(),
            'Charge': form.get_charge(),
            'Length': len(form.seq),
        })
    return seq_props


def calc_bpforms_props_with_rest_api(seqs):
    seq_props = []
    for seq in seqs:
        data = {
            "alphabet": "protein",
            "seq": seq['seq'],
            "circular": False,
            "major_tautomer": False,
            "dearomatize": False,
        }

        response = requests.post(ENDPOINT + '/api/bpform/', json=data)
        response.raise_for_status()
        props = response.json()

        seq_props.append({
            'Species': seq['id'],
            'Formula': props['formula'],
            'Molecular weight': props['mol_wt'],
            'Charge': props['charge'],
            'Length': props['length'],
        })
    return seq_props
