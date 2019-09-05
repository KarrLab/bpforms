""" Curate modified rRNA and tRNA sequences

1. Determine genomic coordinates and non-canonical sequences
   - rRNA: RNAcentral
   - tRNA: Gogakos et al., 2017
2. Curated modifications from MODOMICS
3. Calculated non-canonical sequences for each modified sequence from MODOMICS
4. Aligned genomically-anchored non-canonical sequences and the MODOMICS non-canonical sequences
   rRNA: I computed the alignments with Clustal Omega
   tRNA: I computed the alignments with BioPython pairwise2
5. I lifted over the modifications from the MODOMICS sequences to the genomically-achored sequences
"""

from Bio import pairwise2
from Bio import SeqIO
from Bio.Alphabet import SingleLetterAlphabet
from Bio.Seq import Seq
import bpforms
import pandas

can_monomers = [bpforms.rna_alphabet.monomers.A, 
                bpforms.rna_alphabet.monomers.C,
                bpforms.rna_alphabet.monomers.G,
                bpforms.rna_alphabet.monomers.U]

def liftover_mods(ref_seq, ref_nc_seq, seq, can_monomers=can_monomers):
    ref_form = bpforms.RnaForm().from_str(ref_nc_seq)
    form = bpforms.RnaForm()
    i_nc_nt = 0
    for i_nt, ref_monomer in enumerate(ref_seq):
        if i_nc_nt < len(ref_form.seq):
            ref_nc_monomer = ref_form.seq[i_nc_nt]
        else:
            ref_nc_monomer = None
        if ref_monomer != '-':
            i_nc_nt += 1

        monomer = seq[i_nt]
        if monomer == '-':
            continue
        elif ref_nc_monomer not in can_monomers and monomer == ref_monomer:
            form.seq.append(ref_nc_monomer)
        elif monomer == 'N':
            form.seq.append(bpforms.Monomer(id='N'))
        else:
            form.seq.append(bpforms.rna_alphabet.monomers.get(monomer))

    # verify non-canonical sequences are consistent with the canonical sequences
    assert seq.replace('-', '') == form.get_canonical_seq()

    return form

'''
rRNA
'''
rrna_types = ['5.8S', '18S', '28S']

for rrna_type in rrna_types:
    filename = 'examples/homo_sapiens_rna/{} all seqs.fasta'.format(rrna_type)
    seqs = [str(record.seq) for record in SeqIO.parse(filename, "fasta")]
    ref_nc_seq = seqs[0]
    ref_form = bpforms.RnaForm().from_str(ref_nc_seq)
    ref_seq = seqs[1]
    seqs = seqs[1:]

    # map curated modifications onto sequences
    forms = []
    for seq in seqs:
        forms.append(liftover_mods(ref_seq, ref_nc_seq, seq))

    # save non-canonical sequences
    with open('examples/homo_sapiens_rna/{} nc alignment.txt'.format(rrna_type), 'w') as file:
        for form in forms:
            file.write(str(form) + '\n')

'''
tRNA
'''
trna_mod_rows = pandas.read_excel('examples/homo_sapiens_rna/summary.xlsx', sheet_name='tRNA modifications - MODOMICS')
trna_mods = {}
for _, trna_mod_row in trna_mod_rows.iterrows():
    if trna_mod_row['Organellum'] != 'cytosolic':
        continue
    key = (trna_mod_row['Amino acid'], trna_mod_row['Anticodon (Canonical)'])
    if key not in trna_mods:
        trna_mods[key] = []
    trna_mods[key].append({
        'id': trna_mod_row['Id'],
        'can': trna_mod_row['Sequence (Canonical)'],
        'nc': trna_mod_row['Sequence'],
        })

trna_seq_rows = pandas.read_excel('examples/homo_sapiens_rna/summary.xlsx', 
    sheet_name='tRNA seqs - Gogakos et al.', header=[0, 1])
trna_seqs = []
alphabet = SingleLetterAlphabet()
for _, trna_seq_row in trna_seq_rows.iterrows():
    aa = list(trna_seq_row.items())[0][1]
    anticodon = list(trna_seq_row.items())[1][1]
    key = (aa, anticodon)
    if key not in trna_mods:
        continue

    # find most similar sequences in MODOMICS
    best_id = None
    best_can_seq = None
    best_seq = None
    best_alignment = None
    best_score = -float('inf')
    for trna_mod in trna_mods[key]:
        alignment = pairwise2.align.globalxs(
            Seq(trna_mod['can'], alphabet=alphabet), 
            Seq(trna_seq_row['Mature tRNA', 'Sequence'].replace('T', 'U'), alphabet=alphabet),
            -10, -0.5, one_alignment_only=True)[0]
        score = alignment[2]
        if score > best_score:
            best_id = trna_mod['id']
            best_can_seq = trna_mod['can']
            best_seq = trna_mod['nc']
            best_alignment = alignment
            best_score = score

    # liftover modifications
    form = liftover_mods(best_alignment[0], best_seq, best_alignment[1])
    trna_seqs.append({
        'id': trna_seq_row['Gene/Pre-tRNA', 'Id'],
        'modomics_id': best_id,
        'seq': str(form),
        })

with open('examples/homo_sapiens_rna/trna nc.tsv', 'w') as file:
    file.write('{}\t{}\t{}\n'.format('Pre-tRNA id', 'MODOMICS id', 'Sequence'))
    for trna_seq in trna_seqs:
        file.write('{}\t{}\t{}\n'.format(
            trna_seq['id'], trna_seq['modomics_id'], trna_seq['seq']))
