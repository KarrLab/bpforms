from modomics import plot_codes
from matplotlib import pyplot
import bpforms
import csv
import matplotlib
import os
import sys
sys.path.insert(0, 'examples')

canonical_monomers = [
    bpforms.rna_alphabet.monomers.A,
    bpforms.rna_alphabet.monomers.C,
    bpforms.rna_alphabet.monomers.G,
    bpforms.rna_alphabet.monomers.U,
]

tot_rna = 375000
doubling_time = 45.  # min
half_life = 45.  # min

monomer_freq = {}
can_monomer_freq = {}
monomer_codes = {monomer: code for code, monomer in bpforms.rna_alphabet.monomers.items()}
tot_copies = 0

with open('examples/modomics_trna_copy_numbers.csv', 'r') as file:
    for rna in csv.DictReader(file, dialect='excel'):
        form = bpforms.RnaForm().from_str(rna['Sequence (BpForms)'])
        copies = float(rna['Copies per cell'])

        tot_copies += copies

        for monomer in form.seq:
            if monomer not in canonical_monomers:
                if monomer not in monomer_freq:
                    monomer_freq[monomer] = 0
                monomer_freq[monomer] += copies

                can_code = monomer.get_canonical_code(monomer_codes)
                if can_code not in can_monomer_freq:
                    can_monomer_freq[can_code] = 0
                can_monomer_freq[can_code] += copies

for code in monomer_freq.keys():
    monomer_freq[code] *= (1. + doubling_time / half_life)
for monomer in can_monomer_freq.keys():
    can_monomer_freq[monomer] *= (1. + doubling_time / half_life)

tot_mods = sum(freq for monomer, freq in monomer_freq.items()) * tot_rna / tot_copies

pyplot.style.use('default')
fig, axes = pyplot.subplots(nrows=1, ncols=2, gridspec_kw={'width_ratios': [2.717, 16]})
fig.set_size_inches(6.5, 1.5)

axis1, bars = plot_codes(can_monomer_freq,
                         bpforms.rna_alphabet.monomers,
                         axes[0], ignore_canonical=False, use_rel=False,
                         x_axis_label='Canonical residue', y_axis_label=r'Freq ($10^5$ nt cell cycle$^{-1}$)',
                         axis_font_size=7, tick_font_size=6, font_family='Arial', x_label_pad=16.7)
bars[0].set_color('#e74624')
bars[1].set_color('#2daae1')
bars[2].set_color('#90e227')
bars[3].set_color('#dabe2e')

monomer_code_freq = {monomer_codes[monomer]: freq for monomer, freq in monomer_freq.items()}
axis2, bars = plot_codes(monomer_code_freq,
                         bpforms.rna_alphabet.monomers,
                         axes[1], ignore_canonical=True, use_rel=False,
                         x_axis_label='Modified residue', y_axis_label=r'Freq ($10^5$ nt cell cycle$^{-1}$)',
                         axis_font_size=7, tick_font_size=6, font_family='Arial')
for label, bar in zip(axis2.get_xticklabels(), bars):
    code = label.get_text()
    if code[-1] == 'A':
        bar.set_color('#e74624')
    elif code[-1] == 'C':
        bar.set_color('#2daae1')
    elif code[-1] == 'G':
        bar.set_color('#90e227')
    elif code[-1] == 'U':
        bar.set_color('#dabe2e')

fig.savefig(os.path.join('examples', 'modomics.trna.analysis.svg'),
            transparent=True,
            bbox_inches=matplotlib.transforms.Bbox([[0.45, -0.31], [5.86, 1.35]]))
fig.savefig(os.path.join('examples', 'modomics.trna.analysis.pdf'),
            transparent=True,
            bbox_inches=matplotlib.transforms.Bbox([[0.45, -0.31], [5.86, 1.35]]))
pyplot.close(fig)

top_five = (monomer_code_freq['8U'] +
            monomer_code_freq['9U'] +
            monomer_code_freq['5U'] +
            monomer_code_freq['74U'] +
            monomer_code_freq['7G']) / \
sum(monomer_code_freq.values())
