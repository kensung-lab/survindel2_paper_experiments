from __future__ import print_function

import numpy as np
import vcf, argparse
from collections import defaultdict
from intervaltree import IntervalTree
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

cmd_parser = argparse.ArgumentParser(description='Given comma-separated VCF files and a TRF file, '
                                                 'plot the number of copies being deleted or duplicated.')
cmd_parser.add_argument('vcf_files')
cmd_parser.add_argument('trf_file')
cmd_parser.add_argument('outfile')
cmd_args = cmd_parser.parse_args()

def get_svlen(sv):
    svlen = sv.INFO["SVLEN"]
    if isinstance(svlen, int):
        return abs(svlen)
    else:
        return abs(svlen[0])

trf_ivs = defaultdict(IntervalTree)
with open(cmd_args.trf_file) as trf_f:
    for line in trf_f:
        if line.startswith("#"):
            continue
        sl = line.split()
        chr, start, end, period, perc_match = sl[1], int(sl[2]), int(sl[3]), int(sl[5]), int(sl[8])

        trf_ivs[chr].addi(start, end, (period, perc_match))


n_copies_deldup = list()
vcf_filenames = cmd_args.vcf_files.split(",")
for vcf_filename in vcf_filenames:
    vcf_reader = vcf.Reader(filename=vcf_filename)
    for sv in vcf_reader:
        if sv.FILTER: continue

        start, end = sv.start + 1, sv.end

        reps = trf_ivs[sv.CHROM].overlap(start, end+1)
        svlen = get_svlen(sv)
        min_period = 1000000
        for rep in reps:
            if min_period > rep.data[0]:
                min_period = rep.data[0]
        n_copies_deldup.append(svlen/min_period)

upper_bound = 30
n_bins = upper_bound*10+1

linbins = np.linspace(0, upper_bound, n_bins)
hist, bin_edges = np.histogram(n_copies_deldup, bins=linbins)
colors = ["red" if i%10 == 0 else "blue" for i in range(n_bins-1)]
plt.bar(range(n_bins-1), hist, color=colors)
plt.xscale('linear')
plt.xlabel('Number of copies deleted or duplicated')
plt.ylabel('Count')
plt.legend(handles=[mpatches.Patch(color='red', label='Integer number of copies'),
                    mpatches.Patch(color='blue', label='Non-integer number of copies')])
ticks_pos = range(0, n_bins, 50)
plt.xticks(ticks_pos, [str(int(i/10)) for i in ticks_pos])
plt.savefig(cmd_args.outfile, bbox_inches='tight')
plt.clf()
plt.close()
