from __future__ import print_function
import argparse
import statistics
from collections import defaultdict

import numpy as np
from sklearn.decomposition import PCA
import allel
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

cmd_parser = argparse.ArgumentParser(description='Analyse clustered SVs.')
cmd_parser.add_argument('del_vcf', help='VCF file with deletions.')
cmd_parser.add_argument('dup_vcf', help='VCF file with duplications.')
cmd_parser.add_argument('metadata', help='File with metadata.')
cmd_parser.add_argument('outdir', help='Directory where to output the figures.')
cmd_args = cmd_parser.parse_args()

colors = ['red', 'blue', 'green', 'pink', 'orange', 'purple', 'brown', 'yellow', 'black', 'gray',
          'red', 'blue', 'green', 'pink', 'orange', 'purple', 'brown', 'yellow', 'black', 'gray',
          'red', 'blue', 'green', 'pink', 'orange', 'purple', 'brown', 'yellow', 'black', 'gray',
          'red', 'blue', 'green', 'pink', 'orange', 'purple', 'brown', 'yellow', 'black', 'gray']

class md_entry_t:
    def __init__(self, cohort, race):
        self.race = race
        self.cohort = cohort

# Input
metadata = dict()
with open(cmd_args.metadata, encoding="ISO-8859-1") as metadata_fin:
    for line in metadata_fin.readlines()[1:]:
        sl = line.rstrip().split('\t')
        id, cohort, race = sl[0], sl[1], sl[2]
        metadata[id] = md_entry_t(cohort, race)

del_vcf_header = allel.read_vcf_headers(cmd_args.del_vcf)
dup_vcf_header = allel.read_vcf_headers(cmd_args.dup_vcf)

field_list = ['variants/ID', 'variants/CHROM', 'variants/POS', 'variants/END', 'variants/SVLEN',
              'variants/INS_TYPE', 'calldata/GT']
del_vcf_file = allel.read_vcf(cmd_args.del_vcf, fields=field_list)
dup_vcf_file = allel.read_vcf(cmd_args.dup_vcf, fields=field_list)

def unique_elems(l):
    return list(set(l))
cohorts = sorted(np.array(unique_elems([metadata[sample].cohort for sample in metadata])))
races = np.array(unique_elems([metadata[sample].race for sample in metadata]))

del_gt_array = allel.GenotypeArray(del_vcf_file['calldata/GT'])
dup_gt_array = allel.GenotypeArray(dup_vcf_file['calldata/GT'])

del_allele_counts, dup_allele_counts = del_gt_array.count_alleles(), dup_gt_array.count_alleles()
ac0_count = 0

del_sv_universe_sizes, del_sv_in_pop_sizes, del_allele_freqs = [], [], []
dup_sv_universe_sizes, dup_sv_in_pop_sizes, dup_allele_freqs = [], [], []
for i, coos in enumerate(zip(del_vcf_file['variants/POS'], del_vcf_file['variants/END'],
                             del_vcf_file['variants/SVLEN'])):
    ac = del_allele_counts[i][1]
    if ac == 0:
        ac0_count += 1
        continue
    an = sum(del_allele_counts[i])
    del_sv_universe_sizes.append(abs(coos[2]))
    del_sv_in_pop_sizes.extend([abs(coos[2]) for x in range(ac)])
    del_allele_freqs.append(ac/an)
for i, coos in enumerate(zip(dup_vcf_file['variants/POS'], dup_vcf_file['variants/END'],
                             dup_vcf_file['variants/SVLEN'])):
    ac = dup_allele_counts[i][1]
    if ac == 0:
        ac0_count += 1
        continue
    an = sum(dup_allele_counts[i])
    dup_sv_universe_sizes.append(abs(coos[2]))
    dup_sv_in_pop_sizes.extend([abs(coos[2]) for x in range(ac)])
    dup_allele_freqs.append(ac/an)

# Sizes histograms
logbins = np.logspace(np.log10(50), np.log10(10000), 100)

del_hist, del_bins = np.histogram(del_sv_in_pop_sizes, bins=logbins)
dup_hist, dup_bins = np.histogram(dup_sv_in_pop_sizes, bins=logbins)
plt.plot(del_bins[1:], del_hist, label="Deletions")
plt.plot(dup_bins[1:], dup_hist, label="Duplications")
plt.xscale('log')
plt.xlabel('Length (bp)')
plt.ylabel('Count')
plt.legend()
plt.savefig(cmd_args.outdir + "/sv_sizes_w_ac_hist.png", bbox_inches='tight')
plt.clf()
plt.close()

# AF plots

n_rare = len([i for i in del_allele_freqs if i < 0.01])
print("RARE DELETIONS:", n_rare)
n_singletons = len([i for i in del_allele_counts if i[1] == 1])
print("SINGLETON DELETIONS:", n_singletons)
print("TOT DELETIONS:", len(del_allele_freqs))

n_rare = len([i for i in dup_allele_freqs if i < 0.01])
print("RARE DUPLICATIONS:", n_rare)
n_singletons = len([i for i in dup_allele_counts if i[1] == 1])
print("SINGLETON DUPLICATIONS:", n_singletons)
print("TOT DUPLICATIONS:", len(dup_allele_freqs))

# Count plots

del_counts_by_sample, dup_counts_by_sample = dict(), dict()
for i, sample in enumerate(del_vcf_header.samples):
    a = del_gt_array.values[:,i]
    del_counts_by_sample[sample] = np.count_nonzero(a == 1)
for i, sample in enumerate(dup_vcf_header.samples):
    a = dup_gt_array.values[:, i]
    dup_counts_by_sample[sample] = np.count_nonzero(a == 1)

del_counts_by_cohort, dup_counts_by_cohort = defaultdict(list), defaultdict(list)
for sample in del_counts_by_sample:
    del_counts_by_cohort[metadata[sample].cohort].append(del_counts_by_sample[sample])
plt.xscale('linear')
plt.boxplot([del_counts_by_cohort[cohort] for cohort in cohorts], labels=cohorts)
plt.ylabel('Count')
plt.savefig(cmd_args.outdir + '/DEL_cohort_ACs_box.png', bbox_inches='tight')
plt.clf()
plt.close()

for sample in dup_counts_by_sample:
    dup_counts_by_cohort[metadata[sample].cohort].append(dup_counts_by_sample[sample])
plt.xscale('linear')
plt.boxplot([dup_counts_by_cohort[cohort] for cohort in cohorts], labels=cohorts)
plt.ylabel('Count')
plt.savefig(cmd_args.outdir + '/DUP_cohort_ACs_box.png', bbox_inches='tight')
plt.clf()
plt.close()

# African vs non-African counts
print("AFR DEL:", statistics.median(del_counts_by_cohort["AFR"]))
nonafr = []
for cohort in del_counts_by_cohort:
    if cohort != "AFR":
        nonafr.extend(del_counts_by_cohort[cohort])
print("NON-AFR DEL:", statistics.median(nonafr))

print("AFR DUP:", statistics.median(dup_counts_by_cohort["AFR"]))
nonafr = []
for cohort in dup_counts_by_cohort:
    if cohort != "AFR":
        nonafr.extend(dup_counts_by_cohort[cohort])
print("NON-AFR DUP:", statistics.median(nonafr))

sample_cohort_labels = np.array([metadata[sample].cohort for sample in del_vcf_header.samples])
sample_race_labels = np.array([metadata[sample].race for sample in del_vcf_header.samples])

# PCA plots

def run_pca(gt_array, allele_counts):
    matrix = gt_array.to_n_alt(fill=-1, dtype=float)
    matrix2 = gt_array.to_n_alt(dtype=float)
    no_gt_yes_cn = (matrix == -1) & (matrix2 == 1)  # entries with 1/.
    matrix[no_gt_yes_cn] = 1
    for i, line in enumerate(matrix):
        gt_samples = (allele_counts[i][0]+allele_counts[i][1])/2
        sv_avg_gt = allele_counts[i][1]/gt_samples if gt_samples > 0 else 0
        line[line == -1] = sv_avg_gt
    matrix = matrix.T
    pca = PCA(n_components=3)
    t = pca.fit_transform(matrix)
    return t, pca

def generate_scatter(reduced_matrix, i1, i2, sample_labels, filename, exp_var):
    fig, ax = plt.subplots()
    for i, g in enumerate(np.unique(sample_labels)):
        ix = np.where(sample_labels == g)
        ax.scatter(reduced_matrix[ix, i1], reduced_matrix[ix, i2], label=g, c=colors[i], s=1)
    ax.set_xlabel('PC%d (%.4f)' % (i1+1, exp_var[i1]))
    ax.set_ylabel('PC%d (%.4f)' % (i2+1, exp_var[i2]))
    if len(np.unique(sample_labels)) <= 10:
        ax.legend()
    plt.savefig(filename, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()

def generate_pca_figs(gt_array, vcf_file, prefix, by_cohort = False):
    if gt_array.n_samples >= 3:
        nosexchr_mask = (vcf_file['variants/CHROM'] != "chrX") & (vcf_file['variants/CHROM'] != "chrY")
        nosexchr_gt_values = gt_array.values[nosexchr_mask]
        nosexchr_gt_array = allel.GenotypeArray(nosexchr_gt_values)
        nosexchr_allele_counts = nosexchr_gt_array.count_alleles()

        pca_reduced_matrix, pca = run_pca(nosexchr_gt_array, nosexchr_allele_counts)
        print(pca.explained_variance_ratio_)

        generate_scatter(pca_reduced_matrix, 0, 1, sample_cohort_labels,
                         cmd_args.outdir + ("/%spca12-superpop.png" % prefix),
                         pca.explained_variance_ratio_)

        if by_cohort:
            for cohort in cohorts:
                cohort_idx = np.where(sample_cohort_labels == cohort)[0]
                cohort_gt_array = nosexchr_gt_array.subset(sel1=cohort_idx)
                cohort_allele_counts = cohort_gt_array.count_alleles()
                pca_reduced_matrix, pca = run_pca(cohort_gt_array, cohort_allele_counts)

                generate_scatter(pca_reduced_matrix, 0, 1, sample_race_labels[cohort_idx],
                                cmd_args.outdir + ("/%spca12-%s.png" % (prefix, cohort)), pca.explained_variance_ratio_)

all_gt_array_values = np.concatenate([del_gt_array.values, dup_gt_array.values], axis=0)
all_gt_array = allel.GenotypeArray(all_gt_array_values)
all_vcf_file = dict()
all_vcf_file['variants/CHROM'] = np.append(del_vcf_file['variants/CHROM'], dup_vcf_file['variants/CHROM'])
all_vcf_file['variants/POS'] = np.append(del_vcf_file['variants/POS'], dup_vcf_file['variants/POS'])
all_vcf_file['variants/END'] = np.append(del_vcf_file['variants/END'], dup_vcf_file['variants/END'])
generate_pca_figs(del_gt_array, del_vcf_file, "DEL_", False)
generate_pca_figs(dup_gt_array, dup_vcf_file, "DUP_", False)
generate_pca_figs(all_gt_array, all_vcf_file, "ALL_", True)
