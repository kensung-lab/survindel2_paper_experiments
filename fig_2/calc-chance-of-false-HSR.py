import sys

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pyfaidx, argparse, pysam, ssw
import random
import numpy as np

cmd_parser = argparse.ArgumentParser(description='Given a VCF file with deletions, a BAM file and a TRF file, '
                                                 'classify deletions within tandem repeats.')
cmd_parser.add_argument('bam_file')
cmd_parser.add_argument('trf_file')
cmd_parser.add_argument('reference')
cmd_parser.add_argument('outdir')
cmd_args = cmd_parser.parse_args()

bam_file = pysam.AlignmentFile(cmd_args.bam_file, reference_filename=cmd_args.reference)
reference_fa = pyfaidx.Fasta(cmd_args.reference)
aligner = ssw.Aligner(matrix=ssw.NucleotideScoreMatrix(match=1, mismatch=-4), gap_open=6, gap_extend=1)

def is_right_clipped(r, min_clip=15):
    return r.cigartuples[-1][0] == 4 and r.cigartuples[-1][1] >= min_clip

def is_left_clipped(r, min_clip=15):
    return r.cigartuples[0][0] == 4 and r.cigartuples[0][1] >= min_clip

def is_aln_left_clipped(aln, min_clip):
    return aln.query_begin >= min_clip

def is_aln_right_clipped(aln, min_clip):
    end_len = len(aln.query) - aln.query_end - 1
    return end_len >= min_clip

def is_aln_clipped(aln, min_clip):
    return is_aln_left_clipped(aln, min_clip) or is_aln_right_clipped(aln, min_clip)

repeats = list()
with open(cmd_args.trf_file) as trf_f:
    for line in trf_f:
        if line.startswith("#"):
            continue
        sl = line.split()
        chr, start, end = sl[1], int(sl[2]), int(sl[3])
        if end-start <= 2000 and "_" not in chr:
            repeats.append((chr, start, end))

random.shuffle(repeats)
hsr_generated_indels_fout = open(cmd_args.outdir + "/hsr_generated_indels.sv", "w")

hsr_scorediffs = list()

for i in range(10000):
    reads = list()
    rep_chr, rep_start, rep_end = repeats[i]
    for r in bam_file.fetch(rep_chr, rep_start, rep_end):
        if r.get_tag('AS') < r.query_length:
            reads.append(r)

    if len(reads) > 1000:
        random.shuffle(reads)
        reads = reads[:1000]

    print("REPEAT", i, len(reads), file=sys.stderr)
    ref = str(reference_fa[rep_chr][rep_start-150:rep_end+150])
    for r in reads:
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        if is_left_clipped(r) or is_right_clipped(r): continue

        full_aln = aligner.align(r.query_sequence, ref, revcomp=False)
        best_score, best_start, best_end = 0, 0, 0
        for j in range(15, r.query_length-15+1):
            lseq, rseq = r.query_sequence[:j], r.query_sequence[j:]
            lseq_aln = aligner.align(lseq, ref, revcomp=False)
            if is_aln_clipped(lseq_aln, 1): continue
            rseq_aln = aligner.align(rseq, ref, revcomp=False)
            if is_aln_clipped(rseq_aln, 1): continue
            split_score = lseq_aln.score + rseq_aln.score
            if split_score > best_score:
                best_score = split_score
                best_start, best_end = lseq_aln.reference_end, rseq_aln.reference_begin
        svlen = best_end-best_start-1
        if best_score > full_aln.score and abs(svlen) >= 50:
            if svlen >= 50:
                svtype = "DEL"
            else:
                svtype = "DUP"
                best_start, best_end = best_end, best_start
            print(r.qname, rep_chr, rep_start+best_start, "N", rep_chr, rep_start+best_end, "N", svtype,
                  best_score-full_aln.score, file=hsr_generated_indels_fout)
            hsr_scorediffs.append(best_score-full_aln.score)

linbins = np.linspace(1, max(hsr_scorediffs), max(hsr_scorediffs))
plt.hist(hsr_scorediffs, bins=linbins, align='left')
plt.ylabel('Hidden split reads')
plt.xlabel('HSR-Score')
plt.savefig(cmd_args.outdir + '/hsr_scorediff_hist.png', bbox_inches='tight')
plt.clf()
plt.close()
