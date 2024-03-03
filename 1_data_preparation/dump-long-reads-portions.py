import argparse, pysam, ssw, pyfaidx
import numpy as np
from sklearn.cluster import DBSCAN

cmd_parser = argparse.ArgumentParser(description='Given a region in the genome and a BAM file, output the portion'
                                                 'of the reads covering that region.')
cmd_parser.add_argument('bam_file', help='Input bam file.')
cmd_parser.add_argument('ref', help='Reference fasta.')
cmd_parser.add_argument('region', help='Region.')
cmd_parser.add_argument('svlen', type=int)
cmd_parser.add_argument('out', help='Fasta output.')
cmd_args = cmd_parser.parse_args()

chr, coo = cmd_args.region.split(":")
start, end = map(int, coo.split("-"))

ref_read_len, alt_read_len = end - start, end - start + cmd_args.svlen
ref_reads, alt_reads = [], []

aligner = ssw.Aligner() # a gentler aligner for PB reads than the usual (+1, -4, -6, -1) BWA-style

ref_fa = pyfaidx.Fasta(cmd_args.ref)

reads = list()
dels_inss_bps = list()
qnames = list()
with pysam.AlignmentFile(cmd_args.bam_file) as bam_file:
    for r in bam_file.fetch(chr, start, end):
        del_bps, ins_bps = [0], [0]
        qstart, qend = 0, 0

        if r.is_secondary or r.is_supplementary or r.is_unmapped or r.is_duplicate or not r.query_alignment_sequence: continue
        if r.reference_start <= start and r.cigartuples[-1][0] == 4 and r.cigartuples[-1][1] > 4000 and r.reference_end < end:
            clip = r.query_sequence[r.query_alignment_end:]
            if len(clip) > 5000: clip = clip[:5000]
            aln = aligner.align(clip, ref_fa[chr][end-4000:end+1000].seq)
            if aln.reference_end < 4000: continue
            aln_start_on_ref = end - 4000 + aln.reference_begin
            if aln_start_on_ref-r.reference_end >= 50:
                del_bps.append(aln_start_on_ref-r.reference_end)
            if aln.query_begin >= 50 or r.reference_end-aln_start_on_ref >= 50:
                ins_bps.append(aln.query_begin + max(0, r.reference_end-aln_start_on_ref))
            aln = aligner.align(clip, ref_fa[chr][end-4000:end].seq)
            qend = r.query_alignment_end + aln.query_end
        elif r.reference_end >= end and r.cigartuples[0][0] == 4 and r.cigartuples[0][1] > 4000 and r.reference_start > start:
            clip = r.query_sequence[:r.query_alignment_start]
            if len(clip) > 5000: clip = clip[-5000:]
            aln = aligner.align(clip, ref_fa[chr][start-1000:start+4000].seq)
            if aln.reference_begin > 1000: continue
            aln_end_on_ref = start - 1000 + aln.reference_end
            if r.reference_start-aln_end_on_ref >= 50:
                del_bps.append(r.reference_start-aln_end_on_ref)
            if len(clip)-aln.query_end >= 50 or aln_end_on_ref-r.reference_start >= 50:
                ins_bps.append(len(clip)-aln.query_end + max(0, aln_end_on_ref-r.reference_start))
            aln = aligner.align(clip, ref_fa[chr][start:start+4000].seq)
            qstart = r.query_alignment_start - (len(clip)-aln.query_begin)
        elif r.reference_start > start or r.reference_end <= end: continue

        ap = r.get_aligned_pairs()
        # find starting point
        qstart_i, qend_i = -1, -1
        if qstart == 0:
            for i, p in enumerate(ap):
                if p[1] == start:
                    qstart_i = i
                    break
            while qstart_i < len(ap) and ap[qstart_i][0] is None:
                qstart_i += 1
            qstart = ap[qstart_i][0]

        # find ending point
        if qend == 0:
            for i, p in enumerate(ap):
                if p[1] == end:
                    qend_i = i
                    break
            while qend_i >= 0 and ap[qend_i][0] is None:
                qend_i -= 1
            qend = ap[qend_i][0]

        read_len = qend - qstart
        if read_len < 100: continue

        if qstart != -1 and qend != -1:
            curr_del, curr_ins = 0, 0
            for i in range(qstart_i, qend_i+1):
                if ap[i][0] and curr_del > 0:
                    if curr_del >= 10:
                        del_bps.append(-curr_del)
                    curr_del = 0
                if ap[i][1] and curr_ins > 0:
                    if curr_ins >= 10:
                        ins_bps.append(curr_ins)
                    curr_ins = 0
                if not ap[i][0]:
                    curr_del += 1
                if not ap[i][1]:
                    curr_ins += 1

        dels_inss_bps.append((del_bps, ins_bps))
        reads.append((r.qname, r.query_sequence[qstart:qend]))

if not dels_inss_bps: exit(0)

# cluster points in dels_inss_bps
clustering = DBSCAN(eps=10, min_samples=2).fit(np.array([(sum(x[0]), sum(x[1])) for x in dels_inss_bps]))

svlen_sse = list()
for i in set(clustering.labels_):
    if i == -1: continue
    sse = 0
    for j in np.where(clustering.labels_ == i)[0]:
        read_indels = dels_inss_bps[j][0 if cmd_args.svlen < 0 else 1] # 0 = deletions, 1 = insertions
        closest = min(read_indels, key=lambda x: abs(x - cmd_args.svlen))
        sse += (closest - cmd_args.svlen)**2
    svlen_sse.append((i, sse))
if len(svlen_sse) == 0: exit(0)

# find the cluster with the lowest SSE
cluster = min(svlen_sse, key=lambda x: x[1])[0]

with open(cmd_args.out, "w") as out_fasta:
    for i in np.where(clustering.labels_ == cluster)[0]:
        qname, seq = reads[i]
        print(">%s" % qname, file=out_fasta)
        print(seq, file=out_fasta)
