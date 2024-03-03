import argparse, pysam

cmd_parser = argparse.ArgumentParser(description='Given two BAM files bam_file1 and bam_file2, and another BAM file truth_bam, ' 
                                     'with the alignments of the same reads, output the reads that are better aligned in bam_file1 '
                                     'than in bam_file2, and the contig they are aligned to in bam_file1 and truth_bam are the same.')
cmd_parser.add_argument('bam_file1', help='Input bam file 1.')
cmd_parser.add_argument('bam_file2', help='Input bam file 2.')
cmd_parser.add_argument('truth_bam', help='Bam file where the reads are aligned to the correct contigs.')
cmd_parser.add_argument('sv_breakpoints', help='Exclude reads in that do not overlap one of the breakpoints in truth_bam.')
cmd_parser.add_argument('full_bam', help='Original short reads bam file.')
cmd_parser.add_argument('supporting_pairs_bam', help='File where to write supporting pairs.')
cmd_args = cmd_parser.parse_args()

def get_diffs(r):
    if r.is_unmapped: return 0
    return r.get_tag("NM") - sum([x[1] for x in r.cigartuples if x[0] in (1,2)]) + len(([x[1] for x in r.cigartuples if x[0] in (1,2)]))

def get_longest_clip_len(r):
    if r.is_unmapped: return 0
    left_clip = r.cigartuples[0][1] if r.cigartuples[0][0] == 4 else 0
    right_clip = r.cigartuples[-1][1] if r.cigartuples[-1][0] == 4 else 0
    return max(left_clip, right_clip)

def get_qname(r):
    if r.is_read1:
        return r.qname + "_1"
    else:
        return r.qname + "_2"

sv_bkps = dict()
with open(cmd_args.sv_breakpoints) as f:
    for line in f:
        chr, start, end = line.strip().split()
        sv_bkps[chr] = (int(start), int(end))

with pysam.AlignmentFile(cmd_args.truth_bam) as truth_bam:
    truth_rname = dict()
    for r in truth_bam.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary: continue
        if r.get_tag("AS") < r.query_length: continue
        sv_start, sv_end = sv_bkps[r.reference_name]
        if r.reference_start < sv_start < r.reference_end or r.reference_start < sv_end < r.reference_end:
            truth_rname[get_qname(r)] = r.reference_name

bam2_score, bam2_cigar, bam2_rname, bam2_diffs, bam2_cliplen, bam2_unmapped = dict(), dict(), dict(), dict(), dict(), dict()
with pysam.AlignmentFile(cmd_args.bam_file2) as bam_file2:
    for r in bam_file2.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary: continue
        qname = get_qname(r)
        bam2_score[qname] = r.get_tag("AS")
        bam2_cigar[qname] = r.cigarstring
        bam2_rname[qname] = r.reference_name
        bam2_diffs[qname] = get_diffs(r)
        bam2_cliplen[qname] = get_longest_clip_len(r)
        bam2_unmapped[qname] = r.is_unmapped

better_qnames = set()
better_reads = dict()
with pysam.AlignmentFile(cmd_args.bam_file1) as bam_file1:
    for r in bam_file1.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary: continue
        qname = get_qname(r)
        
        if qname not in truth_rname or truth_rname[qname] != r.reference_name:
            continue

        if r.get_tag("AS") > bam2_score[qname]:
            better_reads[qname] = r
            better_qnames.add(r.qname)

with pysam.AlignmentFile(cmd_args.full_bam) as full_bam, pysam.AlignmentFile(cmd_args.supporting_pairs_bam, "wb", template=full_bam) as supporting_pairs_bam:
    for r in full_bam.fetch(until_eof=True):
        if r.qname in better_qnames:
            supporting_pairs_bam.write(r)
        if r.is_secondary or r.is_supplementary: continue
        qname = get_qname(r)
        if qname in better_reads:
            br = better_reads[qname]
            score = br.get_tag("AS")
            print(br.qname, br.query_sequence[:10], end="\t")
            print(bam2_rname[qname], bam2_cigar[qname], bam2_score[qname], bam2_diffs[qname], bam2_cliplen[qname], bam2_unmapped[qname], end="\t")
            print(br.reference_name, br.cigarstring, score, get_diffs(br), get_longest_clip_len(br), br.is_unmapped, end="\t")
            print(r.reference_name, r.pos, r.cigarstring, r.get_tag("AS"))
