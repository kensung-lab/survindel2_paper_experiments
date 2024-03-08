import argparse, pyfaidx, ssw, pysam

cmd_parser = argparse.ArgumentParser(description='Study reasons for low support of SVs.')
cmd_parser.add_argument('alt_bam', help='BAM file with good pairs mapped to ALT alleles.')
cmd_parser.add_argument('alt_coords', help='File with the coordinates of the SV on the ALT alleles.')
cmd_parser.add_argument('alt', help='Fasta file with the alt alleles.')
cmd_parser.add_argument('ref', help='Fasta file with the reference alleles.')
cmd_parser.add_argument('ref_w_sv', help='Reference fasta file.')
cmd_parser.add_argument('read_len', type=int, help='Read length.')
cmd_args = cmd_parser.parse_args()

alt = pyfaidx.Fasta(cmd_args.alt)
ref = pyfaidx.Fasta(cmd_args.ref)
ref_w_sv = pyfaidx.Fasta(cmd_args.ref_w_sv)

aligner = ssw.Aligner(matrix=ssw.NucleotideScoreMatrix(match=1, mismatch=-4), gap_open=6, gap_extend=1)

bam_in = pysam.AlignmentFile(cmd_args.alt_bam)

# length of the flanking region to test
FLANKING_LEN = 100

def get_diffs(aln):
    return aln.mismatch_count + len([c for c in aln.iter_cigar if c[1] in ('I', 'D')])

def get_left_clip_len(aln):
    return aln.query_begin

def get_right_clip_len(aln):
    return len(aln.query) - aln.query_end - 1

def get_clip_len(aln):
    return max(get_left_clip_len(aln), get_right_clip_len(aln))

alt_coords = dict()
with open(cmd_args.alt_coords) as f:
    for line in f:
        line = line.strip().split()
        alt_coords[line[0]] = (int(line[1]), int(line[2]))

for chr in ref.keys():
    ref_allele = str(ref[chr])
    ref_allele_w_sv = str(ref_w_sv[chr])

    alt_bp1, alt_bp2 = alt_coords[chr]

    supp_bp = 0
    for r in bam_in.fetch(chr, alt_bp1-15, alt_bp1):
        if r.is_secondary or r.is_supplementary: continue
        if r.get_tag("AS") < r.query_length: continue
        if r.reference_start < alt_bp1-15 and r.reference_end > alt_bp1+15:
            supp_bp += 1
        elif r.reference_start < alt_bp2-15 and r.reference_end > alt_bp2+15:
            supp_bp += 1

    alt_reads_pos = set()
    for i in range(alt_bp1-cmd_args.read_len+15, alt_bp1-15):
        alt_reads_pos.add(i)
    for i in range(alt_bp2-cmd_args.read_len+15, alt_bp2-15):
        alt_reads_pos.add(i)
    alt_reads_pos = [i for i in alt_reads_pos if i >= 0 and i < len(alt[chr])-cmd_args.read_len]

    sr_supp, hsr_supp = 0, 0
    for i in alt_reads_pos:
        read = str(alt[chr][i:i+cmd_args.read_len])
        ref_aln = aligner.align(query=read, reference=ref_allele)
        ref_w_sv_aln = aligner.align(query=read, reference=ref_allele_w_sv)
        ref_clip_len, ref_w_sv_clip_len = get_clip_len(ref_aln), get_clip_len(ref_w_sv_aln)
        ref_diffs, ref_w_sv_diffs = get_diffs(ref_aln), get_diffs(ref_w_sv_aln)
        if ref_clip_len >= 15 and ref_w_sv_clip_len == 0:
            sr_supp += 1
        elif ref_diffs > ref_w_sv_diffs:
            hsr_supp += 1

    sr_perc = sr_supp/len(alt_reads_pos)
    hsr_perc = hsr_supp/len(alt_reads_pos)

    exp_sr_reads = int(supp_bp*sr_perc + 0.5)
    exp_hsr_reads = int(supp_bp*hsr_perc + 0.5)

    print(chr, supp_bp, sr_perc, hsr_perc, exp_sr_reads, exp_hsr_reads)
