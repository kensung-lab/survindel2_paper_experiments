from __future__ import print_function
import argparse, pyfaidx, ssw
import math

cmd_parser = argparse.ArgumentParser(description='Classify insertions into INS and DUP.')
cmd_parser.add_argument('orig_fa_seq', help='Fasta file with the original inserted sequence.')
cmd_parser.add_argument('trf_file', help='Output of TRF, run with -ngs.')
cmd_parser.add_argument('reference', help='Reference genome in FASTA format.')
cmd_args = cmd_parser.parse_args()

orig_fa = pyfaidx.Fasta(cmd_args.orig_fa_seq)
reference_fa = pyfaidx.Fasta(cmd_args.reference)

def choose_seq(id, lines):
    ins_seq_len = int(id.split(":")[-1])
    best_seq, max_cov = "", 0
    for line in lines:
        sl = line.split()
        start, end = int(sl[0]), int(sl[1])
        seq = sl[13]
        if end - start > max_cov:
            max_cov = end - start
            best_seq = seq
    if max_cov < ins_seq_len*0.9:
        return ""
    return best_seq

contracted_seqs = dict()

lines = []
ins_id = ""
with open(cmd_args.trf_file) as trf_fin:
    for line in trf_fin:
        if line.startswith("@"):
            if ins_id:
                seq = choose_seq(ins_id, lines)
                if seq:
                    contracted_seqs[ins_id[1:]] = seq
                lines = []
            ins_id = line.rstrip()
        else:
            lines.append(line)
    if ins_id:
        seq = choose_seq(ins_id, lines)
        if seq:
            contracted_seqs[ins_id[1:]] = seq

sw = ssw.Aligner(gap_open=6)
sw.matrix.set_match(1)
sw.matrix.set_mismatch(-4)

for k in orig_fa.keys():
    seq = orig_fa[k]
    contracted_seq = str(seq)
    if k in contracted_seqs:
        contracted_seq = contracted_seqs[k]

    if len(contracted_seq) < 50:
        mult = math.ceil(50.0/len(contracted_seq))
        contracted_seq = mult*contracted_seq

    chr, pos = k.split(":")[-3:-1]
    pos = int(pos)

    # It is a DUP if: we can divide the inserted sequence into A and B, and
    # A is a suffix of ref_seq_left
    # B is a prefix of ref_seq_right
    ref_seq_right = "ACGT"*25 + str(reference_fa[chr][pos:pos+len(contracted_seq)+10])
    p_aln = sw.align("ACGT"*25 + contracted_seq, ref_seq_right, False) # try and force prefixes to align
    if p_aln.query_begin > 5 or p_aln.reference_begin > 5:
        print(k.split(":")[0], "INS")
        continue

    ref_seq_left = str(reference_fa[chr][pos-len(contracted_seq)-10:pos]) + "ACGT"*25
    query = contracted_seq + "ACGT"*25
    s_aln = sw.align(query, ref_seq_left, False)  # try and force suffixes to align
    if s_aln.query_end < len(query)-5 or s_aln.reference_end < len(ref_seq_left)-5:
        print(k.split(":")[0], "INS")
        continue

    query_cov = (p_aln.query_end-p_aln.query_begin-100) + (s_aln.query_end-s_aln.query_begin-100)
    query_cov = query_cov / len(contracted_seq)

    if query_cov >= 0.80:
        print(k.split(":")[0], "DUP")
    else:
        print(k.split(":")[0], "INS")
