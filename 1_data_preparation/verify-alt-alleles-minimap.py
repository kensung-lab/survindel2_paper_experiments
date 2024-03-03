import argparse, pyfaidx, ssw, pysam, os, sys
from collections import defaultdict

cmd_parser = argparse.ArgumentParser(description='Finds if the assembled ALT alleles actually contain the SVs.')
cmd_parser.add_argument('vcf', help='VCF containing the SVs.')
cmd_parser.add_argument('hifi_bam', help='Hifi reads aligned to the reference.')
cmd_parser.add_argument('alt_to_ref_aln', help='Alt contigs aligned to the ref contigs using minimap2, in BAM format.')
cmd_parser.add_argument('ref_fa', help='Reference alleles in fasta.')
cmd_parser.add_argument('alt_fa', help='Alternative alleles in fasta.')
cmd_parser.add_argument('workdir', help='Working directory.')
cmd_args = cmd_parser.parse_args()

aligner = ssw.Aligner(matrix=ssw.NucleotideScoreMatrix(match=1, mismatch=-4), gap_open=6, gap_extend=1)

bam_in = pysam.AlignmentFile(cmd_args.hifi_bam)
alt_to_ref_bam = pysam.AlignmentFile(cmd_args.alt_to_ref_aln)

alt_fa = pyfaidx.Fasta(cmd_args.alt_fa)
ref_fa = pyfaidx.Fasta(cmd_args.ref_fa)

def get_left_clip_len(r):
    if r.cigartuples[0][0] == 4:
        return r.cigartuples[0][1]
    else:
        return 0
    
def get_right_clip_len(r):
    if r.cigartuples[-1][0] == 4:
        return r.cigartuples[-1][1]
    else:
        return 0

# given a cigar and an index, return the index of the query sequence
def get_query_index(cigar, index):
    qindex = 0
    for c in cigar[:index]:
        if c[0] in [0, 1, 4]: qindex += c[1]
    return qindex

def get_ref_index(cigar, index, ref_start):
    rindex = ref_start
    for c in cigar[:index]:
        if c[0] in [0, 2]: rindex += c[1]
    return rindex


vcf_in = pysam.VariantFile(cmd_args.vcf)
reads_for_calling, svs_ok, qidxs, ridxs, closest_lens = [], [], [], [], []
ref_bp_coords_fout = open(os.path.join(cmd_args.workdir, "ref-bp-coords.txt"), "w")
alt_bp_coords_fout = open(os.path.join(cmd_args.workdir, "alt-bp-coords.txt"), "w")
for rec in vcf_in.fetch():
    if rec.id not in ref_fa.keys(): # we excluded this (probably because it was very close to another SV)
        continue

    if rec.id not in alt_fa.keys():
        print(rec.id, "no_alt_allele")
        continue

    alt_to_ref_aln = None
    for r in alt_to_ref_bam.fetch(rec.id):
        if r.qname == rec.id and not r.is_secondary and not r.is_supplementary:
            alt_to_ref_aln = r
            break

    if not alt_to_ref_aln:
        print(rec.id, "no_alt_to_ref_aln")
        continue

    alt_svlens, sv_qidxs, sv_ridxs = [], [], []
    lclip_len, rclip_len = get_left_clip_len(alt_to_ref_aln), get_right_clip_len(alt_to_ref_aln)
    if lclip_len >= rclip_len and lclip_len >= 1000 and 0 < alt_to_ref_aln.reference_start < 32000: # 32000 is due to the limit of ssw
        lc_aln = aligner.align(alt_to_ref_aln.query_sequence[:lclip_len], ref_fa[rec.id][:alt_to_ref_aln.reference_start].seq)
        if lc_aln.query_begin > 0:
            print(rec.id, "lclip_poor_aln")
            continue

        r_del_start, r_del_end = lc_aln.reference_end+1, alt_to_ref_aln.reference_start
        del_len = r_del_end - r_del_start
        q_del_start, q_del_end = lc_aln.query_end+1, alt_to_ref_aln.query_alignment_start
        ins_len = q_del_end - q_del_start
        if del_len >= 50 or ins_len >= 50:
            alt_svlens.append(ins_len-del_len)
            sv_qidxs.append((q_del_start, q_del_end))
            sv_ridxs.append((r_del_start, r_del_end))
    elif rclip_len >= lclip_len and rclip_len >= 1000 and 0 < len(ref_fa[rec.id])-alt_to_ref_aln.reference_end < 32000:
        rc_aln = aligner.align(alt_to_ref_aln.query_sequence[-rclip_len:], ref_fa[rec.id][alt_to_ref_aln.reference_end:].seq)
        if rc_aln.query_end < rclip_len - 1:
            print(rec.id, "rclip_poor_aln")
            continue

        r_del_start, r_del_end = alt_to_ref_aln.reference_end, alt_to_ref_aln.reference_end+rc_aln.reference_begin
        del_len = r_del_end - r_del_start
        q_del_start, q_del_end = alt_to_ref_aln.query_alignment_end, alt_to_ref_aln.query_alignment_end+rc_aln.query_begin
        ins_len = q_del_end - q_del_start
        if del_len >= 50 or ins_len >= 50:
            alt_svlens.append(ins_len-del_len)
            sv_qidxs.append((q_del_start, q_del_end))
            sv_ridxs.append((r_del_start, r_del_end))

    for i, c in enumerate(alt_to_ref_aln.cigartuples):
        if c[0] == 1 and c[1] >= 50:
            alt_svlens.append(c[1])
            qindex = get_query_index(alt_to_ref_aln.cigartuples, i)
            sv_qidxs.append((qindex, qindex+c[1]))
            rindex = get_ref_index(alt_to_ref_aln.cigartuples, i, alt_to_ref_aln.reference_start)
            sv_ridxs.append((rindex, rindex))
        elif c[0] == 2 and c[1] >= 50:
            alt_svlens.append(-c[1])
            qindex = get_query_index(alt_to_ref_aln.cigartuples, i)
            sv_qidxs.append((qindex, qindex))
            rindex = get_ref_index(alt_to_ref_aln.cigartuples, i, alt_to_ref_aln.reference_start)
            sv_ridxs.append((rindex, rindex+c[1]))

    if len(alt_svlens) == 0:
        print(rec.id, "no_svs")
        continue
    
    if len(alt_svlens) > 1:
        alt_svlens, sv_qidxs, sv_ridxs = map(list, zip(*sorted(zip(alt_svlens, sv_qidxs, sv_ridxs), key=lambda x: x[1][0])))
        min_dist = min([sv_qidxs[i][0] - sv_qidxs[i-1][1] for i in range(1, len(sv_qidxs))])
        if min_dist <= 100:
            print(rec.id, "multiple_svs", min_dist, alt_svlens, sv_qidxs)
            continue

    svlen = rec.info["SVLEN"]
    if isinstance(svlen, list) or isinstance(svlen, tuple):
        svlen = svlen[0]
    closest_idx = min(range(len(alt_svlens)), key=lambda i: abs(alt_svlens[i]-svlen))
    if abs(alt_svlens[closest_idx] - svlen) > 100:
        print(rec.id, "svlen_mismatch")
        continue

    alt_len = len(alt_fa[rec.id])
    supp_reads = 0
    for r in bam_in.fetch(rec.id):
        if r.reference_start > 0 or r.reference_end < alt_len: continue
        if any([(c[0] in [1, 2] and c[1] >= 10) for c in r.cigartuples]): continue
        if sum([c[1] for c in r.cigartuples if c[0] == 8])*1000 > alt_len: continue
        reads_for_calling.append(r)
        supp_reads += 1

    if supp_reads < 3:
        print(rec.id, "lt3_supp_reads")
        continue

    svs_ok.append(rec)
    ridxs.append(sv_ridxs[closest_idx])
    qidxs.append(sv_qidxs[closest_idx])
    closest_lens.append(alt_svlens[closest_idx])
    print(rec.id, sv_ridxs[closest_idx][0], sv_ridxs[closest_idx][1], file=ref_bp_coords_fout)
    print(rec.id, sv_qidxs[closest_idx][0], sv_qidxs[closest_idx][1], file=alt_bp_coords_fout)
alt_bp_coords_fout.close()
ref_bp_coords_fout.close()

bam_out = pysam.AlignmentFile(cmd_args.workdir + "/reads_for_calling.bam", "wb", template=bam_in)
for r in reads_for_calling:
    bam_out.write(r)
bam_out.close()

os.system("samtools index " + cmd_args.workdir + "/reads_for_calling.bam")

alt_fa_dir = os.path.abspath(os.path.join(cmd_args.alt_fa, os.pardir))
pepper_deepvariant_outdir = os.path.abspath(cmd_args.workdir) + "/pepper_deepvariant_outdir"
os.makedirs(pepper_deepvariant_outdir, exist_ok=True)
pepper_deepvariant_cmd = ("docker run -v %s:/bam_input:ro -v %s:/ref_input:ro -v %s:/outdir/ "
"kishwars/pepper_deepvariant:r0.8 run_pepper_margin_deepvariant call_variant -b /bam_input/reads_for_calling.bam "
"-f /ref_input/alt.fa -o /outdir -t 100 --hifi --phased_output --skip_final_phased_bam > %s/dv.log 2>&1") \
% (os.path.abspath(cmd_args.workdir), alt_fa_dir, pepper_deepvariant_outdir, cmd_args.workdir)
print(pepper_deepvariant_cmd)
os.system(pepper_deepvariant_cmd)

svs = set()
snv_vcf = pysam.VariantFile(pepper_deepvariant_outdir + "/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.phased.vcf.gz")
alleles1, alleles2, unphased, phased = defaultdict(int), defaultdict(int), defaultdict(int), defaultdict(int)
for rec in snv_vcf.fetch():
    gt = rec.samples[0]["GT"]
    if None not in gt and (gt[0] >= 1 or gt[1] >= 1):
        svs.add(rec.chrom)
        if rec.samples[0].phased:
            phased[rec.chrom] += 1
        else:
            unphased[rec.chrom] += 1
        alleles1[rec.chrom] += gt[0]
        alleles2[rec.chrom] += gt[1]

ref_w_sv_bp_coords_fout = open(os.path.join(cmd_args.workdir, "ref_w_sv-bp-coords.txt"), "w")
with open(cmd_args.workdir + "/ref.fa", "w") as ref_out, open(cmd_args.workdir + "/ref_w_sv.fa", "w") as ref_w_vs_out:
    for sv, closest_len, qidx, ridx in zip(svs_ok, closest_lens, qidxs, ridxs):
        if unphased[sv.id] > 1 or (unphased[sv.id] == 1 and phased[sv.id] > 0) or (alleles1[sv.id] >= 1 and alleles2[sv.id] >= 1):
            print(sv.id, "has_snv", alleles1[sv.id], alleles2[sv.id], unphased[sv.id], phased[sv.id])
        else:
            print(sv.id, "ok")
            print(">" + sv.id, file=ref_out)
            print(ref_fa[sv.id], file=ref_out)
            print(">" + sv.id, file=ref_w_vs_out)
            ref_with_sv = ref_fa[sv.id][:ridx[0]].seq + alt_fa[sv.id][qidx[0]:qidx[1]].seq + ref_fa[sv.id][ridx[1]:].seq
            print(ref_with_sv, file=ref_w_vs_out)
            print(sv.id, ridx[0], ridx[0]+qidx[1]-qidx[0], file=ref_w_sv_bp_coords_fout)
ref_w_sv_bp_coords_fout.close()
