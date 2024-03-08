import argparse
import pysam

def update_vcf_with_ref_alt(vcf_path, ref_fasta_path, output_vcf_path):
    ref_fasta = pysam.FastaFile(ref_fasta_path)
    vcf_in = pysam.VariantFile(vcf_path, "r")
    vcf_out = pysam.VariantFile(output_vcf_path, "w", header=vcf_in.header)

    for record in vcf_in:
        chrom = record.chrom
        start = record.pos
        end = record.stop

        ref_seq = ref_fasta.fetch(chrom, start - 1, end)
        alt_seq = ref_seq[0]

        record.ref = ref_seq
        record.alts = (alt_seq,)

        vcf_out.write(record)

    ref_fasta.close()
    vcf_in.close()
    vcf_out.close()

parser = argparse.ArgumentParser(description='Update VCF REF and ALT fields based on a reference FASTA file.')
parser.add_argument('-v', '--vcf', required=True, help='Input VCF file path')
parser.add_argument('-r', '--ref', required=True, help='Reference FASTA file path')
parser.add_argument('-o', '--output', required=True, help='Output VCF file path')
args = parser.parse_args()

update_vcf_with_ref_alt(args.vcf, args.ref, args.output)
