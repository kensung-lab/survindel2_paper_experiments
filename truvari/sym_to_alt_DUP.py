import argparse
import pysam

def update_vcf_ref_alt_for_dup_ins(vcf_path, ref_fasta_path, output_vcf_path):
    ref_fasta = pysam.FastaFile(ref_fasta_path)
    vcf_in = pysam.VariantFile(vcf_path, "r")
    vcf_out = pysam.VariantFile(output_vcf_path, "w", header=vcf_in.header)

    for record in vcf_in:
        svtype = record.info['SVTYPE']
        chrom = record.chrom
        start = record.pos
        end = record.stop

        if svtype == 'DUP':
            # For duplications, fetch the base before the duplication and the duplicated sequence
            base_before_dup = ref_fasta.fetch(chrom, start - 1, start)  # Fetch the base before the duplication
            dup_seq = ref_fasta.fetch(chrom, start, end)  # Fetch the duplicated sequence
            ref_seq = base_before_dup
            alt_seq = base_before_dup + dup_seq

        elif svtype == 'INS':
            # For insertions, use the SVINSSEQ field, removing any '-' characters
            ins_seq = record.info['SVINSSEQ'].replace('-', '')
            base_before_ins = ref_fasta.fetch(chrom, start - 1, start)  # Fetch the base before the insertion for REF
            alt_seq = base_before_ins + ins_seq  # Append the cleaned insertion sequence to REF for ALT
            ref_seq = base_before_ins

        record.ref = ref_seq
        record.alts = (alt_seq,)

        vcf_out.write(record)

    ref_fasta.close()
    vcf_in.close()
    vcf_out.close()

def main():
    parser = argparse.ArgumentParser(description='Update REF and ALT fields for DUP and INS in a VCF file based on a reference FASTA file.')
    parser.add_argument('-v', '--vcf', required=True, help='Input VCF file path')
    parser.add_argument('-r', '--ref', required=True, help='Reference FASTA file path')
    parser.add_argument('-o', '--output', required=True, help='Output VCF file path')
    args = parser.parse_args()

    update_vcf_ref_alt_for_dup_ins(args.vcf, args.ref, args.output)

if __name__ == "__main__":
    main()
