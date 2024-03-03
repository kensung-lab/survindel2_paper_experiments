VCF=$1
SR_BAM_FILE=$2
REFERENCE=$3
OUTDIR=$4

mkdir $OUTDIR/alt_allele_verification/
minimap2 -ax map-hifi -t 100 $OUTDIR/ref_alleles.fa $OUTDIR/alt_alleles_fa/alt.fa | samtools view -b > $OUTDIR/alt_allele_verification/mm2_alt_to_ref.bam
samtools sort -@ 100 $OUTDIR/alt_allele_verification/mm2_alt_to_ref.bam -o $OUTDIR/alt_allele_verification/mm2_alt_to_ref.cs.bam
samtools index $OUTDIR/alt_allele_verification/mm2_alt_to_ref.cs.bam
python3 verify-alt-alleles-minimap.py $VCF $OUTDIR/alt_alleles_fa/lr.bam $OUTDIR/alt_allele_verification/mm2_alt_to_ref.cs.bam $OUTDIR/ref_alleles.fa $OUTDIR/alt_alleles_fa/alt.fa $OUTDIR/alt_allele_verification/ > $OUTDIR/alt_allele_verification/alt_filter.txt

grep ok $OUTDIR/alt_allele_verification/alt_filter.txt | cut -d" " -f1 > $OUTDIR/alt_allele_verification/alt_filter.filters_ok.txt
rm -rf $OUTDIR/ref_vs_ref_w_sv 
mkdir $OUTDIR/ref_vs_ref_w_sv
mv $OUTDIR/alt_allele_verification/ref*.fa $OUTDIR/ref_vs_ref_w_sv/
grep ">" $OUTDIR/ref_vs_ref_w_sv/ref.fa | cut -c 2- > $OUTDIR/ref_vs_ref_w_sv/ref.qnames

samtools view -H $OUTDIR/alt_alleles_fa/alt.good_pairs.bam > $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.sam
cat $OUTDIR/ref_vs_ref_w_sv/ref.qnames | while read line ; do
    samtools view $OUTDIR/alt_alleles_fa/alt.good_pairs.bam $line >> $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.sam
done
samtools view -b $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.sam > $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.bam
rm $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.sam
samtools sort -n -@ 48 $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.bam -o $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.qs.bam
mv $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.qs.bam $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.bam
samtools fastq -1 $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.1.fq -2 $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.2.fq $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.bam

bwa index $OUTDIR/ref_vs_ref_w_sv/ref.fa
bwa mem -t 64 $OUTDIR/ref_vs_ref_w_sv/ref.fa $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.1.fq $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.2.fq | samtools view -b > $OUTDIR/ref_vs_ref_w_sv/ref.bam
samtools sort -@ 48 $OUTDIR/ref_vs_ref_w_sv/ref.bam -o $OUTDIR/ref_vs_ref_w_sv/ref.cs.bam
mv $OUTDIR/ref_vs_ref_w_sv/ref.cs.bam $OUTDIR/ref_vs_ref_w_sv/ref.bam
samtools index $OUTDIR/ref_vs_ref_w_sv/ref.bam
bwa index $OUTDIR/ref_vs_ref_w_sv/ref_w_sv.fa
bwa mem -t 64 $OUTDIR/ref_vs_ref_w_sv/ref_w_sv.fa $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.1.fq $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.2.fq | samtools view -b > $OUTDIR/ref_vs_ref_w_sv/ref_w_sv.bam
samtools sort -@ 48 $OUTDIR/ref_vs_ref_w_sv/ref_w_sv.bam -o $OUTDIR/ref_vs_ref_w_sv/ref_w_sv.cs.bam
mv $OUTDIR/ref_vs_ref_w_sv/ref_w_sv.cs.bam $OUTDIR/ref_vs_ref_w_sv/ref_w_sv.bam
samtools index $OUTDIR/ref_vs_ref_w_sv/ref_w_sv.bam

python3 find_reads_better_in_bam1.py $OUTDIR/ref_vs_ref_w_sv/ref_w_sv.bam $OUTDIR/ref_vs_ref_w_sv/ref.bam $OUTDIR/alt_alleles_fa/alt.good_pairs.ok_alleles.bam $OUTDIR/alt_allele_verification/alt-bp-coords.txt $SR_BAM_FILE $OUTDIR/supporting_reads.bam > $OUTDIR/supporting_reads
