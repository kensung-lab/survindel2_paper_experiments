VCF=$1
TRF_BED=$2
OUTDIR=$3
READLEN=$4

rm -rf $OUTDIR/classification/
mkdir $OUTDIR/classification/

awk '$2>=10 {print $1"\t"$2-10"\t"$3+10}' $TRF_BED > $OUTDIR/classification/trf-ext.bed

cat $OUTDIR/supporting_reads | awk '($7>=15 && $13==0) && $11>=$18' > $OUTDIR/classification/supporting_reads.strong_sr
cat $OUTDIR/supporting_reads | awk '($6>=3) && $11>=$18' | awk 'NR==FNR{a[$0];next} !($0 in a)' $OUTDIR/classification/supporting_reads.strong_sr - > $OUTDIR/classification/supporting_reads.strong_hsr
cat $OUTDIR/classification/supporting_reads.strong_sr $OUTDIR/classification/supporting_reads.strong_hsr > $OUTDIR/classification/supporting_reads.strong
cat $OUTDIR/supporting_reads | awk '$11>=$18' | awk 'NR==FNR{a[$0];next} !($0 in a)' $OUTDIR/classification/supporting_reads.strong - > $OUTDIR/classification/supporting_reads.weak

cat $OUTDIR/classification/supporting_reads.strong_sr | awk '{print $9}' | sort | uniq -c | awk '$1>=5 {print $2}' > $OUTDIR/classification/strong_support_sr.ids
./filter-vcf-by-ids.sh $VCF $OUTDIR/classification/strong_support_sr $OUTDIR/classification/strong_support_sr.ids 
bedtools intersect -a $OUTDIR/classification/strong_support_sr.vcf.gz -b $OUTDIR/classification/trf-ext.bed -u -f 0.9 -header | bcftools view -Oz -o $OUTDIR/classification/strong_support_sr.inTR.vcf.gz
bedtools intersect -a $OUTDIR/classification/strong_support_sr.vcf.gz -b $OUTDIR/classification/trf-ext.bed -v -f 0.9 -header | bcftools view -Oz -o $OUTDIR/classification/strong_support_sr.noTR.vcf.gz

cat $OUTDIR/classification/supporting_reads.strong_hsr | awk '{print $9}' | sort | uniq -c | awk '$1>=5 {print $2}' | grep -w -v -f $OUTDIR/classification/strong_support_sr.ids > $OUTDIR/classification/strong_support_hsr.ids
./filter-vcf-by-ids.sh $VCF $OUTDIR/classification/strong_support_hsr $OUTDIR/classification/strong_support_hsr.ids 
bedtools intersect -a $OUTDIR/classification/strong_support_hsr.vcf.gz -b $OUTDIR/classification/trf-ext.bed -u -f 0.9 -header | bcftools view -Oz -o $OUTDIR/classification/strong_support_hsr.inTR.vcf.gz
bedtools intersect -a $OUTDIR/classification/strong_support_hsr.vcf.gz -b $OUTDIR/classification/trf-ext.bed -v -f 0.9 -header | bcftools view -Oz -o $OUTDIR/classification/strong_support_hsr.noTR.vcf.gz

cat $OUTDIR/classification/strong_support_sr.ids $OUTDIR/classification/strong_support_hsr.ids > $OUTDIR/classification/strong_support.ids
rm $OUTDIR/classification/strong_support_sr.ids $OUTDIR/classification/strong_support_hsr.ids

cat $OUTDIR/classification/supporting_reads.weak | awk '{print $9}' | sort | uniq -c | awk '$1>=5 {print $2}' | grep -w -v -f $OUTDIR/classification/strong_support.ids > $OUTDIR/classification/weak_support.ids
./filter-vcf-by-ids.sh $VCF $OUTDIR/classification/weak_support $OUTDIR/classification/weak_support.ids
bedtools intersect -a $OUTDIR/classification/weak_support.vcf.gz -b $OUTDIR/classification/trf-ext.bed -u -f 0.9 -header | bcftools view -Oz -o $OUTDIR/classification/weak_support.inTR.vcf.gz
bedtools intersect -a $OUTDIR/classification/weak_support.vcf.gz -b $OUTDIR/classification/trf-ext.bed -v -f 0.9 -header | bcftools view -Oz -o $OUTDIR/classification/weak_support.noTR.vcf.gz

cat $OUTDIR/ref_vs_ref_w_sv/ref.qnames | grep -v -w -f $OUTDIR/classification/strong_support.ids | grep -v -w -f $OUTDIR/classification/weak_support.ids > $OUTDIR/classification/no_support.ids
~/bin/filter-vcf-by-ids.sh $VCF $OUTDIR/classification/no_support $OUTDIR/classification/no_support.ids
bedtools intersect -a $OUTDIR/classification/no_support.vcf.gz -b $OUTDIR/classification/trf-ext.bed -u -f 0.9 -header | bcftools view -Oz -o $OUTDIR/classification/no_support.inTR.vcf.gz
bedtools intersect -a $OUTDIR/classification/no_support.vcf.gz -b $OUTDIR/classification/trf-ext.bed -v -f 0.9 -header | bcftools view -Oz -o $OUTDIR/classification/no_support.noTR.vcf.gz

python3 study-reasons-for-low-support.py $OUTDIR/alt_alleles_fa/alt.good_pairs.bam $OUTDIR/alt_allele_verification/alt-bp-coords.txt $OUTDIR/alt_alleles_fa/alt.fa $OUTDIR/ref_vs_ref_w_sv/ref.fa $OUTDIR/ref_vs_ref_w_sv/ref_w_sv.fa $READLEN > $OUTDIR/classification/expected_support.txt