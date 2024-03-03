READS1_FA=$1
READS2_FA=$2
HIFI_READS=$3
OUTDIR=$4

find $OUTDIR/reads_fa/ -size 0 -print -delete
mkdir $OUTDIR/consensus_fa/
parallel ./generate_consensus.sh ::: $OUTDIR/reads_fa/*.fa

mkdir $OUTDIR/alt_alleles_fa/
for f in $OUTDIR/consensus_fa/*.fa ; do echo ">"`basename $f | cut -d"." -f1` ; cat $f | bioawk -c fastx '{print $seq}' ; done > $OUTDIR/alt_alleles_fa/alt.fa

bwa index $OUTDIR/alt_alleles_fa/alt.fa
bwa mem -t 64 $OUTDIR/alt_alleles_fa/alt.fa $READS1_FA $READS2_FA | samtools view -b -F 4 > $OUTDIR/alt_alleles_fa/alt.bam
samtools sort -@ 48 $OUTDIR/alt_alleles_fa/alt.bam -o $OUTDIR/alt_alleles_fa/alt.cs.bam
mv $OUTDIR/alt_alleles_fa/alt.cs.bam $OUTDIR/alt_alleles_fa/alt.bam
samtools index $OUTDIR/alt_alleles_fa/alt.bam

samtools view -f 2 -F 3840 $OUTDIR/alt_alleles_fa/alt.bam -b > $OUTDIR/alt_alleles_fa/alt.proper_pairs.bam
samtools view --input-fmt-option filter="[AS]>=qlen-10" $OUTDIR/alt_alleles_fa/alt.proper_pairs.bam | cut -f1 | sort | uniq -c | awk '$1==2 {print $2}' > $OUTDIR/alt_alleles_fa/good_pairs.qnames
samtools view -N $OUTDIR/alt_alleles_fa/good_pairs.qnames $OUTDIR/alt_alleles_fa/alt.proper_pairs.bam -b > $OUTDIR/alt_alleles_fa/alt.good_pairs.bam
samtools index $OUTDIR/alt_alleles_fa/alt.good_pairs.bam

samtools view --input-fmt-option filter="sclen==0" $OUTDIR/alt_alleles_fa/alt.proper_pairs.bam | cut -f1 | sort | uniq -c | awk '$1==2 {print $2}' > $OUTDIR/alt_alleles_fa/unclipped.qnames
samtools view -N $OUTDIR/alt_alleles_fa/unclipped.qnames $OUTDIR/alt_alleles_fa/alt.proper_pairs.bam -b > $OUTDIR/alt_alleles_fa/alt.unclipped.bam
samtools index $OUTDIR/alt_alleles_fa/alt.unclipped.bam

samtools view --input-fmt-option filter="[AS]==qlen" $OUTDIR/alt_alleles_fa/alt.good_pairs.bam -b > $OUTDIR/alt_alleles_fa/alt.good_pairs_perfect_reads.bam
samtools index $OUTDIR/alt_alleles_fa/alt.good_pairs_perfect_reads.bam

samtools view --input-fmt-option filter="[AS]==qlen" $OUTDIR/alt_alleles_fa/alt.good_pairs.bam | cut -f1 | sort | uniq -c | awk '$1==2 {print $2}' > $OUTDIR/alt_alleles_fa/perfect_pairs.qnames
samtools view -N $OUTDIR/alt_alleles_fa/perfect_pairs.qnames $OUTDIR/alt_alleles_fa/alt.good_pairs.bam -b > $OUTDIR/alt_alleles_fa/alt.perfect_pairs.bam
samtools index $OUTDIR/alt_alleles_fa/alt.perfect_pairs.bam

samtools depth -a $OUTDIR/alt_alleles_fa/alt.good_pairs.bam > $OUTDIR/alt_alleles_fa/alt.good_pairs.depth
samtools depth -a $OUTDIR/alt_alleles_fa/alt.good_pairs_perfect_reads.bam > $OUTDIR/alt_alleles_fa/alt.good_pairs_perfect_reads.depth
samtools depth -a $OUTDIR/alt_alleles_fa/alt.perfect_pairs.bam > $OUTDIR/alt_alleles_fa/alt.perfect_pairs.depth

minimap2 --MD -ax map-hifi -t 64 --eqx $OUTDIR/alt_alleles_fa/alt.fa $HIFI_READS | samtools view -F 4 -b > $OUTDIR/alt_alleles_fa/lr.bam
samtools sort -@ 48 $OUTDIR/alt_alleles_fa/lr.bam -o $OUTDIR/alt_alleles_fa/lr.cs.bam
mv $OUTDIR/alt_alleles_fa/lr.cs.bam $OUTDIR/alt_alleles_fa/lr.bam
samtools index $OUTDIR/alt_alleles_fa/lr.bam
