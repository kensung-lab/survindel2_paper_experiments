VCF=$1
BAM=$2
REFERENCE=$3
OUTDIR=$4

mkdir -p $OUTDIR/reads_fa/
rm $OUTDIR/ref_alleles.fa
bcftools query -f "%ID %CHROM %POS %END %SVLEN\n" $VCF | awk '$3>2000 {$3-=2000; $4+=2000; print}' | python3 select-nonoverlapping-intervals.py | awk '{print $1,$2":"$3"-"$4,$5}' | while read l; do
	id=`echo $l | cut -d" " -f1`
	reg=`echo $l | cut -d" " -f2`
	svlen=`echo $l | cut -d" " -f3`
	python3 dump-long-reads-portions.py $BAM $REFERENCE $reg $svlen $OUTDIR/reads_fa/$id.fa &
	echo ">"$id >> $OUTDIR/ref_alleles.fa
	samtools faidx $REFERENCE $reg | grep -v ">" >> $OUTDIR/ref_alleles.fa
done
wait
