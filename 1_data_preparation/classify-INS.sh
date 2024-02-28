vcf=$1
outdir=$2
ref=$3
mkdir $outdir
outvcf=${vcf%vcf.gz}DUP.vcf.gz

bcftools view $vcf -H | awk '$5!~/</ {print ">"$3":"$1":"$2":"length($5); print $5}' > $outdir/ins_seq.fa
trf $outdir/ins_seq.fa 2 7 7 80 10 50 500 -h -ngs > $outdir/ins_seq.fa.trf

python ~/surveyor-ins-paper/benchmark/remove-duplicatons/classify-INS.py $outdir/ins_seq.fa $outdir/ins_seq.fa.trf $ref > $outdir/ins_seq.svtype

bcftools view $vcf -h > $outdir/INS.vcf
bcftools view $vcf -H | grep -w -f <(awk '$2=="INS"' $outdir/ins_seq.svtype | cut -d" " -f1) >> $outdir/INS.vcf

bcftools view $vcf -h > $outdir/DUP.vcf
bcftools view $vcf -H | grep -w -f <(awk '$2=="DUP"' $outdir/ins_seq.svtype | cut -d" " -f1) >> $outdir/DUP.vcf

bcftools view $outdir/DUP.vcf -O z -o $outvcf
