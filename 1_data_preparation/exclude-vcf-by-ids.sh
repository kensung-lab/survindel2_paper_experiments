IN_VCF=$1
OUT_PREFIX=$2
ID_FILE=$3

bcftools view -h $IN_VCF > $OUT_PREFIX
bcftools view -H $IN_VCF | grep -w -f $ID_FILE -v >> $OUT_PREFIX
bcftools view $OUT_PREFIX -O z -o $OUT_PREFIX.vcf.gz
rm $OUT_PREFIX
