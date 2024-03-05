# Separate HG002 CNVs into inside tandem repeats and outside tandem repeats. Since duplications are represented as insertions (with the insertion site being a single point, i.e.  POS==END,
# unlike deletions which cover the deleted region), we allow a 10 bp tolerance

bedtools intersect -a ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz -b ../1_data_preparation/simpleRepeat.bed -u -f 0.9 -header | bcftools view -Oz -o NA24385.DEL.inTR.vcf.gz
bedtools intersect -a ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz -b ../1_data_preparation/simpleRepeat.bed -v -f 0.9 -header | bcftools view -Oz -o NA24385.DEL.noTR.vcf.gz
bedtools intersect -a ../1_data_preparation/by-sample/NA24385.DUP.vcf.gz -b <(awk '$2>=10 {print $1"\t"$2-10"\t"$3+10}' ../1_data_preparation/simpleRepeat.bed) -u -header | bcftools view -Oz -o NA24385.DUP.inTR.vcf.gz
bedtools intersect -a ../1_data_preparation/by-sample/NA24385.DUP.vcf.gz -b <(awk '$2>=10 {print $1"\t"$2-10"\t"$3+10}' ../1_data_preparation/simpleRepeat.bed) -v -header | bcftools view -Oz -o NA24385.DUP.noTR.vcf.gz

# Number of deletions in tandem repeats: 5940
bcftools view -H NA24385.DEL.inTR.vcf.gz | wc -l

# Number of deletions outside tandem repeats: 3299
bcftools view -H NA24385.DEL.noTR.vcf.gz | wc -l

# Number of duplications in tandem repeats: 8290
bcftools view -H NA24385.DUP.inTR.vcf.gz | wc -l

# Number of duplications in tandem repeats: 640
bcftools view -H NA24385.DUP.noTR.vcf.gz | wc -l

# Generate Fig. 2b, i.e. histogram with number of TR copies deleted/duplicated by HG002 CNVs
python3 calc-n-copies-deldup-in-TR.py NA24385.DEL.inTR.vcf.gz,NA24385.DUP.inTR.vcf.gz ../1_data_preparation/simpleRepeat.txt NA24385.TR-copies-deleted.pdf

# Fig. 3c - deletions

# SR-, HSR- and Un-supported outside TR: 2196, 112 (67 + 45) and 126
bcftools view -H ../1_data_preparation/hg002-pre-classified/DEL/strong_support_sr.noTR.vcf.gz | wc -l
bcftools view -H ../1_data_preparation/hg002-pre-classified/DEL/strong_support_hsr.noTR.vcf.gz | wc -l
bcftools view -H ../1_data_preparation/hg002-pre-classified/DEL/weak_support.noTR.vcf.gz | wc -l
bcftools view -H ../1_data_preparation/hg002-pre-classified/DEL/no_support.noTR.vcf.gz | wc -l

# SR-, HSR- and Un-supported inside TR: 936, 1251 (650 + 601) and 1529
bcftools view -H ../1_data_preparation/hg002-pre-classified/DEL/strong_support_sr.inTR.vcf.gz | wc -l
bcftools view -H ../1_data_preparation/hg002-pre-classified/DEL/strong_support_hsr.inTR.vcf.gz | wc -l
bcftools view -H ../1_data_preparation/hg002-pre-classified/DEL/weak_support.inTR.vcf.gz | wc -l
bcftools view -H ../1_data_preparation/hg002-pre-classified/DEL/no_support.inTR.vcf.gz | wc -l

# Fig. 3d - duplications

# SR-, HSR- and Un-supported outside TR: 413, 46 (25 + 21) and 37
bcftools view -H ../1_data_preparation/hg002-pre-classified/DUP/strong_support_sr.noTR.vcf.gz | wc -l
bcftools view -H ../1_data_preparation/hg002-pre-classified/DUP/strong_support_hsr.noTR.vcf.gz | wc -l
bcftools view -H ../1_data_preparation/hg002-pre-classified/DUP/weak_support.noTR.vcf.gz | wc -l
bcftools view -H ../1_data_preparation/hg002-pre-classified/DUP/no_support.noTR.vcf.gz | wc -l

# SR-, HSR- and Un-supported inside TR: 1249, 1491 (1004 + 487) and 1560
bcftools view -H ../1_data_preparation/hg002-pre-classified/DUP/strong_support_sr.inTR.vcf.gz | wc -l
bcftools view -H ../1_data_preparation/hg002-pre-classified/DUP/strong_support_hsr.inTR.vcf.gz | wc -l
bcftools view -H ../1_data_preparation/hg002-pre-classified/DUP/weak_support.inTR.vcf.gz | wc -l
bcftools view -H ../1_data_preparation/hg002-pre-classified/DUP/no_support.inTR.vcf.gz | wc -l


