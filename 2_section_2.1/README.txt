# Separate HG002 CNVs into inside tandem repeats and outside tandem repeats. Since duplications are represented as insertions (with the insertion site being a single point, i.e.  POS==END,
# unlike deletions which cover the deleted region), we allow a 10 bp tolerance

bedtools intersect -a ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz -b ../1_data_preparation/simpleRepeat.bed -u -f 0.9 -header | bcftools view -Oz -o NA24385.DEL.inTR.vcf.gz
bedtools intersect -a ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz -b ../1_data_preparation/simpleRepeat.bed -v -f 0.9 -header | bcftools view -Oz -o NA24385.DEL.noTR.vcf.gz
bedtools intersect -a ../1_data_preparation/by-sample/NA24385.DUP.vcf.gz -b <(awk '$2>=10 {print $1"\t"$2-10"\t"$3+10}' ../1_data_preparation/simpleRepeat.bed) -u -header | bcftools view -Oz -o NA24385.DUP.inTR.vcf.gz
bedtools intersect -a ../1_data_preparation/by-sample/NA24385.DUP.vcf.gz -b <(awk '$2>=10 {print $1"\t"$2-10"\t"$3+10}' ../1_data_preparation/simpleRepeat.bed) -v -header | bcftools view -Oz -o NA24385.DUP.noTR.vcf.gz

$ Fig. 2a

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

# Fig. 2c - deletions

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

# Fig. 2d - duplications

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

# Fig. 2e

# HSR-score < 10: 44924
cat ../1_data_preparation/hg002-pre-classified/DEL/supporting_reads.strong_hsr ../1_data_preparation/hg002-pre-classified/DEL/supporting_reads.weak ../1_data_preparation/hg002-pre-classified/DUP/supporting_reads.strong_hsr ../1_data_preparation/hg002-pre-classified/DUP/supporting_reads.weak | awk '$11-$5<10' | wc -l

# HSR-score >=10 and <20: 32184
cat ../1_data_preparation/hg002-pre-classified/DEL/supporting_reads.strong_hsr ../1_data_preparation/hg002-pre-classified/DEL/supporting_reads.weak ../1_data_preparation/hg002-pre-classified/DUP/supporting_reads.strong_hsr ../1_data_preparation/hg002-pre-classified/DUP/supporting_reads.weak | awk '$11-$5>=10 && $11-$5<20' | wc -l

# HSR-score >=20 and <40: 20940
cat ../1_data_preparation/hg002-pre-classified/DEL/supporting_reads.strong_hsr ../1_data_preparation/hg002-pre-classified/DEL/supporting_reads.weak ../1_data_preparation/hg002-pre-classified/DUP/supporting_reads.strong_hsr ../1_data_preparation/hg002-pre-classified/DUP/supporting_reads.weak | awk '$11-$5>=20 && $11-$5<40' | wc -l

# HSR-score >=40: 12667
cat ../1_data_preparation/hg002-pre-classified/DEL/supporting_reads.strong_hsr ../1_data_preparation/hg002-pre-classified/DEL/supporting_reads.weak ../1_data_preparation/hg002-pre-classified/DUP/supporting_reads.strong_hsr ../1_data_preparation/hg002-pre-classified/DUP/supporting_reads.weak | awk '$11-$5>=40' | wc -l

# Fig. 2f

# This script selects 10,000 random repetitive regions, finds HSRs and finds CNVs supported by such HSRs. 
# It will take a while (a few hours) and requires the ILLUMINA_BAM file. I placed the file generated and used in our study in hsr_generated_indels-pre-computed/
# Note that since it uses random sampling, the set of CNVs you will obtain may be different from what was used to generate fig 2f
# If you decide to use the pre-computed file, you can skip the following two commands
mkdir hsr_generated_indels/
python3 calc-chance-of-false-HSR.py $ILLUMINA_BAM ../1_data_preparation/simpleRepeat.txt ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa hsr_generated_indels/

cat hsr_generated_indels-pre-computed/hsr_generated_indels.sv | awk '$9<10' > hsr_generated_indels-pre-computed/hsr_generated_indels.0to10.sv
cat hsr_generated_indels-pre-computed/hsr_generated_indels.sv | awk '$9>=10 && $9<20' > hsr_generated_indels-pre-computed/hsr_generated_indels.10to20.sv
cat hsr_generated_indels-pre-computed/hsr_generated_indels.sv | awk '$9>=20 && $9<40' > hsr_generated_indels-pre-computed/hsr_generated_indels.20to40.sv
cat hsr_generated_indels-pre-computed/hsr_generated_indels.sv | awk '$9>=40' > hsr_generated_indels-pre-computed/hsr_generated_indels.40+.sv

# HSR-score 0-9: 585/8297 deletions TPs and 422/5235 duplications TPs = 1007/13532 TP CNVs
../SurVClusterer/compare-del ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz hsr_generated_indels-pre-computed/hsr_generated_indels.0to10.sv -T ../1_data_preparation/simpleRepeat.bed --report | grep PREC
../SurVClusterer/compare-ins ../1_data_preparation/by-sample/NA24385.INS.vcf.gz hsr_generated_indels-pre-computed/hsr_generated_indels.0to10.sv -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep PREC

# HSR-score 10-19: 299/3247 deletions TPs and 299/2659 duplications TPs = 598/5906 TP CNVs
../SurVClusterer/compare-del ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz hsr_generated_indels-pre-computed/hsr_generated_indels.10to20.sv -T ../1_data_preparation/simpleRepeat.bed --report | grep PREC
../SurVClusterer/compare-ins ../1_data_preparation/by-sample/NA24385.INS.vcf.gz hsr_generated_indels-pre-computed/hsr_generated_indels.10to20.sv -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep PREC

# HSR-score 20-39: 128/802 deletions TPs and 105/935 duplications TPs = 133/1737 TP CNVs
../SurVClusterer/compare-del ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz hsr_generated_indels-pre-computed/hsr_generated_indels.20to40.sv -T ../1_data_preparation/simpleRepeat.bed --report | grep PREC
../SurVClusterer/compare-ins ../1_data_preparation/by-sample/NA24385.INS.vcf.gz hsr_generated_indels-pre-computed/hsr_generated_indels.20to40.sv -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep PREC

# HSR-score 40+: 14/68 deletions TPs and 1/62 duplications TPs = 15/130 TP CNVs
../SurVClusterer/compare-del ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz hsr_generated_indels-pre-computed/hsr_generated_indels.40+.sv -T ../1_data_preparation/simpleRepeat.bed --report | grep PREC
../SurVClusterer/compare-ins ../1_data_preparation/by-sample/NA24385.INS.vcf.gz hsr_generated_indels-pre-computed/hsr_generated_indels.40+.sv -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep PREC
