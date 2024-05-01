# Download hg19 reference
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz

# Recall and precision of callers on the HG002 GIABv0.6 benchmark

# Deletions
for caller in delly manta smoove survindel survindel2.ml ; do 
    echo $caller 
    ../SurVClusterer/compare-del benchmark/HG002_SVs_Tier1_v0.6.DEL.PASS.vcf.gz callers-results/$caller.DEL.vcf.gz -T simple-repeats-hg19.bed --report | grep RECALL ;
    ../SurVClusterer/compare-del benchmark/HG002_SVs_Tier1_v0.6.DEL.PASS.vcf.gz callers-results/$caller.DEL.t1.vcf.gz -T simple-repeats-hg19.bed --report | grep PREC ;
    echo
done

# Duplications
for caller in delly manta smoove survindel survindel2.ml ; do 
    echo $caller 
    caller_dups=callers-results/$caller.DUP+INS.vcf.gz
    if [ ! -f $caller_dups ]; then caller_dups=callers-results/$caller.DUP.vcf.gz; fi
    ../SurVClusterer/compare-ins benchmark/HG002_SVs_Tier1_v0.6.INS.PASS.DUP.vcf.gz $caller_dups -T simple-repeats-hg19.bed -R hg19.fa --report | grep RECALL ;
    ../SurVClusterer/compare-ins benchmark/HG002_SVs_Tier1_v0.6.INS.PASS.vcf.gz callers-results/$caller.DUP.t1.vcf.gz -T simple-repeats-hg19.bed -R hg19.fa --report | grep PREC ;
    echo
done

# Recall and precision of callers on the HG002 Sniffles2 benchmark

# Deletions

for caller in delly manta smoove survindel survindel2.ml ; do 
    echo $caller 
    ../SurVClusterer/compare-del benchmark/HG002-sniffles2.DEL.vcf.gz callers-results/$caller.DEL.vcf.gz -T simple-repeats-hg19.bed --report | grep -v F1;
    echo
done

# Duplications

for caller in delly manta smoove survindel survindel2.ml ; do 
    echo $caller 
    caller_dups=callers-results/$caller.DUP+INS.vcf.gz
    if [ ! -f $caller_dups ]; then caller_dups=callers-results/$caller.DUP.vcf.gz; fi
    ../SurVClusterer/compare-ins benchmark/HG002-sniffles2.DUP.vcf.gz $caller_dups -T simple-repeats-hg19.bed -R hg19.fa --report | grep RECALL ;
    ../SurVClusterer/compare-ins benchmark/HG002-sniffles2.INS.vcf.gz callers-results/$caller.DUP.vcf.gz -T simple-repeats-hg19.bed -R hg19.fa --report | grep PREC ;
    echo
done

# Benchmark with truvari

# the current version, 4.2.1, claims that survindel2 VCF contains no calls and I could not figure why; 3.5.0 does not have this problem
pip install truvari==3.5.0

# Convert from SYM to ALT SV representation, needed by truvari
python3 ../truvari/sym_to_alt_DEL.py -v callers-results/survindel2.ml.DEL.vcf.gz -r hg19.fa -o survindel2.ml.DEL.alt.vcf.gz
python3 ../truvari/sym_to_alt_DEL.py -v callers-results/survindel2.ml.DEL.t1.vcf.gz -r hg19.fa -o survindel2.ml.DEL.t1.alt.vcf.gz

# Reheader on the GIAB benchmark because truvari complains about lack of contig lens
~/bin/bcftools reheader --fai hg19.fa.fai benchmark/HG002_SVs_Tier1_v0.6.DEL.PASS.vcf.gz | bcftools view -Oz -o HG002_SVs_Tier1_v0.6.DEL.PASS.vcf.gz
~/bin/bcftools reheader --fai hg19.fa.fai benchmark/HG002_SVs_Tier1_v0.6.INS.PASS.DUP.vcf.gz | bcftools view -Oz -o HG002_SVs_Tier1_v0.6.INS.PASS.DUP.vcf.gz
~/bin/bcftools reheader --fai hg19.fa.fai benchmark/HG002_SVs_Tier1_v0.6.INS.PASS.vcf.gz | bcftools view -Oz -o HG002_SVs_Tier1_v0.6.INS.PASS.vcf.gz
tabix -p vcf HG002_SVs_Tier1_v0.6.DEL.PASS.vcf.gz
tabix -p vcf HG002_SVs_Tier1_v0.6.INS.PASS.DUP.vcf.gz
tabix -p vcf HG002_SVs_Tier1_v0.6.INS.PASS.vcf.gz

tabix -p vcf benchmark/HG002-sniffles2.DEL.vcf.gz
tabix -p vcf benchmark/HG002-sniffles2.DUP.vcf.gz
tabix -p vcf benchmark/HG002-sniffles2.INS.vcf.gz 

tabix -p vcf callers-results/manta.DEL.vcf.gz 
tabix -p vcf callers-results/manta.DEL.t1.vcf.gz
tabix -p vcf callers-results/manta.DUP+INS.vcf.gz
tabix -p vcf callers-results/manta.DUP.vcf.gz
tabix -p vcf callers-results/manta.DUP.t1.vcf.gz

tabix -p vcf survindel2.ml.DEL.alt.vcf.gz
tabix -p vcf survindel2.ml.DEL.t1.alt.vcf.gz
tabix -p vcf callers-results/survindel2.ml.DUP.vcf.gz
tabix -p vcf callers-results/survindel2.ml.DUP.t1.vcf.gz

# Sniffles benchmark

rm -rf manta-DEL-HG002-sniffles ; truvari bench -b benchmark/HG002-sniffles2.DEL.vcf.gz -c callers-results/manta.DEL.vcf.gz -o manta-DEL-HG002-sniffles --reference hg19.fa --multimatch --sizemax 200000000
echo Manta recall `grep recall manta-DEL-HG002-sniffles/summary.txt | awk '{printf "%.2f",$2}'`
echo Manta precision `grep precision manta-DEL-HG002-sniffles/summary.txt | awk '{printf "%.2f",$2}'`

rm -rf survindel2-DEL-HG002-sniffles ; truvari bench -b benchmark/HG002-sniffles2.DEL.vcf.gz -c survindel2.ml.DEL.alt.vcf.gz -o survindel2-DEL-HG002-sniffles --reference hg19.fa --multimatch --sizemax 200000000
echo SurVIndel2 recall `grep recall survindel2-DEL-HG002-sniffles/summary.txt | awk '{printf "%.2f",$2}'`
echo SurVIndel2 precision `grep precision survindel2-DEL-HG002-sniffles/summary.txt | awk '{printf "%.2f",$2}'`

rm -rf manta-DUP-recall-HG002-sniffles ; truvari bench -b benchmark/HG002-sniffles2.DUP.vcf.gz -c callers-results/manta.DUP+INS.vcf.gz -f hg19.fa -o manta-DUP-recall-HG002-sniffles/ --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
rm -rf manta-DUP-precision-HG002-sniffles ; truvari bench -b benchmark/HG002-sniffles2.INS.vcf.gz -c callers-results/manta.DUP.vcf.gz -f hg19.fa -o manta-DUP-precision-HG002-sniffles/ --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
echo Manta recall `grep recall manta-DUP-recall-HG002-sniffles/summary.txt | awk '{printf "%.2f",$2}'`
echo Manta precision `grep precision manta-DUP-precision-HG002-sniffles/summary.txt | awk '{printf "%.2f",$2}'`

rm -rf survindel2-DUP-recall-HG002-sniffles ; truvari bench -b benchmark/HG002-sniffles2.DUP.vcf.gz -c callers-results/survindel2.ml.DUP.vcf.gz -f hg19.fa -o survindel2-DUP-recall-HG002-sniffles/ --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
rm -rf survindel2-DUP-precision-HG002-sniffles ; truvari bench -b benchmark/HG002-sniffles2.INS.vcf.gz -c callers-results/survindel2.ml.DUP.vcf.gz -f hg19.fa -o survindel2-DUP-precision-HG002-sniffles/ --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
echo SurVIndel2 recall `grep recall survindel2-DUP-recall-HG002-sniffles/summary.txt | awk '{printf "%.2f",$2}'`
echo SurVIndel2 precision `grep precision survindel2-DUP-precision-HG002-sniffles/summary.txt | awk '{printf "%.2f",$2}'`

# GIAB benchmark

rm -rf manta-DEL-recall-HG002-giab ; truvari bench -b HG002_SVs_Tier1_v0.6.DEL.PASS.vcf.gz -c callers-results/manta.DEL.vcf.gz -o manta-DEL-recall-HG002-giab --reference hg19.fa --multimatch --sizemax 200000000
rm -rf manta-DEL-precision-HG002-giab ; truvari bench -b HG002_SVs_Tier1_v0.6.DEL.PASS.vcf.gz -c callers-results/manta.DEL.t1.vcf.gz -o manta-DEL-precision-HG002-giab --reference hg19.fa --multimatch --sizemax 200000000
echo Manta recall `grep recall manta-DEL-recall-HG002-giab/summary.txt | awk '{printf "%.2f",$2}'`
echo Manta precision `grep precision manta-DEL-precision-HG002-giab/summary.txt | awk '{printf "%.2f",$2}'`

rm -rf survindel2-DEL-recall-HG002-giab ; truvari bench -b HG002_SVs_Tier1_v0.6.DEL.PASS.vcf.gz -c survindel2.ml.DEL.alt.vcf.gz -o survindel2-DEL-recall-HG002-giab --reference hg19.fa --multimatch --sizemax 200000000
rm -rf survindel2-DEL-precision-HG002-giab ; truvari bench -b HG002_SVs_Tier1_v0.6.DEL.PASS.vcf.gz -c survindel2.ml.DEL.t1.alt.vcf.gz -o survindel2-DEL-precision-HG002-giab --reference hg19.fa --multimatch --sizemax 200000000
echo SurVIndel2 recall `grep recall survindel2-DEL-recall-HG002-giab/summary.txt | awk '{printf "%.2f",$2}'`
echo SurVIndel2 precision `grep precision survindel2-DEL-precision-HG002-giab/summary.txt | awk '{printf "%.2f",$2}'`

rm -rf manta-DUP-recall-HG002-giab ; truvari bench -b HG002_SVs_Tier1_v0.6.INS.PASS.DUP.vcf.gz -c callers-results/manta.DUP+INS.vcf.gz -o manta-DUP-recall-HG002-giab --reference hg19.fa --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
rm -rf manta-DUP-precision-HG002-giab ; truvari bench -b HG002_SVs_Tier1_v0.6.INS.PASS.vcf.gz -c callers-results/manta.DUP.t1.vcf.gz -o manta-DUP-precision-HG002-giab --reference hg19.fa --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
echo Manta recall `grep recall manta-DUP-recall-HG002-giab/summary.txt | awk '{printf "%.2f",$2}'`
echo Manta precision `grep precision manta-DUP-precision-HG002-giab/summary.txt | awk '{printf "%.2f",$2}'`

rm -rf survindel2-DUP-recall-HG002-giab ; truvari bench -b HG002_SVs_Tier1_v0.6.INS.PASS.DUP.vcf.gz -c callers-results/survindel2.ml.DUP.vcf.gz -o survindel2-DUP-recall-HG002-giab --reference hg19.fa --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
rm -rf survindel2-DUP-precision-HG002-giab ; truvari bench -b HG002_SVs_Tier1_v0.6.INS.PASS.vcf.gz -c callers-results/survindel2.ml.DUP.t1.vcf.gz -o survindel2-DUP-precision-HG002-giab --reference hg19.fa --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
echo SurVIndel2 recall `grep recall survindel2-DUP-recall-HG002-giab/summary.txt | awk '{printf "%.2f",$2}'`
echo SurVIndel2 precision `grep precision survindel2-DUP-precision-HG002-giab/summary.txt | awk '{printf "%.2f",$2}'`
