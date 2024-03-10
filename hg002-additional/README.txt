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
