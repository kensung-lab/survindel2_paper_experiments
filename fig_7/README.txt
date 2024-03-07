# Fig. 7a

mkdir benchmark
for caller in nygc survindel2.ml ; do
    for f in ../1_data_preparation/by-sample/???????.DEL.vcf.gz; do sample=`basename $f | cut -d"." -f1`
        if [ "$sample" == "NA24385" ]; then continue; fi
        echo $sample
        ../SurVClusterer/compare-del $f ../callers-results/$sample/$caller.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed --report
        echo
    done >> benchmark/$caller.DEL.txt
done
# benchmark/${caller}.DEL.txt will contain the sensitivity, precision for each sample

# Fig. 7b

# Download and partition SurVIndel2 and NYGC 1000g datasets
wget https://ftp.ebi.ac.uk/pub/databases/eva/PRJEB71638/1000g-CNV.no_HLA_EBV_decoy.vcf.gz
bcftools view 1000g-CNV.no_HLA_EBV_decoy.vcf.gz -i "SVTYPE=='DEL'" -Oz -o survindel2.DEL.vcf.gz
bcftools view 1000g-CNV.no_HLA_EBV_decoy.vcf.gz -i "SVTYPE!='DEL'" -Oz -o survindel2.DUP.vcf.gz

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz
bcftools view 1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz -i "SVTYPE=='DEL'" -Oz -o nygc.DEL.vcf.gz
bcftools view 1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz -i "SVTYPE=='DUP'" -Oz -o nygc.DUP.vcf.gz

# 117165 out of 247781 deletions in SurVIndel2 are shared with NYGC. Therefore, 247781-117165 = 130616 are private to SurVIndel2
../SurVClusterer/compare-del survindel2.DEL.vcf.gz nygc.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed --report -a

# 84727 out of 90259 deletions in NYGC are shared with SurVIndel2. Therefore, 90259-84727 = 5532 are private to NYGC
../SurVClusterer/compare-del nygc.DEL.vcf.gz survindel2.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed --report -a

# Fig. 7c

for caller in nygc survindel2.ml ; do
    for f in ../1_data_preparation/by-sample/???????.DUP.vcf.gz; do sample=`basename $f | cut -d"." -f1`
        if [ "$sample" == "NA24385" ]; then continue; fi
        echo $sample
        if [ -f ../callers-results//$sample/$caller.DUP+INS.vcf.gz ]; then
            ../SurVClusterer/compare-ins $f ../callers-results/$sample/$caller.DUP+INS.vcf.gz -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep RECALL
        else
            ../SurVClusterer/compare-ins $f ../callers-results/$sample/$caller.DUP.vcf.gz -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep RECALL
        fi
        ../SurVClusterer/compare-ins ../1_data_preparation/by-sample/$sample.INS.vcf.gz ../callers-results//$sample/$caller.DUP.vcf.gz -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep PREC
        echo
    done  >> benchmark/$caller.DUP.txt
done 
# benchmark/${caller}.DUP.txt will contain the sensitivity, precision for each sample

# Fig 7d

# 35219 out of 232193 duplicatons in SurVIndel2 are shared with NYGC. Therefore, 232193-35219 = 196974 are private to SurVIndel2
../SurVClusterer/compare-ins survindel2.DUP.vcf.gz nygc.DUP.vcf.gz -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report -a

# 14294 out of 28242 duplications in NYGC are shared with SurVIndel2. Therefore, 28242-14294 = 13948 are private to NYGC
../SurVClusterer/compare-ins nygc.DUP.vcf.gz survindel2.DUP.vcf.gz -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report -a
