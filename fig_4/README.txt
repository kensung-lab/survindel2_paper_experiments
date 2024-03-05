mkdir benchmark

for caller in delly manta smoove survindel survindel2.ml ; do
    for f in ../1_data_preparation/by-sample/???????.DEL.vcf.gz; do sample=`basename $f | cut -d"." -f1`; echo $sample ; ../SurVClusterer/compare-del $f ../callers-results/$sample/$caller.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed --report; echo ; done >> benchmark/$caller.DEL.txt
done

# benchmark/${caller}.DEL.txt will contain the sensitivity, precision and f1-score for each sample

for caller in delly manta smoove survindel survindel2.ml ; do
    for f in ../1_data_preparation/by-sample/???????.DUP.vcf.gz; do sample=`basename $f | cut -d"." -f1`; 
        echo $sample
        if [ -f ../callers-results//$sample/$caller.DUP+INS.vcf.gz ]; then
            ../SurVClusterer/compare-ins $f ../callers-results//$sample/$caller.DUP+INS.vcf.gz -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep RECALL
        else
            ../SurVClusterer/compare-ins $f ../callers-results//$sample/$caller.DUP.vcf.gz -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep RECALL
        fi
        ../SurVClusterer/compare-ins ../1_data_preparation/by-sample/$sample.INS.vcf.gz ../callers-results//$sample/$caller.DUP.vcf.gz -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep PREC
    done  >> benchmark/$caller.DUP.txt
done 
# for f in benchmark/???????.DUP.HGSVC.vcf.gz; do sample=`basename $f | cut -d"." -f1`; echo $sample ; ~/SurVClusterer/compare-ins $f results/$sample/manta.DUP+INS.vcf.gz -T ~/references/1kg-hg38/simpleRepeat.bed
 -R ~/references/1kg-hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa 200 --report; echo ; done > results/manta.HGSVC2.DUP_DUP+INS.txt &
# for f in benchmark/???????.DUP.HGSVC.vcf.gz; do sample=`basename $f | cut -d"." -f1`; echo $sample ; ~/SurVClusterer/compare-ins $f results/$sample/smoove.DUP.vcf.gz -T ~/references/1kg-hg38/simpleRepeat.bed -R
 ~/references/1kg-hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa 200 --report; echo ; done > results/smoove.HGSVC2.DUP_DUP.txt &
# for f in benchmark/???????.DUP.HGSVC.vcf.gz; do sample=`basename $f | cut -d"." -f1`; echo $sample ; ~/SurVClusterer/compare-ins $f results/$sample/survindel.DUP.vcf.gz -T ~/references/1kg-hg38/simpleRepeat.bed
 -R ~/references/1kg-hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa 200 --report; echo ; done > results/survindel.HGSVC2.DUP_DUP.txt &
# for f in benchmark/???????.DUP.HGSVC.vcf.gz; do sample=`basename $f | cut -d"." -f1`; echo $sample ; ~/SurVClusterer/compare-ins $f results/$sample/survindel2.DUP.vcf.gz -T ~/references/1kg-hg38/simpleRepeat.be
d -R ~/references/1kg-hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa 200 --report; echo ; done > results/survindel2.HGSVC2.DUP_DUP.txt &
# for f in benchmark/???????.DUP.HGSVC.vcf.gz; do sample=`basename $f | cut -d"." -f1`; echo $sample ; ~/SurVClusterer/compare-ins $f results/$sample/survindel2.ml.DUP.vcf.gz -T ~/references/1kg-hg38/simpleRepeat
.bed -R ~/references/1kg-hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa 200 --report; echo ; done > results/survindel2.ml.HGSVC2.DUP_DUP.txt &
# for f in benchmark/???????.DUP.HGSVC.vcf.gz; do sample=`basename $f | cut -d"." -f1`; echo $sample ; ~/SurVClusterer/compare-ins $f results/$sample/nygc.DUP+INS.vcf.gz -T ~/references/1kg-hg38/simpleRepeat.bed 
-R ~/references/1kg-hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa 200 --report; echo ; done > results/nygc.HGSVC2.DUP_DUP+INS.txt &

