mkdir benchmark

# Fig 4a

for caller in delly manta smoove survindel survindel2.ml ; do
    for f in ../1_data_preparation/by-sample/???????.DEL.vcf.gz; do sample=`basename $f | cut -d"." -f1`; echo $sample ; ../SurVClusterer/compare-del $f ../callers-results/$sample/$caller.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed --report; echo ; done >> benchmark/$caller.DEL.txt
done

# benchmark/${caller}.DEL.txt will contain the sensitivity, precision and f1-score for each sample

# Fig 4b

for caller in delly manta smoove survindel survindel2.ml ; do
    for f in ../1_data_preparation/by-sample/???????.DUP.vcf.gz; do sample=`basename $f | cut -d"." -f1`; 
        echo $sample
        if [ -f ../callers-results//$sample/$caller.DUP+INS.vcf.gz ]; then
            ../SurVClusterer/compare-ins $f ../callers-results//$sample/$caller.DUP+INS.vcf.gz -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep RECALL
        else
            ../SurVClusterer/compare-ins $f ../callers-results//$sample/$caller.DUP.vcf.gz -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep RECALL
        fi
        ../SurVClusterer/compare-ins ../1_data_preparation/by-sample/$sample.INS.vcf.gz ../callers-results//$sample/$caller.DUP.vcf.gz -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep PREC
        echo
    done  >> benchmark/$caller.DUP.txt
done 

# benchmark/${caller}.DUP.txt will contain the sensitivity and precision for each sample

# Fig 4c

for caller in delly manta smoove survindel survindel2.ml ; do
    echo $caller, SR support
    ../SurVClusterer/compare-del ../1_data_preparation/hg002-pre-classified/DEL/strong_support_sr.vcf.gz ../callers-results/NA24385/$caller.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed --report | grep RECALL
    echo $caller, strong HSR support
    ../SurVClusterer/compare-del ../1_data_preparation/hg002-pre-classified/DEL/strong_support_hsr.vcf.gz ../callers-results/NA24385/$caller.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed --report | grep RECALL
    echo $caller, weak HSR support
    ../SurVClusterer/compare-del ../1_data_preparation/hg002-pre-classified/DEL/weak_support.vcf.gz ../callers-results/NA24385/$caller.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed --report | grep RECALL
    echo $caller, no support
    ../SurVClusterer/compare-del ../1_data_preparation/hg002-pre-classified/DEL/no_support.vcf.gz ../callers-results/NA24385/$caller.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed --report | grep RECALL
    echo
done 2> /dev/null

# Fig 4d

for caller in delly manta smoove survindel survindel2.ml ; do
    called_dups=../callers-results/NA24385/$caller.DUP+INS.vcf.gz
    if [ ! -f $called_dups ]; then
        called_dups=../callers-results/NA24385/$caller.DUP.vcf.gz
    fi
    echo $caller, SR support
    ../SurVClusterer/compare-ins ../1_data_preparation/hg002-pre-classified/DUP/strong_support_sr.vcf.gz $called_dups -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep RECALL
    echo $caller, strong HSR support
    ../SurVClusterer/compare-ins ../1_data_preparation/hg002-pre-classified/DUP/strong_support_hsr.vcf.gz $called_dups -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep RECALL
    echo $caller, weak HSR support
    ../SurVClusterer/compare-ins ../1_data_preparation/hg002-pre-classified/DUP/weak_support.vcf.gz $called_dups -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep RECALL
    echo $caller, no support
    ../SurVClusterer/compare-ins ../1_data_preparation/hg002-pre-classified/DUP/no_support.vcf.gz $called_dups -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep RECALL
    echo
done 2> /dev/null
