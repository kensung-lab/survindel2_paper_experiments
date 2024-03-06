mkdir benchmark
for caller in delly manta smoove survindel survindel2.ml survindel2.ml-ara ; do
    for f in ../1_data_preparation/ara-benchmark/*.DEL.vcf.gz ; do sample=`basename $f | cut -d"." -f1`; 
        echo $sample;
        ../SurVClusterer/compare-del $f ../callers-results/$sample/$caller.DEL.vcf.gz -T ../1_data_preparation/TAIR10_chr_all.bed --report
        echo
    done > benchmark/$caller.ara.DEL.txt
done

# benchmark/${caller}.ara.DEL.txt will contain the sensitivity, precision and f1-score for each sample
# survindel2.ml was trained on humans, survindel2.ml-ara was trained on arabidopsis

for caller in delly manta smoove survindel survindel2.ml survindel2.ml-ara ; do
    for f in ../1_data_preparation/ara-benchmark/*.DUP.vcf.gz ; do sample=`basename $f | cut -d"." -f1`; 
        echo $sample;
        called_dups=../callers-results/$sample/$caller.DUP+INS.vcf.gz
        if [ ! -f $called_dups ]; then
            called_dups=../callers-results/$sample/$caller.DUP.vcf.gz
        fi
        ../SurVClusterer/compare-ins $f $called_dups -T ../1_data_preparation/TAIR10_chr_all.bed -R ../1_data_preparation/TAIR10_chr_all.fa --report | grep RECALL
        ../SurVClusterer/compare-ins ../1_data_preparation/ara-benchmark/$sample.INS.vcf.gz ../callers-results/$sample/$caller.DUP.vcf.gz -T ../1_data_preparation/TAIR10_chr_all.bed -R ../1_data_preparation/TAIR10_chr_all.fa --report | grep PREC
        echo
    done > benchmark/$caller.ara.DUP.txt
done

# cow reference: GCF_002263795.2
# mouse reference: GCF_000001635.27
# rice reference: GCF_001433935.1
# Download the references from NCBI, unzip them, run trf (default parameters) to obtain the list of repetitive regions, translate them to BED format

# ../1_data_preparation/other-benchmarks/ contains the benchmarks (obtained with sniffles2 and hifi reads)
# ../callers-results/bos_taurus-LIB1 and 2 contains the calls for Cattle-LIB1 and Cattle-LIB2
# ../callers-results/mus_musculus/ contains the calls for Mouse
# ../callers-results/MH63/ contains the calls for Rice MH63
# ../callers-results/ZH97/ contains the calls for Rice ZH97

# To benchmark deletions for an organism and a caller,
../SurVClusterer/compare-del $BENCHMARK $CALLED -T $ORGANISM_TRF_BED --report

# To benchmark duplications for an organism and a caller,
# BENCHMARK_DUP and BENCHMARK_INS are the .DUP.vcf.gz and the .INS.vcf.gz, respectively
# CALLED_INS and CALLED_DUP are the .DUP+INS.vcf.gz and the .DUP.vcf.gz, respectively 
# (for callers that only have DUP.vcf.gz, that will be CALLED_INS as well)
../SurVClusterer/compare-ins $BENCHMARK_DUP $CALLED_INS -T $ORGANISM_TRF_BED -R $ORGANISM_FASTA --report | grep RECALL
../SurVClusterer/compare-ins $BENCHMARK_INS $CALLED_DUP -T $ORGANISM_TRF_BED -R $ORGANISM_FASTA --report | grep PREC
