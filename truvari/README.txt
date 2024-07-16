pip install truvari

# truvari requires that ALT representation (i.e., REF being filled with the deleted sequence), while SurVIndel2 uses the SYM representation (empty REF and <SVTYPE> in ALT)
# we convert SurVIndel2 to the REF representation for the two samples we compare
python3 sym_to_alt_DEL.py -v ../callers-results/HG00512/survindel2.ml.DEL.vcf.gz -r ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o HG00512.survindel2.ml.DEL.alt.vcf.gz
tabix -p vcf HG00512.survindel2.ml.DEL.alt.vcf.gz 

python3 sym_to_alt_DEL.py -v ../callers-results/NA24385/survindel2.ml.DEL.vcf.gz -r ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o HG002.survindel2.ml.DEL.alt.vcf.gz
tabix -p vcf HG002.survindel2.ml.DEL.alt.vcf.gz 

# Compare HG00512 deletions for Manta and SurVIndel2
tabix -p vcf ../1_data_preparation/by-sample/HG00512.DEL.vcf.gz
rm -rf dragen-DEL-HG00512/ ; truvari bench -b ../1_data_preparation/by-sample/HG00512.DEL.vcf.gz -c ../dragen/HG00512.sv.PASS.DEL.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o dragen-DEL-HG00512/ --pick multi --sizemax 200000000 --pctseq 0
echo DRAGEN recall `grep -m 1 \"recall\" dragen-DEL-HG00512/log.txt | awk '{printf "%.2f",$2}'`
echo DRAGEN precision `grep -m 1 \"precision\" dragen-DEL-HG00512/log.txt | awk '{printf "%.2f",$2}'`
rm -rf manta-DEL-HG00512/ ; truvari bench -b ../1_data_preparation/by-sample/HG00512.DEL.vcf.gz -c ../callers-results/HG00512/manta.DEL.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o manta-DEL-HG00512/ --pick multi --sizemax 200000000 --pctseq 0
echo Manta recall `grep -m 1 \"recall\" manta-DEL-HG00512/log.txt | awk '{printf "%.2f",$2}'`
echo Manta precision `grep -m 1 \"precision\" manta-DEL-HG00512/log.txt | awk '{printf "%.2f",$2}'`
rm -rf survindel2-DEL-HG00512/ ; truvari bench -b ../1_data_preparation/by-sample/HG00512.DEL.vcf.gz -c HG00512.survindel2.ml.DEL.alt.vcf.gz -o survindel2-DEL-HG00512/ --reference ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --pick multi --sizemax 20000000 --pctseq 0
echo SurVIndel2 recall `grep -m 1 \"recall\" survindel2-DEL-HG00512/log.txt | awk '{printf "%.2f",$2}'`
echo SurVIndel2 precision `grep -m 1 \"precision\" survindel2-DEL-HG00512/log.txt | awk '{printf "%.2f",$2}'`

# Compare HG002 deletions for Manta and SurVIndel2
tabix -p vcf  ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz
rm -rf manta-DEL-HG002/ ; truvari bench -b ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz -c ../callers-results/NA24385/manta.DEL.vcf.gz -o manta-DEL-HG002/ --reference ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --pick multi --sizemax 200000000 --pctseq 0
echo Manta recall `grep -m 1 \"recall\" manta-DEL-HG002/log.txt | awk '{printf "%.2f",$2}'`
echo Manta precision `grep -m 1 \"precision\" manta-DEL-HG002/log.txt | awk '{printf "%.2f",$2}'`
rm -rf survindel2-DEL-HG002/ ; truvari bench -b ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz -c HG002.survindel2.ml.DEL.alt.vcf.gz -o survindel2-DEL-HG002/ --reference ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --pick multi --sizemax 200000000 --pctseq 0
echo SurVIndel2 recall `grep -m 1 \"recall\" survindel2-DEL-HG002/log.txt | awk '{printf "%.2f",$2}'`
echo SurVIndel2 precision `grep -m 1 \"precision\" survindel2-DEL-HG002/log.txt | awk '{printf "%.2f",$2}'`

tabix -p vcf ../1_data_preparation/by-sample/HG00512.DUP.vcf.gz
tabix -p vcf ../1_data_preparation/by-sample/HG00512.INS.vcf.gz
rm -rf manta-DUP-recall-HG00512/ ; truvari bench -b ../1_data_preparation/by-sample/HG00512.DUP.vcf.gz -c ../callers-results/HG00512/manta.DUP+INS.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o manta-DUP-recall-HG00512/ --pick multi --pctseq 0 --pctsize 0 --dup-to-ins --sizemax 200000000
rm -rf manta-DUP-precision-HG00512/ ; truvari bench -b ../1_data_preparation/by-sample/HG00512.INS.vcf.gz -c ../callers-results/HG00512/manta.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o manta-DUP-precision-HG00512/ --pick multi --pctseq 0 --pctsize 0 --dup-to-ins --sizemax 200000000
echo Manta recall `grep -m 1 \"recall\": manta-DUP-recall-HG00512/log.txt | awk '{printf "%.2f",$2}'`
echo Manta precision `grep -m 1 \"precision\": manta-DUP-precision-HG00512/log.txt | awk '{printf "%.2f",$2}'`
rm -rf survindel2-DUP-recall-HG00512/ ; truvari bench -b ../1_data_preparation/by-sample/HG00512.DUP.vcf.gz -c ../callers-results/HG00512/survindel2.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o survindel2-DUP-recall-HG00512/ --pick multi --pctseq 0 --pctsize 0 --dup-to-ins --sizemax 200000000
rm -rf survindel2-DUP-precision-HG00512/ ; truvari bench -b ../1_data_preparation/by-sample/HG00512.INS.vcf.gz -c ../callers-results/HG00512/survindel2.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o survindel2-DUP-precision-HG00512/ --pick multi --pctseq 0 --pctsize 0 --dup-to-ins --sizemax 200000000
echo SurVIndel2 recall `grep -m 1 \"recall\": survindel2-DUP-recall-HG00512/log.txt | awk '{printf "%.2f",$2}'`
echo SurVIndel2 precision `grep -m 1 \"precision\": survindel2-DUP-precision-HG00512/log.txt | awk '{printf "%.2f",$2}'`

tabix -p vcf ../1_data_preparation/by-sample/NA24385.DUP.vcf.gz
tabix -p vcf ../1_data_preparation/by-sample/NA24385.INS.vcf.gz
rm -rf manta-DUP-recall-HG002/ ; truvari bench -b ../1_data_preparation/by-sample/NA24385.DUP.vcf.gz -c ../callers-results/NA24385/manta.DUP+INS.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o manta-DUP-recall-HG002/ --pick multi --pctseq 0 --pctsize 0 --dup-to-ins --sizemax 200000000
rm -rf manta-DUP-precision-HG002/ ; truvari bench -b ../1_data_preparation/by-sample/NA24385.INS.vcf.gz -c ../callers-results/NA24385/manta.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o manta-DUP-precision-HG002/ --pick multi --pctseq 0 --pctsize 0 --dup-to-ins --sizemax 200000000
echo Manta recall `grep -m 1 \"recall\": manta-DUP-recall-HG002/log.txt | awk '{printf "%.2f",$2}'`
echo Manta precision `grep -m 1 \"precision\": manta-DUP-precision-HG002/log.txt | awk '{printf "%.2f",$2}'`
rm -rf survindel2-DUP-recall-HG002/ ; truvari bench -b ../1_data_preparation/by-sample/NA24385.DUP.vcf.gz -c ../callers-results/NA24385/survindel2.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o survindel2-DUP-recall-HG002/ --pick multi --pctseq 0 --pctsize 0 --dup-to-ins --sizemax 200000000
rm -rf survindel2-DUP-precision-HG002/ ; truvari bench -b ../1_data_preparation/by-sample/NA24385.INS.vcf.gz -c ../callers-results/NA24385/survindel2.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o survindel2-DUP-precision-HG002/ --pick multi --pctseq 0 --pctsize 0 --dup-to-ins --sizemax 200000000
echo SurVIndel2 recall `grep -m 1 \"recall\": survindel2-DUP-recall-HG002/log.txt | awk '{printf "%.2f",$2}'`
echo SurVIndel2 precision `grep -m 1 \"precision\": survindel2-DUP-precision-HG002/log.txt | awk '{printf "%.2f",$2}'`

# Benchmark all HGSVC2 with all callers using truvari
for f in ../1_data_preparation/by-sample/*.vcf.gz ; do tabix -p vcf -f $f ; done
for f in ../callers-results/???????/*.vcf.gz ; do tabix -f -p vcf $f ; done
for f in ../callers-results/???????/ ; do sample=`basename $f`; python3 sym_to_alt_DEL.py -v $f/survindel2.ml.DEL.vcf.gz -r ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o survindel2-for-truvari/$sample.vcf.gz ; done
for f in survindel2-for-truvari/*.vcf.gz ; do tabix -p vcf $f ; done
rm -rf hgsvc2/
mkdir hgsvc2/

# DELETIONS

bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do truvari bench -b ../1_data_preparation/by-sample/$sample.DEL.vcf.gz -c ../dragen/$sample.sv.PASS.DEL.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o hgsvc2/dragen-DEL-$sample/ --pick multi --sizemax 200000000 --pctseq 0 ; done
bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do truvari bench -b ../1_data_preparation/by-sample/$sample.DEL.vcf.gz -c ../callers-results/$sample/manta.DEL.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o hgsvc2/manta-DEL-$sample/ --pick multi --sizemax 200000000 --pctseq 0 ; done 
bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do truvari bench -b ../1_data_preparation/by-sample/$sample.DEL.vcf.gz -c survindel2-for-truvari/$sample.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o hgsvc2/survindel2-DEL-$sample/ --pick multi --sizemax 200000000 --pctseq 0 ; done

bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do
    echo $sample
    echo RECALL: `grep -m 1 \"recall\": hgsvc2/dragen-DEL-$sample/log.txt | awk '{printf "%.2f",$2}'`
    echo PRECISION: `grep -m 1 \"precision\": hgsvc2/dragen-DEL-$sample/log.txt | awk '{printf "%.2f",$2}'`
    echo
done > hgsvc2/dragen.DEL.results
bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do
    echo $sample
    echo RECALL: `grep -m 1 \"recall\": hgsvc2/manta-DEL-$sample/log.txt | awk '{printf "%.2f",$2}'`
    echo PRECISION: `grep -m 1 \"precision\": hgsvc2/manta-DEL-$sample/log.txt | awk '{printf "%.2f",$2}'`
    echo
done > hgsvc2/manta.DEL.results
bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do
    echo $sample
    echo RECALL: `grep -m 1 \"recall\": hgsvc2/survindel2-DEL-$sample/log.txt | awk '{printf "%.2f",$2}'`
    echo PRECISION: `grep -m 1 \"precision\": hgsvc2/survindel2-DEL-$sample/log.txt | awk '{printf "%.2f",$2}'`
    echo
done > hgsvc2/survindel2.DEL.results

# DUPLICATIONS

bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do truvari bench -b ../1_data_preparation/by-sample/$sample.DUP.vcf.gz -c ../dragen/$sample.sv.PASS.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o hgsvc2/dragen-DUP-recall-$sample/ --pick multi --pctseq 0 --pctsize 0 --dup-to-ins --sizemax 200000000 ; done
bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do truvari bench -b ../1_data_preparation/by-sample/$sample.DUP.vcf.gz -c ../callers-results/$sample/manta.DUP+INS.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o hgsvc2/manta-DUP-recall-$sample/ --pick multi --pctseq 0 --pctsize 0 --dup-to-ins --sizemax 200000000 ; done 
bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do truvari bench -b ../1_data_preparation/by-sample/$sample.DUP.vcf.gz -c ../callers-results/$sample/survindel2.ml.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o hgsvc2/survindel2-DUP-recall-$sample/ --pick multi --pctseq 0 --pctsize 0 --dup-to-ins --sizemax 200000000 ; done
bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do truvari bench -b ../1_data_preparation/by-sample/$sample.INS.vcf.gz -c ../dragen/$sample.sv.PASS.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o hgsvc2/dragen-DUP-precision-$sample/ --pick multi --pctseq 0 --pctsize 0 --dup-to-ins --sizemax 200000000 ; done
bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do truvari bench -b ../1_data_preparation/by-sample/$sample.INS.vcf.gz -c ../callers-results/$sample/manta.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o hgsvc2/manta-DUP-precision-$sample/ --pick multi --pctseq 0 --pctsize 0 --dup-to-ins --sizemax 200000000 ; done 
bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do truvari bench -b ../1_data_preparation/by-sample/$sample.INS.vcf.gz -c ../callers-results/$sample/survindel2.ml.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o hgsvc2/survindel2-DUP-precision-$sample/ --pick multi --pctseq 0 --pctsize 0 --dup-to-ins --sizemax 200000000 ; done

bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do
    echo $sample
    echo RECALL: `grep -m 1 \"recall\": hgsvc2/dragen-DUP-recall-$sample/log.txt | awk '{printf "%.2f",$2}'`
    echo PRECISION: `grep -m 1 \"precision\": hgsvc2/dragen-DUP-precision-$sample/log.txt | awk '{printf "%.2f",$2}'`
    echo
done > hgsvc2/dragen.DUP.results
bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do
    echo $sample
    echo RECALL: `grep -m 1 \"recall\": hgsvc2/manta-DUP-recall-$sample/log.txt | awk '{printf "%.2f",$2}'`
    echo PRECISION: `grep -m 1 \"precision\": hgsvc2/manta-DUP-precision-$sample/log.txt | awk '{printf "%.2f",$2}'`
    echo
done > hgsvc2/manta.DUP.results
bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | grep -v NA24385 | while read sample ; do
    echo $sample
    echo RECALL: `grep -m 1 \"recall\": hgsvc2/survindel2-DUP-recall-$sample/log.txt | awk '{printf "%.2f",$2}'`
    echo PRECISION: `grep -m 1 \"precision\": hgsvc2/survindel2-DUP-precision-$sample/log.txt | awk '{printf "%.2f",$2}'`
    echo
done > hgsvc2/survindel2.DUP.results
