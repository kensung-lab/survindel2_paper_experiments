# Download DRAGEN results
bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | while read sample ; do wget https://1000genomes-dragen.s3.amazonaws.com/data/dragen-3.5.7b/hg38_altaware_nohla-cnv-anchored/$sample/$sample.sv.vcf.gz ; done
bcftools query -l ../1_data_preparation/variants_freeze4_sv_insdel_alt.vcf.gz | while read sample ; do wget https://1000genomes-dragen.s3.amazonaws.com/data/dragen-3.5.7b/hg38_altaware_nohla-cnv-anchored/additional_698_related/$sample/$sample.sv.vcf.gz ; done

for f in ???????.sv.vcf.gz ; do bcftools view -f PASS $f -Oz -o ${f%vcf.gz}PASS.vcf.gz ; tabix -p vcf ${f%vcf.gz}PASS.vcf.gz ; done

for f in ???????.sv.PASS.vcf.gz ; do bcftools view -i "SVTYPE=='DEL'" $f -Oz -o ${f%vcf.gz}DEL.vcf.gz ; tabix -p vcf ${f%vcf.gz}DEL.vcf.gz ; done
for f in ???????.sv.PASS.vcf.gz ; do bcftools view -i "SVTYPE=='DUP'" $f -Oz -o ${f%vcf.gz}DUP.vcf.gz ; tabix -p vcf ${f%vcf.gz}DUP.vcf.gz ; done

for f in ???????.sv.PASS.vcf.gz; do sample=`basename $f | cut -d"." -f1`; echo $sample ; 
    ../SurVClusterer/compare-del ../1_data_preparation/by-sample/$sample.DEL.vcf.gz $f -T ../1_data_preparation/simpleRepeat.bed --report
    echo
done > dragen.DEL.results

for f in ???????.sv.PASS.vcf.gz; do sample=`basename $f | cut -d"." -f1`; echo $sample ; 
    ../SurVClusterer/compare-ins ../1_data_preparation/by-sample/$sample.DUP.vcf.gz $f -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep RECALL
    ../SurVClusterer/compare-ins ../1_data_preparation/by-sample/$sample.INS.vcf.gz $sample.sv.PASS.DUP.vcf.gz -T ../1_data_preparation/simpleRepeat.bed -R ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --report | grep PREC
    echo
done > dragen.DUP.results
