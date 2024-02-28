# Download HGSVC2 and the hg38 reference used for the study

wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/integrated_callset/variants_freeze4_sv_insdel_alt.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

# Find insertions due to tandem duplications

./classify-INS.sh variants_freeze4_sv_insdel_alt.vcf.gz temp GRCh38_full_analysis_set_plus_decoy_hla.fa 
rm -rf temp

# Separate CNVs by sample

mkdir by-sample
bcftools query -l variants_freeze4_sv_insdel_alt.vcf.gz | while read sample ; do bcftools view variants_freeze4_sv_insdel_alt.vcf.gz -i "SVTYPE=='DEL'" -s $sample --min-ac=1 -Oz -o by-sample/$sample.DEL.vcf.gz ; done
bcftools query -l variants_freeze4_sv_insdel_alt.vcf.gz | while read sample ; do bcftools view variants_freeze4_sv_insdel_alt.DUP.vcf.gz -i "SVTYPE!='DEL'" -s $sample --min-ac=1 -Oz -o by-sample/$sample.DUP.vcf.gz ; done
