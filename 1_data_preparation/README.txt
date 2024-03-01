# Download HGSVC2 and the hg38 reference used for the study
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/integrated_callset/variants_freeze4_sv_insdel_alt.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

# simpleRepeat.txt was downloaded from the UCSC table browser. We convert it to the BED format since some processes need it
cat simpleRepeat.txt | grep -v "#" | cut -f2,3,4,6,17 | sort -k1,1 -k2,2n > simpleRepeat.bed 

# Find insertions due to tandem duplications
# requires TRF (https://tandem.bu.edu/trf/trf.html) in the PATH, python, and python packages pyfaidx and ssw
./classify-INS.sh variants_freeze4_sv_insdel_alt.vcf.gz temp GRCh38_full_analysis_set_plus_decoy_hla.fa 
rm -rf temp

# The result should have 42219 lines
bcftools view variants_freeze4_sv_insdel_alt.DUP.vcf.gz -H | wc -l 
# 42219

# Separate CNVs by sample
mkdir by-sample
bcftools query -l variants_freeze4_sv_insdel_alt.vcf.gz | while read sample ; do bcftools view variants_freeze4_sv_insdel_alt.vcf.gz -i "SVTYPE=='DEL'" -s $sample --min-ac=1 -Oz -o by-sample/$sample.DEL.vcf.gz ; done
bcftools query -l variants_freeze4_sv_insdel_alt.vcf.gz | while read sample ; do bcftools view variants_freeze4_sv_insdel_alt.DUP.vcf.gz -i "SVTYPE!='DEL'" -s $sample --min-ac=1 -Oz -o by-sample/$sample.DUP.vcf.gz ; done
