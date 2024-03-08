# Download HGSVC2 and the hg38 reference used for the study
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/integrated_callset/variants_freeze4_sv_insdel_alt.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

# Unzip Arabidopsis reference
gunzip -c TAIR10_chr_all.fa.gz > TAIR10_chr_all.fa

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
bcftools query -l variants_freeze4_sv_insdel_alt.vcf.gz | while read sample ; do bcftools view variants_freeze4_sv_insdel_alt.vcf.gz -i "SVTYPE!='DEL'" -s $sample --min-ac=1 -Oz -o by-sample/$sample.INS.vcf.gz ; done


# The following will divide the calls into SR-supported, HSR-supported and without support. 
# IMPORTANT
# It is rather labour and computationally intense, depending on the available resources
# For this reason, I have also included the pre-classified VCFs in the subfolder hg002-pre-classified. 
# You can proceed to Sections 2 and without performing these steps

# It requires:
# python, along with packages scikit-learn, ssw, pysam, pyfaidx, numpy
# spoa (https://github.com/rvaser/spoa, version 4.0.9 was used in the study) executable on PATH
# samtools, bwa (0.7.17-r1188 was used in the study) and minimap2 (2.24 was used in the stusy) on PATH

# We will need 4 things:
# 1: raw HiFi reads in fq format (path stored in HIFI_FQ)
# For the study, we used the following accessions
#SRR10382244
#SRR10382245
#SRR10382246
#SRR10382247
#SRR10382248
#SRR10382249
# 2: same reads, mapped to hg38 with minimap2, producing a sorted and indexed (samtools sort and index) BAM file (HIFI_BAM)
# 3: raw Illumina paired reads (we used reads with accession SRR1766442 to SRR1766486). Stored in $ILLUMINA_FQ_1 and $ILLUMINA_FQ_2
# 4: same reads, mapped to hg38 with bwa mem, sorted and indexed (ILLUMINA_BAM)

./1-extract-reads.sh by-sample/NA24385.DEL.vcf.gz $HIFI_BAM GRCh38_full_analysis_set_plus_decoy_hla.fa hg002-DEL-HGSVC
./2-generate-consensus-and-BAM.sh $ILLUMINA_FQ_1 $ILLUMINA_FQ_2 $HIFI_FQ hg002-DEL-HGSVC/
./3-filter-and-find-reads-supporting-SV.sh by-sample/NA24385.DEL.vcf.gz $ILLUMINA_BAM GRCh38_full_analysis_set_plus_decoy_hla.fa hg002-DEL-HGSVC/
./4-classify-SVs.sh by-sample/NA24385.DEL.vcf.gz simpleRepeat.bed hg002-DEL-HGSVC/ 148
# 148 is the length of the Illumina reads

# this step may print an error for a couple of duplications, it is fine
./1-extract-reads.sh by-sample/NA24385.DUP.vcf.gz $HIFI_BAM GRCh38_full_analysis_set_plus_decoy_hla.fa hg002-DUP-HGSVC
./2-generate-consensus-and-BAM.sh $ILLUMINA_FQ_1 $ILLUMINA_FQ_2 $HIFI_FQ hg002-DUP-HGSVC/
./3-filter-and-find-reads-supporting-SV.sh by-sample/NA24385.DUP.vcf.gz $ILLUMINA_BAM GRCh38_full_analysis_set_plus_decoy_hla.fa hg002-DUP-HGSVC/
./4-classify-SVs.sh by-sample/NA24385.DUP.vcf.gz simpleRepeat.bed hg002-DUP-HGSVC/ 148

# The classified DEL (DUP) will be inside hg002-DEL(DUP)-HGSVC/classification/
