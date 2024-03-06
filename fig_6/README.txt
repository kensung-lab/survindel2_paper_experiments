# Download 1000g CNVs
wget https://ftp.ebi.ac.uk/pub/databases/eva/PRJEB71638/1000g-CNV.no_HLA_EBV_decoy.vcf.gz

# Divide DELs and DUPs
bcftools view 1000g-CNV.no_HLA_EBV_decoy.vcf.gz -i "SVTYPE=='DEL'" -Oz -o DEL.vcf.gz
bcftools view 1000g-CNV.no_HLA_EBV_decoy.vcf.gz -i "SVTYPE!='DEL'" -Oz -o DUP.vcf.gz

# Run the analysis script
# Requires scikit-learn, scikit-allel, numpy and matplotlib
python3 analysis_vcf.py DEL.vcf.gz DUP.vcf.gz 1000g.metadata.tsv figs/
# figs/ will contain the panels of Fig. 6, S16 and S17

# note that the number may be VERY slightly different from the ones in the paper
# this is because variants in some of contigs of the reference used had to be removed because not supported by EBI (e.g., HLA)
# The difference is very minor: 247,822 deletions and 232,227 duplications reported in the paper, 247,781 and 232,193 in the EBI version
