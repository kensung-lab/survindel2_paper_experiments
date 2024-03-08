# the current version, 4.2.1, claims that survindel2 VCF contains no calls and I could not figure why; 3.5.0 does not have this problem
pip install truvari==3.5.0

# truvari requires that ALT representation (i.e., REF being filled with the deleted sequence), while SurVIndel2 uses the SYM representation (empty REF and <SVTYPE> in ALT)
# we convert SurVIndel2 to the REF representation for the two samples we compare
python3 sym_to_alt_DEL.py -v ../callers-results/HG00512/survindel2.ml.DEL.vcf.gz -r ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o HG00512.survindel2.ml.DEL.alt.vcf.gz
tabix -p vcf HG00512.survindel2.ml.DEL.alt.vcf.gz 

python3 sym_to_alt_DEL.py -v ../callers-results/NA24385/survindel2.ml.DEL.vcf.gz -r ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o HG002.survindel2.ml.DEL.alt.vcf.gz
tabix -p vcf HG002.survindel2.ml.DEL.alt.vcf.gz 

# Compare HG00512 deletions for Manta and SurVIndel2
tabix -p vcf ../1_data_preparation/by-sample/HG00512.DEL.vcf.gz
rm -rf manta-DEL-HG00512/ ; truvari bench -b ../1_data_preparation/by-sample/HG00512.DEL.vcf.gz -c ../callers-results/HG00512/manta.DEL.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o manta-DEL-HG00512/ --multimatch --sizemax 200000000
echo Manta recall `grep recall manta-DEL-HG00512/summary.txt | awk '{print int($2*100)/100}'`
echo Manta precision `grep prec manta-DEL-HG00512/summary.txt | awk '{print int($2*100)/100}'`
rm -rf survindel2-DEL-HG00512/ ; truvari bench -b ../1_data_preparation/by-sample/HG00512.DEL.vcf.gz -c HG00512.survindel2.ml.DEL.alt.vcf.gz -o survindel2-DEL-HG00512/ --reference ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --multimatch --sizemax 200000000
echo SurVIndel2 recall `grep recall survindel2-DEL-HG00512/summary.txt | awk '{print int($2*100)/100}'`
echo SurVIndel2 precision `grep prec survindel2-DEL-HG00512/summary.txt | awk '{print int($2*100)/100}'`

# Compare HG002 deletions for Manta and SurVIndel2
tabix -p vcf  ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz
rm -rf manta-DEL-HG002/ ; truvari bench -b ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz -c ../callers-results/NA24385/manta.DEL.vcf.gz -o manta-DEL-HG002/ --reference ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --multimatch --sizemax 200000000
echo Manta recall `grep recall manta-DEL-HG002/summary.txt | awk '{print int($2*100)/100}'`
echo Manta precision `grep prec manta-DEL-HG002/summary.txt | awk '{print int($2*100)/100}'`
rm -rf survindel2-DEL-HG002/ ; truvari bench -b ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz -c HG002.survindel2.ml.DEL.alt.vcf.gz -o survindel2-DEL-HG002/ --reference ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa --multimatch --sizemax 200000000
echo SurVIndel2 recall `grep recall survindel2-DEL-HG002/summary.txt | awk '{print int($2*100)/100}'`
echo SurVIndel2 precision `grep prec survindel2-DEL-HG002/summary.txt | awk '{print int($2*100)/100}'`

tabix -p vcf ../1_data_preparation/by-sample/HG00512.DUP.vcf.gz
tabix -p vcf ../1_data_preparation/by-sample/HG00512.INS.vcf.gz
rm -rf manta-DUP-recall-HG00512/ ; truvari bench -b ../1_data_preparation/by-sample/HG00512.DUP.vcf.gz -c ../callers-results/HG00512/manta.DUP+INS.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o manta-DUP-recall-HG00512/ --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
rm -rf manta-DUP-precision-HG00512/ ; truvari bench -b ../1_data_preparation/by-sample/HG00512.INS.vcf.gz -c ../callers-results/HG00512/manta.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o manta-DUP-precision-HG00512/ --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
echo Manta recall `grep recall manta-DUP-recall-HG00512/summary.txt | awk '{print int($2*100)/100}'`
echo Manta precision `grep prec manta-DUP-precision-HG00512/summary.txt | awk '{print int($2*100)/100}'`
rm -rf survindel2-DUP-recall-HG00512/ ; truvari bench -b ../1_data_preparation/by-sample/HG00512.DUP.vcf.gz -c ../callers-results/HG00512/survindel2.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o survindel2-DUP-recall-HG00512/ --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
rm -rf survindel2-DUP-precision-HG00512/ ; truvari bench -b ../1_data_preparation/by-sample/HG00512.INS.vcf.gz -c ../callers-results/HG00512/survindel2.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o survindel2-DUP-precision-HG00512/ --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
echo SurVIndel2 recall `grep recall survindel2-DUP-recall-HG00512/summary.txt | awk '{print int($2*100)/100}'`
echo SurVIndel2 precision `grep prec survindel2-DUP-precision-HG00512/summary.txt | awk '{print int($2*100)/100}'`

tabix -p vcf ../1_data_preparation/by-sample/NA24385.DUP.vcf.gz
tabix -p vcf ../1_data_preparation/by-sample/NA24385.INS.vcf.gz
rm -rf manta-DUP-recall-HG002/ ; truvari bench -b ../1_data_preparation/by-sample/NA24385.DUP.vcf.gz -c ../callers-results/NA24385/manta.DUP+INS.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o manta-DUP-recall-HG002/ --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
rm -rf manta-DUP-precision-HG002/ ; truvari bench -b ../1_data_preparation/by-sample/NA24385.INS.vcf.gz -c ../callers-results/NA24385/manta.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o manta-DUP-precision-HG002/ --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
echo Manta recall `grep recall manta-DUP-recall-HG002/summary.txt | awk '{print int($2*100)/100}'`
echo Manta precision `grep prec manta-DUP-precision-HG002/summary.txt | awk '{print int($2*100)/100}'`
rm -rf survindel2-DUP-recall-HG002/ ; truvari bench -b ../1_data_preparation/by-sample/NA24385.DUP.vcf.gz -c ../callers-results/NA24385/survindel2.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o survindel2-DUP-recall-HG002/ --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
rm -rf survindel2-DUP-precision-HG002/ ; truvari bench -b ../1_data_preparation/by-sample/NA24385.INS.vcf.gz -c ../callers-results/NA24385/survindel2.DUP.vcf.gz -f ../1_data_preparation/GRCh38_full_analysis_set_plus_decoy_hla.fa -o survindel2-DUP-precision-HG002/ --multimatch --pctsim 0 --pctsize 0 --dup-to-ins --sizemax 200000000
echo SurVIndel2 recall `grep recall survindel2-DUP-recall-HG002/summary.txt | awk '{print int($2*100)/100}'`
echo SurVIndel2 precision `grep prec survindel2-DUP-precision-HG002/summary.txt | awk '{print int($2*100)/100}'`
