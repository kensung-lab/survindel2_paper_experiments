# HG00512

# Find matches between HGSVC2 and each caller's calls
../SurVClusterer/compare-del ../1_data_preparation/by-sample/HG00512.DEL.vcf.gz HG00512/delly.DEL.raw.vcf.gz -T ../1_data_preparation/simpleRepeat.bed | grep -v NONE | cut -d" " -f1,2 > HG00512/delly.comp 
../SurVClusterer/compare-del ../1_data_preparation/by-sample/HG00512.DEL.vcf.gz HG00512/manta.DEL.raw.vcf.gz -T ../1_data_preparation/simpleRepeat.bed | grep -v NONE | cut -d" " -f1,2 > HG00512/manta.comp 
../SurVClusterer/compare-del ../1_data_preparation/by-sample/HG00512.DEL.vcf.gz HG00512/smoove.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed | grep -v NONE | cut -d" " -f1,2 > HG00512/smoove.comp
../SurVClusterer/compare-del ../1_data_preparation/by-sample/HG00512.DEL.vcf.gz HG00512/survindel.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed | grep -v NONE | cut -d" " -f1,2 > HG00512/survindel.comp
../SurVClusterer/compare-del ../1_data_preparation/by-sample/HG00512.DEL.vcf.gz HG00512/survindel2.ml.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed | grep -v NONE | cut -d" " -f1,2 > HG00512/survindel2.comp

python3 roc.py --caller_vcf HG00512/delly.DEL.raw.vcf.gz HG00512/manta.DEL.raw.vcf.gz HG00512/smoove.DEL.vcf.gz HG00512/survindel.DEL.vcf.gz HG00512/survindel2.ml.DEL.vcf.gz --benchmark_vcf ../1_data_preparation/by-sample/HG00512.DEL.vcf.gz --equivalents HG00512/delly.comp HG00512/manta.comp HG00512/smoove.comp HG00512/survindel.comp HG00512/survindel2.comp --labels Delly Manta Smoove SurVIndel SurVIndel2 --output HG00512_roc_curve.png 


# HG002

# Find matches between HGSVC2 and each caller's calls
../SurVClusterer/compare-del ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz HG002/delly.DEL.raw.vcf.gz -T ../1_data_preparation/simpleRepeat.bed | grep -v NONE | cut -d" " -f1,2 > HG002/delly.comp 
../SurVClusterer/compare-del ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz HG002/manta.DEL.raw.vcf.gz -T ../1_data_preparation/simpleRepeat.bed | grep -v NONE | cut -d" " -f1,2 > HG002/manta.comp 
../SurVClusterer/compare-del ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz HG002/smoove.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed | grep -v NONE | cut -d" " -f1,2 > HG002/smoove.comp
../SurVClusterer/compare-del ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz HG002/survindel.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed | grep -v NONE | cut -d" " -f1,2 > HG002/survindel.comp
../SurVClusterer/compare-del ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz HG002/survindel2.ml.DEL.vcf.gz -T ../1_data_preparation/simpleRepeat.bed | grep -v NONE | cut -d" " -f1,2 > HG002/survindel2.comp

python3 roc.py --caller_vcf HG002/delly.DEL.raw.vcf.gz HG002/manta.DEL.raw.vcf.gz HG002/smoove.DEL.vcf.gz HG002/survindel.DEL.vcf.gz HG002/survindel2.ml.DEL.vcf.gz --benchmark_vcf ../1_data_preparation/by-sample/NA24385.DEL.vcf.gz --equivalents HG002/delly.comp HG002/manta.comp HG002/smoove.comp HG002/survindel.comp HG002/survindel2.comp --labels Delly Manta Smoove SurVIndel SurVIndel2 --output HG002_roc_curve.png
