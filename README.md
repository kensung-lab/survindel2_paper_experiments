# survindel2_paper_experiments

This repository contains instructions and code to replicate the figures in the SurVIndel2 paper.

Download and compile SurVClusterer in this folder, i.e.
```
git clone https://github.com/Mesh89/SurVClusterer
cd SurVClusterer/
./build_htslib.sh
cmake -DCMAKE_BUILD_TYPE=Release . && make
```

Then, follow the instructions in 1_data_preparation/README.txt, since they are required to generate figures from all sections.

After that, you can follow the README in whichever figure you are interested in replicating.
