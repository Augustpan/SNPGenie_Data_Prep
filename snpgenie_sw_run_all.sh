#! /bin/bash

source /home/dell2/Miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate renv_yfpan

for line in $(cat manifest_snpgenie)
do
    echo $line
    Rscript /home/dell2/yfpan/SNPGenie/SNPGenie_sliding_windows.R \
        results/$line/codon_results.txt \
        N \
        S \
        40 \
        1 \
        1000 \
        100 \
        NONE \
        40 \
        > log_files/$line.log
done
