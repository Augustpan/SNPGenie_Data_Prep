#! /bin/bash

source /home/dell2/Miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate renv_yfpan

mkdir sw_log_files
rm -rf sw_log_files/*

for line in $(cat manifest_snpgenie)
do
    echo $line
    Rscript /home/dell2/yfpan/SNPGenie/SNPGenie_sliding_windows.R \
        results/$line/codon_results.txt \
        N \
        S \
        40 \
        1 \
        1 \
        100 \
        NONE \
        40 \
        > sw_log_files/$line.log
done
