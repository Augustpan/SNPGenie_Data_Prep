#! /bin/bash

# File:         snp_call_varscan.sh
# Date:         2021-08-24 14:00
# Author:       Yuanfei Pan
# Description:  SNP calls from sorted *.bam files with VarScan2.

source /home/dell2/Miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate bio_yfpan

cd /mnt/BACKUP2/BACKUP_TAN_20210325/Ecology

mkdir varscan_vcfs

for line in `cat manifest_vcf`
do
    sp=(${line//@/ })
    idx=${sp[0]}
    ref=${sp[1]}
    echo $idx;
    samtools mpileup -B -d 100000 -f $ref output/$idx.sorted.bam \
    | varscan mpileup2snp \
        --min-var-freq 0.01 \
        --min-avg-qual 20 \
        --min-coverage 10 \
        --p-value 0.01 \
        --output-vcf 1 > varscan_vcfs/$idx.vcf
done;
