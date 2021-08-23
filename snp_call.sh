# File:         snp_call.sh
# Date:         2021-08-23 8:00AM
# Author:       Yuanfei Pan
# Description:  SNP calls from sorted *.bam files.

#! /bin/bash

source /home/dell2/Miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate bio_yfpan

cd /mnt/BACKUP2/BACKUP_TAN_20210325/Ecology

mkdir new_vcfs

for line in `cat manifest_vcf`
do
    sp=(${line//@/ })
    idx=${sp[0]}
    ref=${sp[1]}
    echo $idx;
    bcftools mpileup -d 100000 -f $ref -a FORMAT/AD,FORMAT/SP,INFO/AD -Ou output/$idx.sorted.bam | bcftools call --ploidy 1 -vmO z -o new_vcfs/$idx.vcf
done;
