# File:         mapping_all.sh
# Date:         2021-08-23 8:00AM
# Author:       Yuanfei Pan
# Description:  Mapping reads onto reference sequences using bowtie2.

#! /bin/bash

source /home/dell2/Miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate bio_yfpan

cd /mnt/BACKUP2/BACKUP_TAN_20210325/Ecology/norRNA_file

for line in `cat ../manifest`
do
    sp=(${line//@/ })
    file=${sp[0]}
    ref=${sp[1]}
    newname=${sp[2]}
    echo $newname;

    bowtie2 --local --threads 48 -x $ref -U $file -S ../output/$newname.sam; 
    samtools view -@ 36 -bSF4 ../output/$newname.sam > ../output/$newname.bam; 
    samtools sort -@ 36 ../output/$newname.bam > ../output/$newname.sorted.bam; 
    samtools index ../output/$newname.sorted.bam; 
    samtools idxstats ../output/$newname.sorted.bam > ../output/$newname.txt; 
    rm ../output/$newname.sam;
done;
