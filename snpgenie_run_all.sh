#! /bin/bash

mkdir results
for line in $(cat manifest_snpgenie)
do
    /home/dell2/yfpan/SNPGenie/snpgenie.pl \
        --minfreq=0.01 \
        --snpreport=snpgenie_input_files/$line.vcf \
        --vcfformat=4 \
        --fastafile=snpgenie_input_files/$line.fasta \
        --gtffile=snpgenie_input_files/$line.gtf \
        --outdir results/$line
done