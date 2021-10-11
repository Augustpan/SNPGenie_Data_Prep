#! /bin/bash

mkdir results
rm -rf results/*

for line in $(cat manifest_snpgenie)
do
    /home/dell2/yfpan/SNPGenie/snpgenie.pl \
        --minfreq=0.01 \
        --snpreport=snpgenie_input/$line.vcf \
        --vcfformat=4 \
        --fastafile=snpgenie_input/$line.fasta \
        --gtffile=snpgenie_input/$line.gtf \
        --outdir results/$line
done
