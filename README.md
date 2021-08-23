# SNPGenie_Data_Prep

## Step 1 - Mapping reads onto genomes

use `mapping_all.sh`

## Step 2 - Variant calling

use `snp_call.sh`

## Step 3 - Gathering data for SNPGenie

### Annotation files

The folder `annotation_files/` should contain at least `feature table` and `reference sequence`. 

1. `reference sequence`
Names of `reference sequence` files **must** starts with "sequences_", e.g. "sequence_Astroviridae_complete.txt". These reference sequence files should be orgnized in fasta format, while file extension doesn't matter.

2. `feature table`
Names of `feature table` files **must** starts with "Feature table file-", e.g. "sequence_Astroviridae_complete.txt". They should be in NCBI recognizable feature table format.

### VCF files

Variant calling files should be generated using `snp_call.sh` and the results should be placed into `vcf_files/`.

### Run `data_prep.py`

If all data are properly gathered and orgnized as described before, a simple run should get everything ready for runing SNPGenie pipeline.

## Step 4 - Runing SNPGenie

use `snpgenie_run_all.sh`