# SNPGenie_Data_Prep

## Preparing raw data

### annotation files

The folder `annotation_files/` should contain at least `feature table` and `reference sequence`. 

1. `reference sequence`
Names of `reference sequence` files **must** starts with "sequences_", e.g. "sequence_Astroviridae_complete.txt". These reference sequence files should be orgnized in fasta format, while file extension doesn't matter.

2. `feature table`
Names of `feature table` files **must** starts with "Feature table file-", e.g. "sequence_Astroviridae_complete.txt". They should be in NCBI recognizable feature table format.

### VCF files

`variant calls` should be generated 