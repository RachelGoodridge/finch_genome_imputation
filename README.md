# Finch Genome Imputation
GLIMPSE pipeline for genomic imputation and variant calling

## Requirements for Use
- install snakemake package through bioconda
- environment.yaml : must be inside environment folder
- fastp.yaml : must be inside environment folder
- glimpse.yaml : must be inside environment folder
- config.yaml : must be inside project folder and contain full paths
- samples.csv : must contain full paths
- average_depth.py
- reference genome
- reference panel
- genetic maps : must be named "chr1.gmap", "chr2.gmap", etc.
- Snakefile

## Instructions for Running
- install snakemake or activate the environment containing snakemake
- change to the directory containing the snakefile
- run with the command : snakemake --sdm conda -d /path/to/project -c [# of cores]
- to clean up temp files : snakemake --delete-temp-output -d /path/to/project