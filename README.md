# Finch Genome Imputation
GLIMPSE pipeline for genomic imputation and variant calling

## Steps for Installation
1. install snakemake : conda create -c conda-forge -c bioconda -n snakemake "snakemake=9.1.10" "python==3.11.4"
2. activate the snakemake environment : conda activate snakemake
3. install plugins (in the same environment) : conda install -c bioconda -c conda-forge snakemake-executor-plugin-slurm snakemake-storage-plugin-fs
4. install biopython (in the same environment) : conda install -c conda-forge biopython=1.85

## Requirements for Use
- environments in the "envs" folder
    - environment.yaml
    - fasterq.yaml
    - fastp.yaml
    - glimpse.yaml
- profiles/config.yaml : must be inside a folder called "profiles" and contain the slurm account, if using slurm
- project/config.yaml : must be inside project folder and contain full paths
- samples.csv : must contain a unqiue identifier for each sample (unless sample has more than one run) and full path to the data, data may be stored remotely
- average_depth.py
- reference genome : a .fna file
- reference panel : a .vcf or .bcf file
- genetic maps : must be named "chr1.gmap", "chr2.gmap", etc. and contained in a single folder, names must exactly match as in reference genome and reference panel
- Snakefile

## Instructions for Running
1. activate the environment containing snakemake and the plugins : conda activate snakemake
2. change to the directory containing the snakefile : cd path/to/snakefile
3. option to create environments before running the pipeline (otherwise they will be created in the next step) : snakemake --sdm conda --conda-create-envs-only -d path/to/project
4. run with one of the following commands
    - without SLURM : snakemake --sdm conda -d path/to/project -c #cores
    - with SLURM : snakemake --workflow-profile path/to/profiles -d path/to/project --scheduler greedy
5. to clean up temp files afterwards, if needed : snakemake --delete-temp-output -d path/to/project