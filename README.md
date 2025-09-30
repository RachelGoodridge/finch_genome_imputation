# Finch Genome Imputation
GLIMPSE pipeline for genomic imputation and variant calling

## Requirements for Use
- install snakemake : conda create -c bioconda -c conda-forge -n snakemake snakemake=9.1.10
- install plugins (in the same environment) : conda install -c bioconda -c conda-forge snakemake-executor-plugin-slurm snakemake-storage-plugin-fs
- install biopython (in the same environment) : conda install -c conda-forge biopython=1.85
- environment.yaml : must be inside environment folder
- fasterq.yaml : must be inside environment folder
- fastp.yaml : must be inside environment folder
- glimpse.yaml : must be inside environment folder
- profiles/config.yaml : must be inside a folder called "profiles" and contain the slurm account, if using slurm
- project/config.yaml : must be inside project folder and contain full paths
- samples.csv : must contain a unqiue identifier for each sample and full path to the data, data may be stored remotely
- average_depth.py
- reference genome : a .fna file
- reference panel : a .vcf file
- genetic maps : must be named "chr1.gmap", "chr2.gmap", etc. and contained in a single folder
- Snakefile

## Instructions for Running
1. activate the environment containing snakemake and the plugins : conda activate snakemake
2. change to the directory containing the snakefile : cd path/to/snakefile
3. option to create environments before running the pipeline (otherwise they will be created in the next step) :  snakemake --sdm conda --conda-create-envs-only -d path/to/project
4. run with one of the following commands
    - without SLURM : snakemake --sdm conda -d path/to/project -c #cores
    - with SLURM : snakemake --workflow-profile path/to/profiles -d path/to/project --scheduler greedy
5. to clean up temp files afterwards, if needed : snakemake --delete-temp-output -d path/to/project