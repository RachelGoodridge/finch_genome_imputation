configfile: "config.yaml"
import pandas as pd

# register storage provider
if config["input_type"] == "remote_fs":
    storage:
        provider="fs"

# get a list of sample names from the .csv file
def list_samples():
    # read in the sample sheet
    reads = pd.read_csv(config["csv"])
    reads["run_num"] = reads["run_num"].astype(int)
    names = set(reads["name"].to_list())
    # check if there are any with duplicate sample names
    dups = reads[reads.duplicated(subset=["name", "run_num"])]
    if not dups.empty:
        raise Exception("Samples with the same identifier must have different run numbers.")
    return(names)

# find the corresponding srr given the sample name, if pulling from NCBI
def list_srr(wc):
    reads = pd.read_csv(config["csv"])
    reads["run_num"] = reads["run_num"].astype(int).astype(str)
    srr = reads[(reads["name"]==wc.sample) & (reads["run_num"]==wc.run)]["srr"].values[0]
    return(srr)

# find the corresponding read1 given the sample name
def list_read1(wc):
    if (config["input_type"] == "local") or (config["input_type"] == "remote_fs"):
        reads = pd.read_csv(config["csv"])
        reads["run_num"] = reads["run_num"].astype(int).astype(str)
        read1 = reads[(reads["name"]==wc.sample) & (reads["run_num"]==wc.run)]["read1"].values[0]
        if config["input_type"] == "remote_fs":
            return(storage.fs(read1))
        return(read1)
    elif config["input_type"] == "NCBI":
        return(f"raw_NCBI/{wc.sample}_run{wc.run}_R1.fastq.gz")

# find the corresponding read2 given the sample name
def list_read2(wc):
    if (config["input_type"] == "local") or (config["input_type"] == "remote_fs"):
        reads = pd.read_csv(config["csv"])
        reads["run_num"] = reads["run_num"].astype(int).astype(str)
        read2 = reads[(reads["name"]==wc.sample) & (reads["run_num"]==wc.run)]["read2"].values[0]
        if config["input_type"] == "remote_fs":
            return(storage.fs(read2))
        return(read2)
    elif config["input_type"] == "NCBI":
        return(f"raw_NCBI/{wc.sample}_run{wc.run}_R2.fastq.gz")

# figure out how many runs of each sample and list .bam files
def get_runs_bam(wc):
    reads = pd.read_csv(config["csv"])
    reads["run_num"] = reads["run_num"].astype(int)
    num_runs = reads[reads["name"] == wc.sample]["run_num"].tolist()
    return expand("results/1_mapped/{sample}_run{run}.bam", sample=wc.sample, run=num_runs)

# figure out how many runs of each sample and list .bai files
def get_runs_bai(wc):
    reads = pd.read_csv(config["csv"])
    reads["run_num"] = reads["run_num"].astype(int)
    num_runs = reads[reads["name"] == wc.sample]["run_num"].tolist()
    return expand("results/1_mapped/{sample}_run{run}.bam.bai", sample=wc.sample, run=num_runs)

# get a list of chunk ids from the chunks.txt file given the chromosome and list .vcfs to merge
def extract_ids_vcf(wc):
    chunks = checkpoints.chunk.get(chr_name=wc.chr_name).output[0]
    with open(chunks, "r") as file:
        ids = ["{:02d}".format(int(line.split()[0])) for line in file]
    return expand("results/3_imputed/{sample}/{chr_name}/{id_num}.vcf.gz", sample=wc.sample, chr_name=wc.chr_name, id_num=ids)

# get a list of chunk ids from the chunks.txt file given the chromosome and list .vcfs to merge
def extract_ids_csi(wc):
    chunks = checkpoints.chunk.get(chr_name=wc.chr_name).output[0]
    with open(chunks, "r") as file:
        ids = ["{:02d}".format(int(line.split()[0])) for line in file]
    return expand("results/3_imputed/{sample}/{chr_name}/{id_num}.vcf.gz.csi", sample=wc.sample, chr_name=wc.chr_name, id_num=ids)

# rule to call all the other rules
rule all:
    input:
        "results/depth_stats.csv",
        "results/final_output_filled.vcf.gz"

# index the reference genome
rule ref_bwa_index:
    input:
        config["ref"]  # reference genome
    output:
        expand("{ref}.{ext}", ref=config["ref"], ext=["amb", "ann", "bwt", "pac", "sa", "fai"])
    conda:
        "envs/environment.yaml"
    shell:
        """
        bwa index {input}
        samtools faidx {input}
        """

# get the data from NCBI if needed
rule get_data_ncbi:
    output:
        r1 = temp("raw_NCBI/{sample}_run{run}_R1.fastq.gz"),
        r2 = temp("raw_NCBI/{sample}_run{run}_R2.fastq.gz")
    threads: 8
    params:
        srr = list_srr  # get the SRR id corresponding to the sample
    conda:
        "envs/fasterq.yaml"
    shell:
        """
        prefetch --max-size 1T {params.srr}
        fasterq-dump {params.srr}/{params.srr}.sra -S -e {threads} -O raw_NCBI
        mv raw_NCBI/{params.srr}_1.fastq raw_NCBI/{wildcards.sample}_run{wildcards.run}_R1.fastq
        mv raw_NCBI/{params.srr}_2.fastq raw_NCBI/{wildcards.sample}_run{wildcards.run}_R2.fastq
        pigz -p {threads} raw_NCBI/{wildcards.sample}_run{wildcards.run}_R*.fastq
        rm -rf {params.srr}
        """

# clean up the reads with adapter trimming
rule fastp:
    input:
        r1 = list_read1,  # forward read
        r2 = list_read2  # reverse read
    output:
        r1 = temp("results/0_trimmed/{sample}_run{run}_R1.fastq.gz"),
        r2 = temp("results/0_trimmed/{sample}_run{run}_R2.fastq.gz"),
        summ = "logs/0_trimmed/{sample}_run{run}.fastp.out"
    threads: 8
    conda:
        "envs/fastp.yaml"
    log:
        "logs/0_trimmed/{sample}_run{run}.txt"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -w {threads} --detect_adapter_for_pe -j {output.summ} -h /dev/null &> {log}"

# map reads to the finch genome
rule bwa_mem_map:
    input:
        ref = config["ref"],  # reference genome
        ref_index = expand("{ref}.{ext}", ref=config["ref"], ext=["amb", "ann", "bwt", "pac", "sa", "fai"]),
        r1 = "results/0_trimmed/{sample}_run{run}_R1.fastq.gz",
        r2 = "results/0_trimmed/{sample}_run{run}_R2.fastq.gz"
    output:
        bam = temp("results/1_mapped/{sample}_run{run}.bam"),  # alignment file
        bam_index = temp("results/1_mapped/{sample}_run{run}.bam.bai")
    threads: 8
    params:
        rg = lambda wc: r"'@RG\tID:id\tSM:{sm}\tLB:lib\tPL:ILLUMINA'".format(sm=wc.sample)
    conda:
        "envs/environment.yaml"
    log:
        "logs/1_mapped/{sample}_run{run}.txt"
    shell:  # sort and convert to bam
        """
        bwa mem -M -t {threads} -R {params.rg} {input.ref} {input.r1} {input.r2} 2> {log} | samtools sort -o {output.bam}
        samtools index {output.bam} -o {output.bam_index}
        """

# merge bam files for multiple runs of the same sample
rule merge_bam:
    input:
        bam = get_runs_bam,
        bam_index = get_runs_bai
    output:
        bam = temp("results/1.0_merged/{sample}.bam"),
        bam_index = temp("results/1.0_merged/{sample}.bam.bai")
    threads: 8
    conda:
        "envs/environment.yaml"
    shell:
        """
        samtools merge --threads {threads} -o {output.bam} {input.bam}
        samtools index {output.bam} -o {output.bam_index}
        """

# separate into chromosomes
rule separate_bam:
    input:
        bam = "results/1.0_merged/{sample}.bam",  # alignment file
        bam_index = "results/1.0_merged/{sample}.bam.bai"
    output:
        bam = temp("results/1.1_split/{sample}/{chr_name}.bam"),
        bam_index = temp("results/1.1_split/{sample}/{chr_name}.bam.bai")
    conda:
        "envs/environment.yaml"
    shell:
        """
        samtools view -b {input.bam} {wildcards.chr_name} -o {output.bam}
        samtools index {output.bam} -o {output.bam_index}
        """

# estimate the sequencing depth by chromosome
rule read_depth:
    input:
        bam = "results/1.1_split/{sample}/{chr_name}.bam",
        bam_index = "results/1.1_split/{sample}/{chr_name}.bam.bai"
    output:
        "results/1.2_depth/{sample}/{chr_name}.txt"
    conda:
        "envs/environment.yaml"
    shell:
        "samtools depth {input.bam} -o {output}"

# calculate the average depth per sample and also per chromosome
rule avg_depth:
    input:
        expand("results/1.2_depth/{sample}/{chr_name}.txt", sample=list_samples(), chr_name=config["chrs"])
    output:
        "results/depth_stats.csv",
        "results/avg_depth_sample.csv",
        "results/avg_depth_chr.csv"
    script:
        "scripts/average_depth.py"

# index the entire panel first
rule index_panel:
    input:
        config["panel"]  # reference panel
    output:
        config["panel"] + ".csi"
    threads: 8
    conda:
        "envs/environment.yaml"
    shell:
        "bcftools index --threads {threads} -f {input}"

# split the reference panel into chromosomes
rule panel_split:
    input:
        vcf = config["panel"],  # reference panel
        csi = config["panel"] + ".csi"
    output:
        vcf = "panels/panel_{chr_name}/finches_panel_{chr_name}.vcf.gz",
        csi = "panels/panel_{chr_name}/finches_panel_{chr_name}.vcf.gz.csi"
    threads: 8
    conda:
        "envs/environment.yaml"
    shell:
        """
        bcftools view {input.vcf} --regions {wildcards.chr_name} -Oz -o {output.vcf}
        bcftools index --threads {threads} -f {output.vcf} -o {output.csi}
        """

# setup reference panel list of sites, extract, format, and zip
rule panel_sites_list:
    input:
        vcf = "panels/panel_{chr_name}/finches_panel_{chr_name}.vcf.gz",
        csi = "panels/panel_{chr_name}/finches_panel_{chr_name}.vcf.gz.csi"
    output:
        vcf = "panels/panel_{chr_name}/finches_panel_{chr_name}_sites.vcf.gz",
        tsv = "panels/panel_{chr_name}/finches_panel_{chr_name}_sites.tsv.gz"
    conda:
        "envs/environment.yaml"
    shell:
        """
        bcftools view -G -m 2 -M 2 -v snps {input.vcf} -Oz -o {output.vcf}
        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.vcf} | bgzip -c > {output.tsv}
        """

 # index the reference panel list of sites
rule index_panel_sites:
    input:
        vcf = "panels/panel_{chr_name}/finches_panel_{chr_name}_sites.vcf.gz",
        tsv = "panels/panel_{chr_name}/finches_panel_{chr_name}_sites.tsv.gz"
    output:
        vcf = "panels/panel_{chr_name}/finches_panel_{chr_name}_sites.vcf.gz.tbi",
        tsv = "panels/panel_{chr_name}/finches_panel_{chr_name}_sites.tsv.gz.tbi"
    conda:
        "envs/environment.yaml"
    shell:
        """
        tabix -s1 -b2 -e2 {input.vcf}
        tabix -s1 -b2 -e2 {input.tsv}
        """

# list potential sites for each read, call genotypes & generate likelihoods, index
rule sites_mpileup:
    input:
        ref = config["ref"],  # reference genome
        ref_index = expand("{ref}.{ext}", ref=config["ref"], ext=["amb", "ann", "bwt", "pac", "sa", "fai"]),
        vcf = "panels/panel_{chr_name}/finches_panel_{chr_name}_sites.vcf.gz",
        vcf_index = "panels/panel_{chr_name}/finches_panel_{chr_name}_sites.vcf.gz.tbi",
        tsv = "panels/panel_{chr_name}/finches_panel_{chr_name}_sites.tsv.gz",
        tsv_index = "panels/panel_{chr_name}/finches_panel_{chr_name}_sites.tsv.gz.tbi",
        bam = "results/1.1_split/{sample}/{chr_name}.bam",
        bam_index = "results/1.1_split/{sample}/{chr_name}.bam.bai"
    output:
        vcf = temp("results/2_mpileup/{sample}/{chr_name}.vcf.gz"),
        csi = temp("results/2_mpileup/{sample}/{chr_name}.vcf.gz.csi")
    threads: 8
    conda:
        "envs/environment.yaml"
    log:
        "logs/2_mpileup/{sample}/{chr_name}.txt"
    shell:
        """
        bcftools mpileup --threads {threads} -f {input.ref} -I -E -a 'FORMAT/DP' -T {input.vcf} {input.bam} -Ou 2> {log} |\
            bcftools call -Aim -C alleles -T {input.tsv} -Oz -o {output.vcf}
        bcftools index -f {output.vcf}
        """

# need to split into chromosome chunks before runnning GLIMPSE
checkpoint chunk:
    input:
        vcf = "panels/panel_{chr_name}/finches_panel_{chr_name}_sites.vcf.gz",
        vcf_index = "panels/panel_{chr_name}/finches_panel_{chr_name}_sites.vcf.gz.tbi"
    output:
        "panels/panel_{chr_name}/chunks.txt"
    conda:
        "envs/glimpse.yaml"
    log:
        "logs/0.1_chunk/chunks_{chr_name}.txt"
    shell:
        "GLIMPSE_chunk --input {input.vcf} --region {wildcards.chr_name} --window-size 2000000 --buffer-size 200000 --output {output} &> {log}"

# use GLIMPSE to do phasing and imputation refining genotype likelihoods
rule imputation:
    input:
        vcf = "results/2_mpileup/{sample}/{chr_name}.vcf.gz",
        csi = "results/2_mpileup/{sample}/{chr_name}.vcf.gz.csi",
        panel = "panels/panel_{chr_name}/finches_panel_{chr_name}.vcf.gz",  # reference panel
        panel_index = "panels/panel_{chr_name}/finches_panel_{chr_name}.vcf.gz.csi",
        gmap = config["gmap"] + "/{chr_name}.gmap",  # genetic map
        chunks = "panels/panel_{chr_name}/chunks.txt"
    output:
        vcf = temp("results/3_imputed/{sample}/{chr_name}/{id_num}.vcf.gz"),
        csi = temp("results/3_imputed/{sample}/{chr_name}/{id_num}.vcf.gz.csi")
    conda:
        "envs/glimpse.yaml"
    log:
        "logs/3_imputed/{sample}/{chr_name}/{id_num}.txt"
    shell:
        """
        LINE=$(awk "NR == $((10#{wildcards.id_num} + 1))" {input.chunks})
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)
        GLIMPSE_phase --main 15 --input {input.vcf} --reference {input.panel} --map {input.gmap} --input-region ${{IRG}} \
            --output-region ${{ORG}} --output {output.vcf} &> {log}
        bcftools index -f {output.vcf}
        """

# merge the imputation chunks back together
rule ligate:
    input:
        chunks = "panels/panel_{chr_name}/chunks.txt",
        vcf = extract_ids_vcf,
        csi = extract_ids_csi
    output:
        lst = temp("results/4_ligated/{sample}/{chr_name}.txt"),
        vcf = temp("results/4_ligated/{sample}/{chr_name}.vcf.gz"),
        csi = temp("results/4_ligated/{sample}/{chr_name}.vcf.gz.csi")
    conda:
        "envs/glimpse.yaml"
    log:
        "logs/4_ligated/{sample}/{chr_name}.txt"
    shell:
        """
        ls {input.vcf} > {output.lst}
        GLIMPSE_ligate --input {output.lst} --output {output.vcf} &> {log}
        bcftools index -f {output.vcf}
        """

# generate haplotype calls by sampling haplotype estimates
rule sample_hap:
    input:
        lst = "results/4_ligated/{sample}/{chr_name}.txt",
        vcf = "results/4_ligated/{sample}/{chr_name}.vcf.gz",
        csi = "results/4_ligated/{sample}/{chr_name}.vcf.gz.csi"
    output:
        vcf = temp("results/5_phased/{sample}/{chr_name}.vcf.gz"),
        csi = temp("results/5_phased/{sample}/{chr_name}.vcf.gz.csi")
    conda:
        "envs/glimpse.yaml"
    log:
        "logs/5_phased/{sample}/{chr_name}.txt"   
    shell:
        """
        GLIMPSE_sample --input {input.vcf} --solve --output {output.vcf} &> {log}
        bcftools index -f {output.vcf}
        """

# carry over the genotype probabilities from the imputation step
rule add_GP:
    input:
        ligate_lst = "results/4_ligated/{sample}/{chr_name}.txt",
        ligate_vcf = "results/4_ligated/{sample}/{chr_name}.vcf.gz",
        ligate_csi = "results/4_ligated/{sample}/{chr_name}.vcf.gz.csi",
        phase_vcf = "results/5_phased/{sample}/{chr_name}.vcf.gz",
        phase_csi = "results/5_phased/{sample}/{chr_name}.vcf.gz.csi"
    output:
        vcf = temp("results/5.1_addGP/{sample}/{chr_name}.vcf.gz"),
        csi = temp("results/5.1_addGP/{sample}/{chr_name}.vcf.gz.csi")
    conda:
        "envs/environment.yaml"
    shell:
        """
        bcftools annotate -a {input.ligate_vcf} -c FORMAT/GP {input.phase_vcf} -Oz -o {output.vcf}
        bcftools index -f {output.vcf}
        """

# merge all the chromosomes back together for each sample
rule merge_chrs:
    input:
        vcf = expand("results/5.1_addGP/{{sample}}/{chr_name}.vcf.gz", chr_name=config["chrs"]),
        csi = expand("results/5.1_addGP/{{sample}}/{chr_name}.vcf.gz.csi", chr_name=config["chrs"])
    output:
        vcf = "results/6_merged/{sample}.vcf.gz",
        csi = "results/6_merged/{sample}.vcf.gz.csi"
    conda:
        "envs/environment.yaml"
    shell:
        """
        bcftools concat {input.vcf} -Oz -o {output.vcf}
        bcftools index -f {output.vcf}
        """

# merge all the vcf files together into one final output
rule merge_all:
    input:
        vcf = expand("results/6_merged/{sample}.vcf.gz", sample=list_samples()),
        csi = expand("results/6_merged/{sample}.vcf.gz.csi", sample=list_samples())
    output:
        vcf = temp("results/final_output.vcf.gz"),
        filled = "results/final_output_filled.vcf.gz",
        csi = "results/final_output_filled.vcf.gz.csi",
        stats = "results/final_output_filled.stats"
    conda:
        "envs/environment.yaml"
    shell:
        """
        bcftools merge {input.vcf} -Oz -o {output.vcf}
        bcftools +fill-tags {output.vcf} -Oz -o {output.filled} -- -t all
        bcftools index -f {output.filled}
        bcftools stats {output.filled} > {output.stats}
        """