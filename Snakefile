configfile: "config.yaml"
import pandas as pd
import os

# register storage provider
if config["reads_remote_fs"]:
    storage:
        provider="fs"

# get a list of sample names from the .csv file
def list_samples():
    reads = pd.read_csv(config["csv"])
    names = reads["name"].tolist()
    return(names)

# find the corresponding read1 given the sample name
def list_read1(wc):
    reads = pd.read_csv(config["csv"])
    read1 = reads[reads["name"]==wc.sample]["read1"].values[0]
    if config["reads_remote_fs"]:
        return(storage.fs(read1))
    return(read1)

# find the corresponding read2 given the sample name
def list_read2(wc):
    reads = pd.read_csv(config["csv"])
    read2 = reads[reads["name"]==wc.sample]["read2"].values[0]
    if config["reads_remote_fs"]:
        return(storage.fs(read2))
    return(read2)

# get a list of chunk ids from the chunks.txt file given the chromosome and list .vcfs to merge
def extract_ids_vcfs(wc):
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

# get a list of chromosomes from the reference genome and list .vcfs to merge
def list_chrs_vcf(wc):
    cmd = f"grep '^>chr' {config['ref']} | cut -d ' ' -f 2"
    output = os.popen(cmd).read().splitlines()
    chrs = [out.replace(">", "") for out in output]
    return expand("results/5_phased/{sample}/{chr_name}.vcf.gz", sample=wc.sample, chr_name=chrs)

# get a list of chromosomes from the reference genome and list .csi files
def list_chrs_csi(wc):
    cmd = f"grep '^>chr' {config['ref']} | cut -d ' ' -f 2"
    output = os.popen(cmd).read().splitlines()
    chrs = [out.replace(">", "") for out in output]
    return expand("results/5_phased/{sample}/{chr_name}.vcf.gz.csi", sample=wc.sample, chr_name=chrs)

# get a list of chromosomes from the reference genome
def list_chrs():
    cmd = f"grep '^>chr' {config['ref']} | cut -d ' ' -f 2"
    output = os.popen(cmd).read().splitlines()
    chrs = [out.replace(">", "") for out in output]
    return(chrs)

# rule to call all the other rules
rule all:
    input:
        "results/depth_stats.csv",
        "results/avg_depth_sample.csv",
        "results/avg_depth_chr.csv",
        "results/final_output.vcf.gz",
        "results/final_output_filled.vcf.gz",
        "results/final_output_filled.vcf.gz.csi"

# index the reference genome
rule ref_bwa_index:
    input:
        config["ref"]  # reference genome
    output:
        expand("{ref}.{ext}", ref=config["ref"], ext=["amb", "ann", "bwt", "pac", "sa", "fai"])
    conda:
        config["envs"] + "/environment.yaml"
    shell:
        """
        bwa index {input}
        samtools faidx {input}
        """

# clean up the reads with adapter trimming
rule fastp:
    input:
        r1 = list_read1,  # forward read
        r2 = list_read2  # reverse read
    output:
        r1 = "results/0_trimmed/{sample}_R1.fastq.gz",
        r2 = "results/0_trimmed/{sample}_R2.fastq.gz",
        summ = "logs/0_trimmed/{sample}.fastp.out"
    threads: 8
    conda:
        config["envs"] + "/fastp.yaml"
    log:
        "logs/0_trimmed/{sample}.txt"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -w {threads} --detect_adapter_for_pe -j {output.summ} -h /dev/null &> {log}"

# map reads to the finch genome
rule bwa_mem_map:
    input:
        ref = config["ref"],  # reference genome
        ref_index = expand("{ref}.{ext}", ref=config["ref"], ext=["amb", "ann", "bwt", "pac", "sa", "fai"]),
        r1 = "results/0_trimmed/{sample}_R1.fastq.gz",
        r2 = "results/0_trimmed/{sample}_R2.fastq.gz"
    output:
        "results/1_mapped/{sample}.bam"  # alignment file
    threads: 8
    params:
        rg = lambda wc: r"'@RG\tID:id\tSM:{sm}\tLB:lib\tPL:ILLUMINA'".format(sm=wc.sample)
    conda:
        config["envs"] + "/environment.yaml"
    log:
        "logs/1_mapped/{sample}.txt"
    shell:  # sort and convert to bam
        "bwa mem -M -t {threads} -R {params.rg} {input.ref} {input.r1} {input.r2} 2> {log} | samtools sort -o {output}"

# index the sorted bam file
rule bam_samtools_index:
    input:
        "results/1_mapped/{sample}.bam"  # alignment file
    output:
        "results/1_mapped/{sample}.bam.bai"
    conda:
        config["envs"] + "/environment.yaml"
    shell:
        "samtools index {input} -o {output}"

# separate into chromosomes
rule separate_bam:
    input:
        bam = "results/1_mapped/{sample}.bam",  # alignment file
        bam_index = "results/1_mapped/{sample}.bam.bai"
    output:
        bam = "results/1.1_split/{sample}/{chr_name}.bam",
        bam_index = "results/1.1_split/{sample}/{chr_name}.bam.bai"
    conda:
        config["envs"] + "/environment.yaml"
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
        config["envs"] + "/environment.yaml"
    shell:
        "samtools depth {input.bam} -o {output}"

# calculate the average depth per sample and also per chromosome
rule avg_depth:
    input:
        expand("results/1.2_depth/{sample}/{chr_name}.txt", sample=list_samples(), chr_name=list_chrs())
    output:
        "results/depth_stats.csv",
        "results/avg_depth_sample.csv",
        "results/avg_depth_chr.csv"
    script:
        config["script"]

# index the entire panel first
rule index_panel:
    input:
        config["panel"]  # reference panel
    output:
        config["panel"] + ".csi"
    threads: 8
    conda:
        config["envs"] + "/environment.yaml"
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
        config["envs"] + "/environment.yaml"
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
        config["envs"] + "/environment.yaml"
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
        config["envs"] + "/environment.yaml"
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
        vcf = "results/2_mpileup/{sample}/{chr_name}.vcf.gz",
        csi = "results/2_mpileup/{sample}/{chr_name}.vcf.gz.csi"
    threads: 8
    conda:
        config["envs"] + "/environment.yaml"
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
        config["envs"] + "/glimpse.yaml"
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
        vcf = "results/3_imputed/{sample}/{chr_name}/{id_num}.vcf.gz",
        csi = "results/3_imputed/{sample}/{chr_name}/{id_num}.vcf.gz.csi"
    conda:
        config["envs"] + "/glimpse.yaml"
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
        vcf = extract_ids_vcfs,
        csi = extract_ids_csi
    output:
        lst = "results/4_ligated/{sample}/{chr_name}.txt",
        vcf = "results/4_ligated/{sample}/{chr_name}.vcf.gz",
        csi = "results/4_ligated/{sample}/{chr_name}.vcf.gz.csi"
    conda:
        config["envs"] + "/glimpse.yaml"
    log:
        "logs/4_ligated/{sample}/{chr_name}.txt"
    shell:
        #wc -l $CHUNKS | awk '{print $1}'
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
        vcf = "results/5_phased/{sample}/{chr_name}.vcf.gz",
        csi = "results/5_phased/{sample}/{chr_name}.vcf.gz.csi"
    conda:
        config["envs"] + "/glimpse.yaml"
    log:
        "logs/5_phased/{sample}/{chr_name}.txt"   
    shell:
        """
        GLIMPSE_sample --input {input.vcf} --solve --output {output.vcf} &> {log}
        bcftools index -f {output.vcf}
        """

# merge all the chromosomes back together for each sample
rule merge_chrs:
    input:
        vcf = list_chrs_vcf,
        csi = list_chrs_csi
    output:
        vcf = "results/6_merged/{sample}.vcf.gz",
        csi = "results/6_merged/{sample}.vcf.gz.csi"
    conda:
        config["envs"] + "/environment.yaml"
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
        vcf = "results/final_output.vcf.gz",
        filled = "results/final_output_filled.vcf.gz",
        csi = "results/final_output_filled.vcf.gz.csi"
    conda:
        config["envs"] + "/environment.yaml"
    shell:
        """
        bcftools merge {input.vcf} -Oz -o {output.vcf}
        bcftools +fill-tags {output.vcf} -Oz -o {output.filled} -- -t all
        bcftools index -f {output.filled}
        """