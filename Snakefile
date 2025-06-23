import os
configfile: "config.yaml"

SAMPLES = config["samples"]
READS_DIR = config["reads_dir"]
GENOME_INDEX = config["genome"]["index"]

rule all:
    input:
        expand("results/narrowpeaks/bigwig/{sample}.bw", sample=SAMPLES.keys())

rule fastqc:
    input:
        R1 = lambda wildcards: os.path.join(READS_DIR, SAMPLES[wildcards.sample]["R1"]),
        R2 = lambda wildcards: os.path.join(READS_DIR, SAMPLES[wildcards.sample]["R2"])
    output:
        "results/narrowpeaks/fastqc/{sample}_R1_fastqc.html",
        "results/narrowpeaks/fastqc/{sample}_R2_fastqc.html"
    shell:
        "fastqc {input.R1} {input.R2} -o results/narrowpeaks/fastqc"

rule align:
    input:
        R1 = lambda wildcards: os.path.join(READS_DIR, SAMPLES[wildcards.sample]["R1"]),
        R2 = lambda wildcards: os.path.join(READS_DIR, SAMPLES[wildcards.sample]["R2"])
    output:
        temp("results/narrowpeaks/aligned/{sample}.sam")
    shell:
        "bowtie2 -x {GENOME_INDEX} -1 {input.R1} -2 {input.R2} -S {output}"

rule sam_to_bam:
    input:
        "results/narrowpeaks/aligned/{sample}.sam"
    output:
        "results/narrowpeaks/bam/{sample}.bam"
    shell:
        "samtools view -bS {input} > {output}"

rule sort_index_bam:
    input:
        "results/narrowpeaks/bam/{sample}.bam"
    output:
        bam = "results/narrowpeaks/bam/{sample}.sorted.bam",
        bai = "results/narrowpeaks/bam/{sample}.sorted.bam.bai"
    shell:
        """
        samtools sort {input} -o {output.bam}
        samtools index {output.bam}
        """

rule bam_to_bigwig:
    input:
        "results/narrowpeaks/bam/{sample}.sorted.bam"
    output:
        "results/narrowpeaks/bigwig/{sample}.bw"
    shell:
        "bamCoverage -b {input} -o {output} --normalizeUsing RPKM"
