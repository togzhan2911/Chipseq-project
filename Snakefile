import os
configfile: "config.yaml"

SAMPLES = config["samples"]
READS_DIR = config["reads_dir"]
GENOME_INDEX = config["genome"]["index"]
GENOME_FASTA = config["genome"]["fasta"]

rule all:
    input:
        expand("results/narrowpeaks/bigwig/{sample}_{type}.bw", sample=SAMPLES.keys(), type=["ip", "input"]),
        expand("results/narrowpeaks/aligned/{sample}_{type}.bowtie2.log", sample=SAMPLES.keys(), type=["ip", "input"]),
        expand("results/narrowpeaks/bam/{sample}_{type}.flagstat.txt", sample=SAMPLES.keys(), type=["ip", "input"]),
        expand("results/narrowpeaks/macs2/{sample}_peaks.narrowPeak", sample=SAMPLES.keys()),
        expand("results/narrowpeaks/genrich/{sample}_genrich_peaks.narrowPeak", sample=SAMPLES.keys()),
        expand("results/narrowpeaks/fastqc/{sample}_{type}.done", sample=SAMPLES.keys(), type=["ip", "input"]),
        expand("results/narrowpeaks/peaks/{sample}_peaks.fasta", sample=SAMPLES.keys())

rule fastqc:
    input:
        R1 = lambda wc: os.path.join(READS_DIR, SAMPLES[wc.sample][wc.type]["R1"]),
        R2 = lambda wc: os.path.join(READS_DIR, SAMPLES[wc.sample][wc.type]["R2"])
    output:
        touch("results/narrowpeaks/fastqc/{sample}_{type}.done")
    params:
        outdir = "results/narrowpeaks/fastqc"
    shell:
        """
        fastqc {input.R1} {input.R2} -o {params.outdir} && touch {output}
        """

rule align:
    input:
        R1 = lambda wc: os.path.join(READS_DIR, SAMPLES[wc.sample][wc.type]["R1"]),
        R2 = lambda wc: os.path.join(READS_DIR, SAMPLES[wc.sample][wc.type]["R2"])
    output:
        sam = "results/narrowpeaks/aligned/{sample}_{type}.sam",
        log = "results/narrowpeaks/aligned/{sample}_{type}.bowtie2.log"
    shell:
        """
        bowtie2 -x {GENOME_INDEX} -1 {input.R1} -2 {input.R2} \
        -S {output.sam} 2> {output.log}
        """

rule sam_to_bam:
    input:
        "results/narrowpeaks/aligned/{sample}_{type}.sam"
    output:
        "results/narrowpeaks/bam/{sample}_{type}.bam"
    wildcard_constraints:
        sample = "(?!.*\\.sorted).*",
        type = "(ip|input)"
    shell:
        """
        samtools view -bS {input} > {output}
        """

rule sort_index_bam:
    input:
        "results/narrowpeaks/bam/{sample}_{type}.bam"
    output:
        bam = "results/narrowpeaks/bam/{sample}_{type}.sorted.bam",
        bai = "results/narrowpeaks/bam/{sample}_{type}.sorted.bam.bai"
    shell:
        """
        samtools sort {input} -o {output.bam}
        samtools index {output.bam}
        """

rule bam_to_bigwig:
    input:
        "results/narrowpeaks/bam/{sample}_{type}.sorted.bam"
    output:
        "results/narrowpeaks/bigwig/{sample}_{type}.bw"
    shell:
        """
        bamCoverage -b {input} -o {output} --normalizeUsing RPKM
        """

rule bam_stats:
    input:
        "results/narrowpeaks/bam/{sample}_{type}.sorted.bam"
    output:
        "results/narrowpeaks/bam/{sample}_{type}.flagstat.txt"
    shell:
        """
        samtools flagstat {input} > {output}
        """

rule macs2_callpeak:
    input:
        chip = "results/narrowpeaks/bam/{sample}_ip.sorted.bam",
        ctrl = "results/narrowpeaks/bam/{sample}_input.sorted.bam"
    output:
        "results/narrowpeaks/macs2/{sample}_peaks.narrowPeak"
    params:
        name = lambda wc: wc.sample
    shell:
        """
        macs2 callpeak -t {input.chip} -c {input.ctrl} \
        -f BAMPE -g 1.8e7 -n {params.name} \
        --outdir results/narrowpeaks/macs2 --keep-dup all
        """

rule genrich_callpeak:
    input:
        chip = "results/narrowpeaks/bam/{sample}_ip.sorted.bam",
        ctrl = "results/narrowpeaks/bam/{sample}_input.sorted.bam"
    output:
        "results/narrowpeaks/genrich/{sample}_genrich_peaks.narrowPeak"
    shell:
        """
        genrich -t {input.chip} -c {input.ctrl} -o {output} -f narrow
        """

rule extract_peak_sequences:
    input:
        genome = GENOME_FASTA,
        peaks = "results/narrowpeaks/macs2/{sample}_peaks.narrowPeak"
    output:
        fasta = "results/narrowpeaks/peaks/{sample}_peaks.fasta"
    shell:
        """
        bedtools getfasta -fi {input.genome} -bed {input.peaks} -fo {output.fasta}
        """

