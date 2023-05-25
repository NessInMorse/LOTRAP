configfile: "env.yml"

folder = config["FOLDER"]
specimen = config["SPECIMEN"]
read_id = config["READS_ID"]
reference_id = config["REFERENCE_ID"]


rule all:
  input:
    f"{folder}{specimen}/",
    f"{folder}{specimen}/analysis/{specimen}_fastq_analysis.txt",
    f"{folder}{specimen}/reference/{specimen}.fasta",
    f"{folder}{specimen}/mapping/{specimen}.sam",
    f"{folder}{specimen}/mapping/{specimen}.bam",
    f"{folder}{specimen}/mapping/{specimen}.bai",
    f"{folder}{specimen}/variant_call/{specimen}.vcf",
    f"{folder}{specimen}/variant_call/{specimen}.vcf.gz",
    f"{folder}{specimen}/consensus/{specimen}.fa",
    f"{folder}{specimen}/analysis/analysis_vcf_{specimen}.txt",
    directory(f"{folder}{specimen}/variant_plots/")

rule install_fastq:
    output:
        f"{folder}{specimen}/reads/{specimen}.fastq"
    conda:
        "env.yml"
    shell:
        """
        fastq-dump {read_id} -O {folder}{specimen}/reads/ &&
        mv {folder}{specimen}/reads/{read_id}.fastq {folder}{specimen}/reads/{specimen}.fastq
        """

rule install_reference:
    output:
        f"{folder}{specimen}/reference/{specimen}.fasta"
    conda:
        "env.yml"
    shell:
        """
        efetch -db nuccore -id {reference_id} -format fasta > {folder}{specimen}/reference/{specimen}.fasta
        """


rule analyse_reads:
    input: f"{folder}{specimen}/reads/{specimen}.fastq"
    output: f"{folder}{specimen}/analysis/{specimen}_fastq_analysis.txt"
    conda:
        "env.yml"
    shell:
        "julia ./scripts/analyser.jl {folder}{specimen}/reads/{specimen}.fastq {folder}{specimen}/analysis/{specimen}_fastq_analysis.txt"


rule map_reads:
    input: f"{folder}{specimen}/reference/{specimen}.fasta",
           f"{folder}{specimen}/reads/{specimen}.fastq" 
    output:
        f"{folder}/{specimen}/mapping/{specimen}.sam"
    conda:
        "env.yml"
    threads: workflow.cores
    shell:
        "minimap2 -k 15 -w 30 -t {threads} -a {folder}{specimen}/reference/{specimen}.fasta {folder}{specimen}/reads/{specimen}.fastq > {folder}{specimen}/mapping/{specimen}.sam"

rule sam_bam:
    input:  
        f"{folder}{specimen}/mapping/{specimen}.sam"
    output:
        f"{folder}{specimen}/mapping/{specimen}.bam"
    conda:
        "env.yml"
    threads: workflow.cores * 0.25
    shell:
        """
        samtools view -bS {folder}{specimen}/mapping/{specimen}.sam > {folder}{specimen}/mapping/{specimen}.bam
        samtools sort -@ {threads} {folder}{specimen}/mapping/{specimen}.bam > {folder}{specimen}/mapping/file_sorted_bam
        mv {folder}{specimen}/mapping/file_sorted_bam {folder}{specimen}/mapping/{specimen}.bam
        """

rule create_bai:
    input:
        f"{folder}{specimen}/mapping/{specimen}.bam"
    output:
        f"{folder}{specimen}/mapping/{specimen}.bai"
    conda:
        "env.yml"
    threads: workflow.cores * 0.25
    shell:
        """
        samtools index -@ {threads} {folder}{specimen}/mapping/{specimen}.bam > {folder}{specimen}/mapping/{specimen}.bai
        """

rule create_vcf:
    input:
        f"{folder}{specimen}/reference/{specimen}.fasta",
        f"{folder}{specimen}/mapping/{specimen}.bam"
    output:
        f"{folder}{specimen}/variant_call/{specimen}.vcf"
    conda:
        "env.yml"
    shell:
        """
        bcftools mpileup -Ov -o {folder}{specimen}/variant_call/{specimen}.vcf -f {folder}{specimen}/reference/{specimen}.fasta {folder}{specimen}/mapping/{specimen}.bam 
        """

rule create_vcf_zipped:
    input:
        f"{folder}{specimen}/variant_call/{specimen}.vcf"
    output:
        f"{folder}{specimen}/variant_call/{specimen}.vcf.gz"
    conda:
        "env.yml"
    shell:
        """
        bgzip {folder}{specimen}/variant_call/{specimen}.vcf
        tabix -p vcf {folder}{specimen/variant_call}/{specimen}.vcf.gz
        """

rule create_consensus:
    input:
        f"{folder}{specimen}/reference/{specimen}.fasta",
        f"{folder}{specimen}/variant_call/{specimen}.vcf.gz"
    output:
        f"{folder}{specimen}/consensus/{specimen}.fa"
    conda:
        "env.yml"
    shell:
        """
        cat {folder}{specimen}/reference/{specimen}.fasta | bcftools consensus {folder}{specimen}/variant_call/{specimen}.vcf.gz > {folder}{specimen}/consensus/{specimen}.fa
        """

rule analyse_vcf:
    input:
        f"{folder}{specimen}/variant_call/{specimen}.vcf"
    output:
        f"{folder}{specimen}/analysis/analysis_vcf_{specimen}.txt"
    conda:
        "env.yml"
    shell:
        """
        julia scripts/vcf_analysis.jl {folder}{specimen}/variant_call/{specimen}.vcf {folder}{specimen}/analysis/analysis_vcf_{specimen}.txt
        """

rule plot_vcf:
    input:
        f"{folder}{specimen}/variant_call/{specimen}.vcf",
    output:
        directory(f"{folder}{specimen}/variant_plots/")
    shell:
        """
        mkdir {folder}{specimen}/variant_plots/
        julia scripts/vcf_plot.jl {folder}{specimen}/variant_call/{specimen}.vcf {folder}{specimen}/variant_plots/
        """