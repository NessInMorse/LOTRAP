configfile: "env.yml"

folder = config["FOLDER"]
specimen = config["SPECIMEN"]
read_id = config["READS_ID"]
reference_id = config["REFERENCE_ID"]


rule all:
  input:
    f"{folder}{specimen}/",
    expand(f"{folder}/{specimen}/" + "{output_folders}/", output_folders=["reads", "ref_genome", "mapping", "consensus", "variant_call"]),
    f"{folder}{specimen}/analyses/{specimen}_fastq_analysis.txt",
    f"{folder}{specimen}/reference/{specimen}.fasta"

rule install_fastq:
    output:
        f"{folder}{specimen}/reads/{specimen}.fastq"
    shell:
        """
        fastq-dump {read_id} -O {folder}{specimen}/reads/ &&
        mv {folder}{specimen}/reads/{read_id}.fastq {folder}{specimen}/reads/{specimen}.fastq
        """

rule install_reference:
    output:
        f"{folder}{specimen}/reference/{specimen}.fasta"
    shell:
        """
        efetch -db nuccore -id {reference_id} -format fasta > {folder}{specimen}/reference/{specimen}.fasta
        """


rule analyse_reads:
    input: f"{folder}{specimen}/reads/{specimen}.fastq"
    output: f"{folder}{specimen}/analyses/{specimen}_fastq_analysis.txt"
    conda:
        "env.yml"
    shell:
        "julia ./scripts/analyser.jl {folder}{specimen}/reads/{specimen}.fastq {folder}{specimen}/analyses/{specimen}_fastq_analysis.txt"


rule map_reads:
    input: f"{folder}{specimen}/reference/{specimen}.fasta",
           f"{folder}{specimen}/reads/{specimen}.fastq" 
    output:
        f"{folder}/{specimen}/mapping/{specimen}.sam"
    shell:
        "minimap2 -k 15 -w 30 -t 8 -a {folder}{specimen}/reference/{specimen}.fasta {folder}{specimen}/reads/{specimen}.fastq > {folder}{specimen}/mapping/{specimen}.sam"

rule sam_bami:
    input:  "media/ness/PortableSSD/pipelines/genomes/lutra_lutra.sam"
    output:
        "file.bam",
    shell:
        """
        samtools view -bS alignment.sam > alignment.bam
        samtools sort -@ 8 file.bam > file_sorted_bam
        mv file_sorted_bam file.bam
        """
