folder = "/media/ness/PortableSSD/pipelines/"
specimen = "enhydra_lutris_kenyoni"

reference_id = ""
reads_id = ""
max_corecount = 8


rule all:
  input:
    expand("output/lutra_lutra_1_{seq}.txt", seq=["nanopore", "pac_bio"])


rule start:
    output: "{folder}{specimen}/reads",
            "{folder}{specimen}/ref_genome",
            "{folder}{specimen}/variant_call",
            "{folder}{specimen}/mapping",
            "{folder}{specimen}/consensus"
    shell:
        """
        mkdir {folder}{specimen}/reads
        mkdir {folder}{specimen}/ref_genome
        mkdir {folder}{specimen}/variant_call
        mkdir {folder}{specimen}/mapping
        mkdir {folder}{specimen}/consensus
        """

rule analyse_reads:
    input: "input/lutra_lutra_1_{seq}.fastq"
    output: "output/lutra_lutra_1_{seq}.txt"
    conda:
        "env.yml"
    shell:
        "julia ./scripts/analyser.jl {input} {output}"

rule map_reads:
    input: "/media/ness/PortableSSD/pipelines/ref_genome/lutra_lutra_genome.fna",
           "/media/ness/PortableSSD/pipelines/reads/lutra_lutra_nanopore.fastq" 
    output:
        "media/ness/PortableSSD/pipelines/genomes/lutra_lutra.sam"
    shell:
        "minimap2 -k 15 -w 30 -t 8 -a /media/ness/PortableSSD/pipelines/ref_genomes/lutra_lutra_genome.fna /media/ness/PortableSSD/pipelines/reads/lutra_lutra_nanopore.fastq > media/ness/PortableSSD/pipelines/genomes/lutra_lutra.sam"

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
