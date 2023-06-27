configfile: "env.yml"


folder = "/media/ness/PortableSSD/pipelines/"
specimen = "lutra_lutra"
read_id = "ERR3313341"
reference_name = "lutra lutra"


rule all:
  input:
    f"{folder}{specimen}/",
    f"{folder}{specimen}/notes/",
    f"{folder}{specimen}/analysis/{specimen}_fastq_analysis.txt",
    f"{folder}{specimen}/analysis/{specimen}_readqc.png",
    f"{folder}{specimen}/analysis/{specimen}_readqc.html",
    f"{folder}{specimen}/reference/{specimen}.fasta",
    f"{folder}{specimen}/mapping/{specimen}.sam",
    f"{folder}{specimen}/mapping/{specimen}.bam",
    f"{folder}{specimen}/mapping/{specimen}.bai",
    f"{folder}{specimen}/variant_call/{specimen}_cp.vcf",
    f"{folder}{specimen}/variant_call/{specimen}.vcf.gz",
    f"{folder}{specimen}/consensus/{specimen}.fa",
    f"{folder}{specimen}/analysis/analysis_vcf_{specimen}.txt",
    f"{folder}{specimen}/variant_plots/",
    f"{folder}{specimen}/script_times/",
    f"{folder}{specimen}/script_times/analyser.tsv",
    f"{folder}{specimen}/script_times/fastq_analyser_plot.tsv",
    f"{folder}{specimen}/script_times/vcf_analysis.tsv",
    f"{folder}{specimen}/script_times/vcf_plot.tsv",
    f"{folder}{specimen}/script_times/summary.tsv"
    

rule create_notes:
    output:
        directory(f"{folder}{specimen}/notes/")
    shell:
        """
        mkdir {folder}{specimen}/notes/
        echo "species: {specimen}" > {folder}{specimen}/notes/notes.txt
        echo "read id: {read_id}" >> {folder}{specimen}/notes/notes.txt
        echo "reference name: {reference_name}" >> {folder}{specimen}/notes/notes.txt
        """

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
        cd {folder}{specimen}/reference/
        datasets download genome taxon {reference_name} --reference --filename reference.zip
        gunzip reference.zip
        folder=$(ls reference/ncbi_dataset/data | egrep "\.[0-9]" | awk '{{print($NF)}}')
        file=$(ls reference/ncbi_dataset/data/$folder)
        mv reference/ncbi_dataset/data/$folder/$file {specimen}.fasta
        """


rule analyse_reads:
    input: f"{folder}{specimen}/reads/{specimen}.fastq"
    output: f"{folder}{specimen}/analysis/{specimen}_fastq_analysis.txt",
            f"{folder}{specimen}/script_times/analyser.tsv"
    conda:
        "env.yml"
    shell:
        "julia ./scripts/analyser.jl {folder}{specimen}/reads/{specimen}.fastq {folder}{specimen}/analysis/{specimen}_fastq_analysis.txt {folder}{specimen}/script_times/analyser.tsv"


rule plot_quality_reads:
    input: 
        f"{folder}{specimen}/reads/{specimen}.fastq"
    output:
        f"{folder}{specimen}/analysis/{specimen}_readqc.png",
        f"{folder}{specimen}/analysis/{specimen}_readqc.html",
        f"{folder}{specimen}/script_times/fastq_analyser_plot.tsv"
    shell:
        """
        julia ./scripts/fastq_analyser_plot.jl {folder}{specimen}/reads/{specimen}.fastq {folder}{specimen}/analysis/{specimen}_readqc {folder}{specimen}/script_times/fastq_analyser_plot.tsv
        """


rule map_reads:
    input: f"{folder}{specimen}/reference/{specimen}.fasta",
           f"{folder}{specimen}/reads/{specimen}.fastq" 
    output:
        f"{folder}{specimen}/mapping/{specimen}.sam"
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
    threads:
        workflow.cores * 0.25
    shell:
        """
        bcftools mpileup -Ov --threads {threads} -f {folder}{specimen}/reference/{specimen}.fasta {folder}{specimen}/mapping/{specimen}.bam | bcftools call -mv -Ov --threads {threads} > {folder}{specimen}/variant_call/{specimen}.vcf
        """

rule create_vcf_zipped:
    input:
        f"{folder}{specimen}/variant_call/{specimen}.vcf"
    output:
        f"{folder}{specimen}/variant_call/{specimen}_cp.vcf",
        f"{folder}{specimen}/variant_call/{specimen}.vcf.gz"
    conda:
        "env.yml"
    shell:
        """
        cp {folder}{specimen}/variant_call/{specimen}.vcf {folder}{specimen}/variant_call/{specimen}_cp.vcf
        bgzip {folder}{specimen}/variant_call/{specimen}.vcf
        tabix -p vcf {folder}{specimen}/variant_call/{specimen}.vcf.gz
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
        f"{folder}{specimen}/variant_call/{specimen}_cp.vcf"
    output:
        f"{folder}{specimen}/analysis/analysis_vcf_{specimen}.txt",
        f"{folder}{specimen}/script_times/vcf_analysis.tsv"
    conda:
        "env.yml"
    shell:
        """
        julia scripts/vcf_analysis.jl {folder}{specimen}/variant_call/{specimen}_cp.vcf {folder}{specimen}/analysis/analysis_vcf_{specimen}.txt {folder}{specimen}/script_times/vcf_analysis.tsv
        """

rule plot_vcf:
    input:
        f"{folder}{specimen}/variant_call/{specimen}_cp.vcf",
    output:
        directory(f"{folder}{specimen}/variant_plots/"),
        f"{folder}{specimen}/script_times/vcf_plot.tsv"
    shell:
        """
        mkdir {folder}{specimen}/variant_plots/
        julia scripts/vcf_plot.jl {folder}{specimen}/variant_call/{specimen}_cp.vcf {folder}{specimen}/variant_plots/ {folder}{specimen}/script_times/vcf_plot.tsv
        """

rule script_time:
    input:
        f"{folder}{specimen}/script_times/analyser.tsv",
        f"{folder}{specimen}/script_times/fastq_analyser_plot.tsv",
        f"{folder}{specimen}/script_times/vcf_analysis.tsv",
        f"{folder}{specimen}/script_times/vcf_plot.tsv",
    output:
        f"{folder}{specimen}/script_times/summary.tsv"
    shell:
        """
        julia scripts/summary.jl {folder}{specimen}/script_times/ {folder}{specimen}/script_times/summary.tsv
        """