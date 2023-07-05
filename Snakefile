configfile: "env.yml"


folder = "/media/ness/PortableSSD/pipelines/"
specimen = "Salmonella_enterica"
read_id = "SRR24883431"
reference_name = "28901"


rule all:
  input:
    f"{folder}{specimen}/analysis/{specimen}_fastq_analysis.txt",
    f"{folder}{specimen}/analysis/{specimen}_readqc.png",
    f"{folder}{specimen}/analysis/{specimen}_readqc.html",
    f"{folder}{specimen}/analysis/{specimen}_gc_distribution.png",
    f"{folder}{specimen}/analysis/{specimen}_gc_distribution.html",
    f"{folder}{specimen}/reference/{specimen}.fasta",
    f"{folder}{specimen}/mapping/{specimen}.sam",
    f"{folder}{specimen}/mapping/{specimen}.bam",
    f"{folder}{specimen}/mapping/{specimen}.bai",
    f"{folder}{specimen}/variant_call/{specimen}_cp.vcf",
    f"{folder}{specimen}/variant_call/{specimen}.vcf.gz",
    f"{folder}{specimen}/consensus/{specimen}.fa",
    f"{folder}{specimen}/analysis/analysis_vcf_{specimen}.txt",
    f"{folder}{specimen}/script_times/analyser.tsv",
    f"{folder}{specimen}/script_times/fastq_analyser_plot.tsv",
    f"{folder}{specimen}/script_times/vcf_analysis.tsv",
    f"{folder}{specimen}/script_times/vcf_plot.tsv",
    f"{folder}{specimen}/script_times/summary.tsv"


rule create_notes:
    output: notes = f"{folder}{specimen}/notes/notes.txt"
    shell:
        """
        echo "species: {specimen}" > {output.notes}
        echo "read id: {read_id}" >> {output.notes}
        echo "reference name: {reference_name}" >> {output.notes}
        """

rule install_fastq:
    output:
        reads = f"{folder}{specimen}/reads/{specimen}.fastq"
    conda:
        "env.yml"
    shell:
        """
        fastq-dump {read_id} -O {folder}{specimen}/reads/ &&
        mv {folder}{specimen}/reads/{read_id}.fastq {output.reads}
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
        unzip reference.zip
        folder=$(ls ncbi_dataset/data | egrep "\.[0-9]" | awk '{{print($NF)}}')
        file=$(ls ncbi_dataset/data/$folder)
        mv ncbi_dataset/data/$folder/$file {specimen}.fasta
        """


rule analyse_reads:
    input: reads = f"{folder}{specimen}/reads/{specimen}.fastq"
    output: fastq_analysis = f"{folder}{specimen}/analysis/{specimen}_fastq_analysis.txt",
            time_file = f"{folder}{specimen}/script_times/analyser.tsv"
    conda:
        "env.yml"
    shell:
        """
        julia ./scripts/analyser.jl {input.reads} {output.fastq_analysis} {output.time_file}
        """


rule plot_quality_reads:
    input: 
        reads = f"{folder}{specimen}/reads/{specimen}.fastq"
    output:
        f"{folder}{specimen}/analysis/{specimen}_readqc.png",
        f"{folder}{specimen}/analysis/{specimen}_readqc.html",
        f"{folder}{specimen}/analysis/{specimen}_gc_distribution.png",
        f"{folder}{specimen}/analysis/{specimen}_gc_distribution.html",
        f"{folder}{specimen}/script_times/fastq_analyser_plot.tsv"
    shell:
        """
        julia ./scripts/fastq_analyser_plot.jl {input.reads} {folder}{specimen}/analysis/{specimen}_readqc {folder}{specimen}/analysis/{specimen}_gc_distribution {folder}{specimen}/script_times/fastq_analyser_plot.tsv
        """

rule map_reads:
    input: reference = f"{folder}{specimen}/reference/{specimen}.fasta",
           reads = f"{folder}{specimen}/reads/{specimen}.fastq" 
    output:
           samfile = f"{folder}{specimen}/mapping/{specimen}.sam"
    conda:
        "env.yml"
    threads: workflow.cores
    shell:
        "minimap2 -k 15 -w 30 -t {threads} -a {input.reference} {input.reads} > {output.samfile}"

rule sam_bam:
    input:  
        samfile = f"{folder}{specimen}/mapping/{specimen}.sam"
    output:
        bamfile = f"{folder}{specimen}/mapping/{specimen}.bam"
    conda:
        "env.yml"
    threads: workflow.cores * 0.25
    shell:
        """
        samtools view -bS {input.samfile} > {output.bamfile}
        samtools sort -@ {threads} {output.bamfile} > {folder}{specimen}/mapping/file_sorted_bam
        mv {folder}{specimen}/mapping/file_sorted_bam {output.bamfile}
        """

rule create_bai:
    input:
        bamfile = f"{folder}{specimen}/mapping/{specimen}.bam"
    output:
        baifile = f"{folder}{specimen}/mapping/{specimen}.bai"
    conda:
        "env.yml"
    threads: workflow.cores * 0.25
    shell:
        """
        samtools index -@ {threads} {input.bamfile} > {output.baifile}
        """

rule create_vcf:
    input:
        reference = f"{folder}{specimen}/reference/{specimen}.fasta",
        bamfile = f"{folder}{specimen}/mapping/{specimen}.bam"
    output:
        variant_call = f"{folder}{specimen}/variant_call/{specimen}.vcf"
    conda:
        "env.yml"
    threads:
        workflow.cores * 0.25
    shell:
        """
        bcftools mpileup -Ov -f {input.reference} {input.bamfile} | bcftools call -mv -Ov > {output.variant_call}
        """

rule create_vcf_zipped:
    input:
        variant_call = f"{folder}{specimen}/variant_call/{specimen}.vcf"
    output:
        variant_copy = f"{folder}{specimen}/variant_call/{specimen}_cp.vcf",
        variant_zipped = f"{folder}{specimen}/variant_call/{specimen}.vcf.gz"
    conda:
        "env.yml"
    shell:
        """
        cp {input.variant_call} {output.variant_copy}
        bgzip {input.variant_call}
        tabix -p vcf {output.variant_zipped}
        """

rule create_consensus:
    input:
        reference = f"{folder}{specimen}/reference/{specimen}.fasta",
        variant_zipped = f"{folder}{specimen}/variant_call/{specimen}.vcf.gz"
    output:
        consensus = f"{folder}{specimen}/consensus/{specimen}.fa"
    conda:
        "env.yml"
    shell:
        """
        cat {input.reference} | bcftools consensus {input.variant_zipped} > {output.consensus}
        """

rule analyse_vcf:
    input:
        variant_copy = f"{folder}{specimen}/variant_call/{specimen}_cp.vcf"
    output:
        vcf_text = f"{folder}{specimen}/analysis/analysis_vcf_{specimen}.txt",
        script_time = f"{folder}{specimen}/script_times/vcf_analysis.tsv"
    conda:
        "env.yml"
    shell:
        """
        julia scripts/vcf_analysis.jl {input.variant_copy} {output.vcf_text} {output.script_time}
        """

rule plot_vcf:
    input:
        variant_copy = f"{folder}{specimen}/variant_call/{specimen}_cp.vcf",
    output:
        script_time = f"{folder}{specimen}/script_times/vcf_plot.tsv"
    shell:
        """
        julia scripts/vcf_plot.jl {input.variant_copy} {folder}{specimen}/analysis/ {output.script_time}
        """

rule script_time:
    input:
        f"{folder}{specimen}/script_times/analyser.tsv",
        f"{folder}{specimen}/script_times/fastq_analyser_plot.tsv",
        f"{folder}{specimen}/script_times/vcf_analysis.tsv",
        f"{folder}{specimen}/script_times/vcf_plot.tsv",
    output:
        summary_time = f"{folder}{specimen}/script_times/summary.tsv"
    shell:
        """
        julia scripts/summary.jl {folder}{specimen}/script_times/ {output.summary_time}
        """