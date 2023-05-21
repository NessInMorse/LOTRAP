Little OTteR Assembly Pipeline (LOTRAP)
-------
A small pipeline created for the assembly of Otter genomes using reference guided assembly and longread sequences

current functionality:
* Ability to install reads and reference genomes using ID's
* A simple fastq analysis to view average quality and length of reads
* A script to perform VCF analysis 
* A line split program

Packages:
* samtools
* minimap2
* bcftools
* tabix
* bgzip
* fastq-dump
* efetch

Planned functionality:
* A protein predictor
* Interactive heatmap plot for the SNPs of the  VCF-file
* creating of a cladogram of multiple species
* annotation of contigs & visualisation of the annotation