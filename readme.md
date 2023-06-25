Little OTteR Assembly Pipeline (LOTRAP)
-------
A small pipeline created for the assembly of Otter genomes using reference guided assembly and longread sequences

current functionality:
* Ability to install reads and reference genomes using ID's
* A simple fastq analysis to view average quality and length of reads
* A script to perform VCF analysis 
* A line split program
* An interactive heatmap plot for the SNPs of the  VCF-file
* An interactive plot for the boxplot of the quality of each of the basepairs

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
* creating of a cladogram of multiple species
* annotation of contigs & visualisation of the annotation

Performance:
On the following data of the Streptococcus ruminantium
reads:      1.5Gb
reference:  2.1Mb

real    71m54.665s
user    51m29.261s
sys     23m7.101s

This was done using the time commando before the snakemake commando call on a Swift 3 Acer laptop.
Most of the runtime is from:
* Installing the sequences
* Mapping the reads onto the reference
* Creating the bam and VCF-file
