Little OTteR Assembly Pipeline (LOTRAP)
-------
A small pipeline created for the assembly of Otter genomes using reference guided assembly and longread sequences

current functionality:
* A non-functioning snakemake file of the pipeline.
* various scripts to do FASTQ analysis, VCF analysis and split files based on lines.

Packages:
* samtools
* minimap2
* bcftools
* tabix
* bgzip

Planned functionality:
* A protein predictor
* Interactive heatmap plot for the SNPs of the  VCF-file
* creating of a cladogram of multiple species
* annotation of contigs & visualisation of the annotation