# Isoform Project
**Goal:** quantify and describe the presence of different gene isoforms (variant transcripts of a gene produced by alternative splicing) in scRNA-seq data

[R Notebook for Visualization & Quantification](https://github.com/daviemel/isoform_project/blob/main/sample_viz.Rmd):
* Creates random sampling of FASTQ files for cells within a specific subclass (for input in [STAR script](https://github.com/daviemel/isoform_project/blob/main/STAR_align_l23.sh))
* Merges BAM files (the outputs of [STAR script](https://github.com/daviemel/isoform_project/blob/main/STAR_align_l23.sh)) for later visualization and quantification at the subclass-level
* Uses Bioconductor functions to set up and plot a figure visualizing the following: 1) read coverage, 2) junctions (i.e., sashimi), and 3) gene isoform models

[Example STAR Script for Read Alignment](https://github.com/daviemel/isoform_project/blob/main/STAR_align_l23.sh)
* Given FASTQ files, processed genome index directory, and genome annotation file (.gtf or .gff), creates scripts which will produce BAMs
