# Simple ChIP-seq Pipeline

This is a simple ChIP-seq analysis pipeline built with Nextflow. It uses local files and the following tools:
* FastQC
* BWA
* Samtools
* MACS2
* deepTools

## 1. Setup

This pipeline does not download data. You must prepare your local files first.

### Genome

1.  Place your `hg38.fa` file in the `genomee/` folder.
2.  Build the BWA index: `bwa index genomee/hg38.fa`

### Samples

1.  Place your paired-end `_1.fastq.gz` and `_2.fastq.gz` files in the `dataa/` folder.

## 2. Dependencies

This pipeline requires all tools to be in your `PATH`. A Conda environment is recommended:

```bash
mamba create -n chipseq_env -c bioconda -c conda-forge \
  openjdk nextflow bwa samtools fastqc macs2 deeptools
mamba activate chipseq_env
