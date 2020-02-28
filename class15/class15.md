Class 15: Transcriptomics and the Analysis of RNA-Seq Data
================

## Background

The data for this hands-on session comes from a published RNA-seq
experiment where airway smooth muscle cells were treated with
dexamethasone, a synthetic glucocorticoid steroid with anti-inflammatory
effects [Himes et
al. 2014](http://www.ncbi.nlm.nih.gov/pubmed/24926665).

Glucocorticoids are used, for example, by people with asthma to reduce
inflammation of the airways. The anti-inflammatory effects on airway
smooth muscle (ASM) cells has been known for some time but the
underlying molecular mechanisms are unclear.

## 1\. Bioconductor and DESeq2 Setup

Bioconductor is a large repository and resource for R packages that
focus on analysis of high-throughput genomic data.

To install Bioconductor, we need to use the following commands in
console: `install.packages("BiocManager")` `BiocManager::install()`

We will also need DESeq2: `BiocManager::install("DESeq2")`

### Side-note: aligning reads to a reference genome

The computational analysis of an RNA-seq experiment begins from the
FASTQ files that contain the nucleotide sequence of each read and a
quality score at each position. These reads must first be aligned to a
reference genome or transcriptome. The output of this alignment step is
commonly stored in a file format called SAM/BAM.

### DESeq2 required inputs

As input, the DESeq2 package expects **(1)** a data.frame of **count
data** (as obtained from RNA-seq or another high-throughput sequencing
experiment) and **(2)** a second data.frame with information about the
samples - often called sample metadata (or `colData` in DESeq2-speak
because it supplies metadata/information about the columns of the
countData matrix).
