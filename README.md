# ATAC-seq pipeline with snakemake
###  Created 18-Feb-2021
###  Author: Brian Wray


## Before the pipeline, here are considerations for experiment design:

- two or more biological replicates
- each replicate has 25 million non-duplicate, non-mitochondrial aligned reads for 
  single-end sequencing and 50 million for paired-ended sequencing
- typically, no need for “input”
- use as few PCR cycles as possible when constructing the library
- paired-end sequencing is preferred

## Basic steps of ATACseq
[source] (https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/)

### Reads to Peaks
- trim reads with cutadapt
- Alignment with bowtie2
- sort with samtools
- filter out mito genes with samtools
- mark pcr duplicates with Picard's MarkDuplicates
- remove duplicates with samtools
- remove multi-mapped reads with samtools
- Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, 
  duplicates (-F 1804) with samtools
- Retain properly paired reads -f 2 with samtools
- (optional) check correlations between BAM files with deeptools multiBamSummary
- (optional) merge bams with samtools
- (optional) analyses that require single-base resolution e.g TF motif footprinting) shift reads with bedtools/awk
- (optional) filter for nucleosome-free regions
- Call peaks with MACS2 (use -B parameter for creating browser tracks)
- Annotate peaks with ChIPseeker or Homer

### QC
- Count fragments in bam files. Should have at least 25 million non-duplicate, non-mito fragments
- Check alignment rate, should be > 95%
- Measure library complexity with [Non-Redundant Fraction (NRF) and PCR Bottlenecking Coefficients 1 and 2, or PBC1 and PBC2](https://www.encodeproject.org/data-standards/terms/#library).
- Calculate IDR values ([Irreproducible Discovery Rate](https://www.encodeproject.org/data-standards/terms/#concordance)). Rescue and self consistency rations should be < 2
- Calculate TSS enrichment and FRiP score
- Plot fragment size distribution

