# Overview of mapping sequences to a reference genome

* Teaching: 30 minutes

#### Objectives

* Understand the fundamental process of mapping reads to a reference sequence.
* Understand the output files from mapping tools.
* Learn how to import mapping results into GUI programs such as `Geneious` for further analysis.

#### Keypoints

* Reads can be mapped to any reference sequence, for example whole genomes or individual genes, to generate a consensus sequence or view the coverage in particular regions.
* Mapping alignments are stored in the `sam` (text, human-readable) or `bam` (binary, compressed) formats.
* Different mapping tools work in different ways and make different assumptions when working with your data.
  * Sometimes, a given tool will work better with a given dataset.
  * It is recommended to try a few tools when working with new datasets to see which one produces better results.

---

## Contents

1. [Introduction to mapping](#introduction-to-mapping)
1. [Overview of the mapping process](#overview-of-the-mapping-process)
1. [Storing read mapping information](#storing-read-mapping-information)
1. [Different tools for the job](#different-tools-for-the-job)
1. [Obtaining a reference sequence](#obtaining-a-reference-sequence)

---

## Introduction to mapping

Sequencing produces a collection of sequences without genomic context - as the sequencing process reads genomic fragments in a random order we do not know the region of the genome to which each sequence corresponds. Mapping the reads of an experiment to a reference genome or a reference gene is a key step for resolving this issue. After mapping, the reads are assigned to a specific location in the genome/gene and a consensus file can be generated that reconstructs the sequence of that region on the organism that provided the reads.

During sequencing, errors are introduced, such as incorrect nucleotides being called. Sequencing errors might bias the analysis and can lead to a misinterpretation of the data and erroneous mapping. It is good practice to filter out low quality data from our dataset before mapping, as we learned in the previous session.

Read mapping is the process of aligning reads to a reference sequence. A mapping tool takes as input a reference genome and a set of reads. Its aim is to align each read in the set of reads on the reference genome, allowing mismatches, indels and clipping of some short fragments on the two ends of the reads.

---

## Overview of the mapping process

When mapping sequences, the input consists of a set of reads and a reference genome. when mapping reads to a reference the matches can be a perfect alignment between the read and reference, but often the match will be inexact. This is generally due to one of two reasons:

1. The reference genome and the sequenced organisms are not genetically identical, either due to some individual-specific mutations or potentially that there is not a closely related reference genome for use in the mapping process.
1. Sequencing error.

![Mapping overview](../img/level2_31_mapping.png)

The figure above gives an example of the ways in which reads can map with mismatches to a reference. The first read is aligned at position 100 and the alignment has two mismatches. The second read is aligned at position 114. It is a local alignment with clippings on the left and right. The third read is aligned at position 123. It consists of a 2-base insertion and a 1-base deletion.

When working with mapped reads it is our job to examine these areas of mismatch and evaluate whether or not these represent real sequence variation.

---

## Storing read mapping information

Mapping tools usually produce their output in the `sam` format. This file stores the read sequences, a flag on whether they have been aligned to a reference sequence, and if so, the position on the reference sequence at which they have been aligned. You can view this file using any text viewer, although owing to the file size `less` is generally the best choice. Should you choose to do so, a `sam` file is basically a table of the alignment between each read and the reference and consists of the following columns:

|Column|Field|Description|
|:---|:---|:---|
|1|QNAME|Query template (read) name|
|2|FLAG|Bitwise flag|
|3|RNAME|References sequence name|
|4|POS|1- based leftmost mapping position|
|5|MAPQ|Mapping quality|
|6|CIGAR|CIGAR string|
|7|RNEXT|Ref. name of the mate/next read|
|8|PNEXT|Position of the mate/next read|
|9|TLEN|Observed template length|
|10|SEQ|Segment sequence|
|11|QUAL|ASCII representation of Phred-scaled base quality|

The full specification for the file format [can be found here](https://samtools.github.io/hts-specs/SAMv1.pdf).

However, it is not usually necessary to inspect the `sam` file directly - there is a lot of information in here and unless you are looking to extract specific information from the alignment it is just an intermediate file in the workflow. In order to save disk space, and to prepare the file for downstream analysis, `sam` files are ussually sorted and compressed.

Sorting the mapping information is an important prerequisite for performing certain downstream processes. Not every tool we use requires reads to be sorted, but it can be frustrating having to debug the instances where read sorting matters, so we typically just get it done as soon as possible and then we don't have to worry about it again. Reads will initially be mapped in an unsorted order, as they are added to the `sam` file in more or less the same order as they are encountered in the input fastq files. When sorted the reads are ordered by the position of their first mapped nucleotide, as exemplified below:

```
# Unsorted reads
Ref: REFERENCECONTIG
Map: --------ECONT--
Map: REFE-----------
Map: --FERENCECO----
Map: -----------NTIG

# Sorted reads
Ref: REFERENCECONTIG
Map: REFE-----------
Map: --FERENCECO----
Map: --------ECONT--
Map: -----------NTIG
```

Compressing the file to the `bam` format is an important step as `sam` files can be massive and our storage capacity on NeSI is limited.

The final output (`bam`) files can be used as input for downstream processeses in our command line work, or can be imported into GUI programs for a visual analysis. Although we rely on NeSI for performing mapping in these exercises, we will use `Geneious` for visualising the results.

There are several ways to import `bam` files into `Geneious` either by using the `Import Files...` option from the `File` menu, or simply dragging and dropping the files into the `Geneious` document table. In either case, it is important that **both** the `bam` file and the reference sequence against which the reads were mapped are imported. If you import the `bam` without the reference, `Geneious` will not know the content of the original reference sequence.

---

## Different tools for the job

When it comes to mapping reads to a reference sequence, there is a vast list of tools available with different strengths and weaknesses. As this is a constantly expanding area of the literature we are going to focus on two mapping tools for these exercises but for your specific application they may not be the appropriate choice. If you are in a position where you need to perform reference mapping there are two main questions to consider when deciding which mapping tool to use.

1. Are you working with short (Illumina) or long (Nanopore, PacBio) sequencing data?
   1. The heuristics used for accelerating the mapping process differ between short and long sequences, due to both the differences in length and expected average read quality between the platforms.
1. Are you mapping DNA or RNA sequence data to your reference?
   1. Most of the time your reference will be a genome sequence.
   1. If you are working with RNA sequence data from eukaryotic (or some viral) species then some portion of reads will most likely have undergone splicing after transcription.
   1. The reference genome contains the full gene sequences not just the exonic regions so a mapping tool which is aware of gene splicing and can map around a spliced intron sequences will achieve a higher mapping rate than a mapping tool which only considers direct matches to the reference sequence.

As a starting point, some of the most popular mapping tools and their uses include:

* [minimap2](https://lh3.github.io/minimap2/) ([Li, 2018](https://doi.org/10.1093/bioinformatics/bty191))
  * Short and long read mapping
  * Specific mapping models for Nanopore and PacBio sequence profiles
  * Option for splice-aware mapping
* [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) ([Landmead & Salzberg, 2012](https://doi.org/10.1038/nmeth.1923))
  * Short read mapping
  * DNA mapping applications
* [HISAT2](http://daehwankimlab.github.io/hisat2/manual/) ([Kim *et al.*, 2019](https://doi.org/10.1038/s41587-019-0201-4))
  * Short read mapping
  * RNA (splice-aware) mapping, haplotype resolution
* [bwa-mem](https://bio-bwa.sourceforge.net/bwa.shtml) ([Li, 2013](https://arxiv.org/abs/1303.3997))
  * Short read mapping
  * DNA mapping applications
  * Long read mapping settings exist, but have been superceded by `minimap2`

---

## Obtaining a reference sequence

In order to perform read mapping we need to have a downloaded copy of the sequence we wish to map against. Generally this will be either a completed or near-complete reference genome of the species/strain of interest, or a closely related organism. In this exercise, however, we are going to use a single gene so that when we later visualise the results in `Geneious` the display will be easy to navigate.

Performing the reference download will not be performed in this section and you have been provided with a reference sequence for this exercise.

Navigate to the `/nesi/project/nesi03181/phel/USERNAME/level2/mapping/` folder and check the contents of the `references/` folder to begin.


```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/level2/mapping/
$ ls references/
```

There are two sequence files in this directory. The larger of the two is a copy of the representative genome for *Mycoplasmopsis bovis* 8790 (accession NZ_LAUS00000000). Finding the correct reference sequence is less bioinformatics and more a matter of subject expertise so we are not covering the details of how to perform this work. The commands used to produce these reference sequences are provided below in case this is something you wish to repeat in your own work.

> <details>
> <summary>Downloading the reference sequence</summary>
>
> If you are interested in how this sequence was obtained, it was downloaded from the NCBI GenBank database using the `Entrez Direct` toolkit provided by NCBI ([https://www.ncbi.nlm.nih.gov/books/NBK179288/](https://www.ncbi.nlm.nih.gov/books/NBK179288/)). Once the genome was obtained the coordinates of the 16S ribosomal RNA sequence were found by browsing the annotation information for the genome ([available here](https://www.ncbi.nlm.nih.gov/nuccore/1632735866?report=genbank)) by simply searching the web page for the text '16S ribosomal RNA'.
>
> The coordinate information for this sequence if given by the line entry `rRNA 41621..43134` which we can use in combination with `seqmagick` to extract just the region of interest for our mapping example.
>
> ```
> $ module load entrez-direct/13.3 seqmagick/0.8.4-gimkl-2020a-Python-3.8.2
> $
> $ efetch -format fasta -db sequences -id NZ_LAUS01000004 > Mbovis_87900.genome.fna
> $ seqmagick convert --cut 41621-43134 Mbovis_87900.genome.fna Mbovis_87900.16S_rRNA.fna
> ```
>
> </details>

---
