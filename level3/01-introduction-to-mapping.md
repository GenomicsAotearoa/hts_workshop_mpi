# Mapping sequences to a reference genome

* Teaching: 20 minutes

#### Objectives

* Know how to map reads to a reference sequence.
* Understand the output files from mapping tools.
* Learn how to import mapping results into GUI programs such as `Geneious` for further analysis.

#### Keypoints

* Reads can be mapped to any reference sequence, for example whole genomes or individual genes, to generate a consensus or view the coverage in particular regions.
* It is good practice to generate an index of your reference sequence prior to mapping - this will speed up the mapping process if you are performing multiple mapping operations.
* Mapping alignments are stored in the `sam` (text, human-readable) or `bam` (binary, compressed) formats.
* Different mapping tools work in different ways and make different assumptions when working with your data.
  * Sometimes, a given tool will work better with a given dataset.
  * It is recommended to try a few tools when working with new datasets to see which one produces better results.

---

## Contents

1. [Introduction to mapping](#introduction-to-mapping)
1. [Overview of the mapping process](#overview-of-the-mapping-process)
1. [Storing read mapping information](#storing-read-mapping-information)
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

![Mapping overview](../img/03_mapping.png)

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
|4|POS|1- based leftmost mapping position|
|5|MAPQ|Mapping quality|
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

The final output (`bam`) files can be used as input for downstream processeses in our command line work, or can be imported into GUI programs for a visual analysis. We will use `Geneious` for this purpose. There are several ways to import `bam` files into `Geneious` either by using the `Import Files...` option from the `File` menu, or simply dragging and dropping the files into the `Geneious` document table.

In either case, it is important that **both** the `bam` file and the reference sequence against which the reads were mapped are imported. If you import the `bam` without the reference, `Geneious` will not know the content of the original reference sequence.

---

## Obtaining a reference sequence

In the exercises that follow we will be working with a single gene as the reference sequence, rathre than the full genome. We have previously covered a method of obtaining full genome sequences using `wget` (see <a href="../level2/06-homework.md#exercise-3---writing-your-own-slurm-script" target="_blank">the module 2 homework</a>), but there is an easier way to obtain a reference if it is just a single sequence that we want.

If we know the accession number of a sequence to download, it can be retrieved using the `Entrez Direct` utilities ([https://www.ncbi.nlm.nih.gov/books/NBK179288/](https://www.ncbi.nlm.nih.gov/books/NBK179288/)) on NeSI quite simply. For example, to obtain the complete sequence for the *Mycoplasma bovis* CQ-W70 genome (accession CP005933) we can use the following command: 

```bash
$ module load entrez-direct/13.3

# Download the raw fasta file
$ efetch -format fasta -db sequences -id CP005933 > CP005933.fna

# Download the annotated GenBank file
$ efetch -format gb -db sequences -id CP005933 > CP005933.gb
```

To find the correct accession number, you will still need to perform your own digging to find a trusted reference. This is just a means to speed up the process of getting the data onto NeSI once you have chosen the correct sequence for use.

--

[Next lesson](02-illlumina-mapping.md)
