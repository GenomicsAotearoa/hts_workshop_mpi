# 3.2 - Filtering Nanopore data

!!! clock "time"

    * Teaching: 10 minutes
    * Exercises: 10 minutes

!!! circle-info "Learning objectives"

    **Objectives**

    * Be able to perform quality filtering of Nanopore data to remove short reads and adapter sequences

    **Key points**

    * Sequences can be filtered to remove low quality reads and short reads using tools such as `Nanofilt`
    * It is important to use tools designed specifically for long read data when performing assessment and trimming of Nanopore reads
    * Identifying and removing residual sequencing adapters is not always necessary, but tools exist for it if required

---

## Identifying and removing adapters

Before we start worring about the raw sequencing quality of our reads, there is one factor to consider. When creating HTS libraries for sequencing, barcode and adpater sequences are appended to the ends of the template DNA to facilitate sequencing. These are an artificial construct, not reflective of the read biology of our sample, but are sequenced the same as the template DNA by the sequence device.

Adapters may be found at the start of end of a sequence, as with Illumina reads, but an issue unique to Nanopore sequencing is the occurence of *mid-sequence* adapters. This occurs when a sequence moves through the sequencing pore and is immediately followed by another read. If this happens quickly enough that the base calling tool can realise the first sequence has ended, the the new sequence is appended to tail of the first.

When this phenomenon happens you end up with a chimeric fusion of the first and second sequence, separated by the adapter of the second sequence. While mid-sequence adapaters artifacts are rare (typically less than 1% of a sequencing library) they are very easy to screen for, and it is important to identify and correct them in order to get the best classification or assembly result possible.

Like quality filtering, adapater removal (including mid-sequence adapters) can be performed during base calling, but if it was not, we can use the tool `porechop` to identify and remove adapaters:

!!! terminal "code"

    ```bash
    module load Porechop/0.2.4-gimkl-2020a-Python-3.8.2
    porechop --threads 2 -i reads/nanopore.fq.gz -o results/nanopore.porechop.fq
    ```

    ??? success "Output"

        ```bash
        Trimming adapters from read ends
             SQK-NSK007_Y_Top: AATGTACTTCGTTCAGTTACGTATTGCT
          SQK-NSK007_Y_Bottom: GCAATACGTAACTGAACGAAGT
             1D2_part_2_start: CTTCGTTCAGTTACGTATTGCTGGCGTCTGCTT
               1D2_part_2_end: CACCCAAGCAGACGCCAGCAATACGTAACT

        3,722 / 3,722 (100.0%)

        3,120 / 3,722 reads had adapters trimmed from their start (68,630 bp removed)
          834 / 3,722 reads had adapters trimmed from their end (9,680 bp removed)


        Splitting reads containing middle adapters
        3,722 / 3,722 (100.0%)

        20 / 3,722 reads were split based on middle adapters
        ```

Watch the output as `porechop` runs, and see how many residual adapters were found in your mock data set. One of the really nice features of `porechop` is that we do not need to know which adapter sequences were used in library construction - the tool has a built-in reference set of the common adapter sequences and reads are screened for all by default.

!!! question "Exercise"

    From the output of `porechop`, answere the following questions:

    1. Which adapter sets were used to produce this library?
    1. How many sequences were flagged as carrying a forward or reverse adapter sequence?
    1. How many sequences had mid-sequence adapaters which needed removal?

    ??? circle-check "Solution"

        1. The adapters from the **SQK-NSK007** and **1D2** sequencing kits were found in this library. 
        1. 3,120 sequences were found to have a forward adapter attached, and 834 sequences had an end adapter attached.
        1. 20 of the reads were found ot have a mid-sequence adapter.

---

## Removing short and low quality reads

Finally, we need to remove low quality data from the sequencing library. It can also be helpful to remove short reads from the data set at this point, as these can negatively impact the results of your downstream analysis in some situations.For example, if you are trying to assemble a genome from the data, the short reads may interfere with the assembly process and result in a worse genome that you would get if they were screened rom the library beforehand.

`NanoFilt` is a tool that can filter reads by quality score and length. It is also capable of cropping a specified number of bases from the start or the end of a read, which can be useful sometimes as read quality is not uniform along the sequences - in Nanopore data it is common to see that the average quality score in the first few dozen nucleotides is significantly worse than the rest of the sequence.

Unload `porechop`, and load the `NanoFilt` module to:

* Remove all reads with quality scores under 15
* Remove all reads shorter than 500 bp
* Trim the first 50 nucleotides off all reads

!!! terminal "code"

    ```bash
    module purge
    module load nanofilt/2.6.0-gimkl-2020a-Python-3.8.2

    NanoFilt -q 15 -l 500 --headcrop 50 < results/nanopore.porechop.fq > results/nanopore.qc.fq
    ```

!!! note "Note"

    `NanoFilt` behaves a bit strangely in the way it takes input data and writes output. Rather than specifying files with flags such as `-i` and `-o`, it reads a data from a stream called `stdin`, and writes to a separate channel called `stdout`.

    Working with these channels is part of the Level 2 training, for now you just need to know that the `<` character sets the `stdin` input, and `>` specifies where the `stdout` data will go.

We can get a quick estimate for how much data passed our quality requirements by checking the number of lines in the output file. If you have read the [fastq description file](../supplementary/fastq_format.md) you will know how data is recorded in the fastq file format. We can quickly use the `wc` command to count the number of lines in the input and output files, to compare how many reads we started with and how many remain after filtering:

!!! terminal "code"

    ```bash
    wc -l results/*
    ```

    ??? success "Output"

        ```bash
        14940 results/nanopore.nanopore_porechop.fq
        5252 results/nanopore.qc.fq
        20192 total
        ```

A different more informative way of looking at our data before and after filtering is to use a tool called Seqkit. Seqkit has a range of features but the `stats` command gives really good summary statistics for sequencing datasets. Here we are using the flag `-a` to return all the available stats. 

!!! terminal "code"

    ```bash
    module purge
    module load SeqKit/2.4.0

    seqkit stats -a results/*.fq
    ```

    ??? success "Output"

        ```bash
        |file|format|type|num_seqs|sum_len|min_len|avg_len|max_len|Q1|Q2|Q3|sum_gap|N50|Q20(%)|Q30(%)|GC(%)|
        |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
        |results/nanopore.porechop.fq|FASTQ|DNA|3,735|45,214,457|6|12,105.6|82,753|2,723.5|8,489|17,927|0|21,118|69.79|50.97|30.83|
        |results/nanopore.qc.fq|FASTQ|DNA|1,313|16,675,742|500|12,700.5|77,866|3,658|9,245|18,869|0|21,365|78.07|58.66|31.03|
        ```

How many sequences did we lose in the process?

---
