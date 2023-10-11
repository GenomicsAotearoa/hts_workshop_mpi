# 3.3 - Mapping Oxford Nanopore sequences to a reference

!!! clock "time"

    * Teaching: 10 minutes
    * Exercises: 5 minutes

!!! circle-info "Objectives and Key points"

    #### Objectives
    
    * Use `minimap2` to index a reference genome and map long-read data produced with Oxford Nanopore sequencing technology to a reference genome.
    
    #### Keypoints
    
    * Understand how to index a reference sequence for mapping with `minimap2`, and know when you may need to perform this step.
    * Understand how to apply `minimap2` to map a set of DNA paired-end reads to the reference.

---

## Indexing the reference sequence

Similar to `bowtie2`, it is possible to create a pre-computed index file for performing mapping using `minimap2`. However this is not strictly necessary - `minimap2` will automatically read the reference file and determine if it is indexed or not. If the file is in fasta format, `minimap2` will produce an index file on the fly for performing the alignment operations. If the index file is already indexed, this step is skipped.

??? question "When should we index the refernece file?"

    Whether or not you need to manually index your reference depends on the situation. The indexing process in `minimap2` is very fast and for quick, one-off applications you can probably skip it. If you are producing an index which is going to be reused many times it might be worth creating a common index file which can be recycled between mapping applications. This might particularly be desirable if you application is for some routine diagnostic purpose, so that you have a single mapping index for all analyses for quality tradking purposes.

To perform mapping, navigate to the `/nesi/project/nesi03181/phel/USERNAME/level2/mapping/` directory again, and map the smaller fasta file using `minimap2`:

!!! terminal "code"

    ```bash
    module purge
    module load minimap2/2.24-GCC-11.3.0

    minimap2 -d references/Mbovis_87900.genome.mmi references/Mbovis_87900.genome.fna
    ```

Note that this time we provide the target path for the index file as a named parameter *before* providing the path to the fasta file to be indexed. We also have specified an extension for the output file (`.mmi`). Unlike `bowtie2`, which splits the indexing information over a number of files as size dictates, `minimap2` contains all indexing information in a single file. This can be handy for portability, and can also be convenient when writing commands as it is easier to use tab-completion to get to the index file directly.

---

## Mapping reads with `minimap2`

Similar to `bowtie2` there are a number of pre-configured settings for mapping with `minimap2`.

For most cases, the most important parameter to be aware of is the `-x` parameter. This is the toggle for applying appropriate parameters for applying short- or long-read mapping, with mapping profiles for PacBio and Nanopore sequence reads. The data we are working with in this exercise is based on the Oxford Nanopore MinION technology so this is the mapping preset we will apply.

!!! terminal "code"

    ```bash
    minimap2 -ax map-ont references/Mbovis_87900.genome.mmi reads/Mbovis_87900.nanopore.fq.gz > Mbovis_87900.genome.nanopore.sam
    ```

---

## Differences in sequence quality between short- and long-read platforms

The nature of long-read sequencing is that it is inherently more error-prone than short-read sequencing. This is due to the fundamental difference in approach between long-read platforms, which sequence individual nucleic acid sequences, and short-read platforms which perform PCR amplification to create a clonal population of sequences from which sequencing signals are generated.

While the error rates of Nanopore sequences are drastically better than they were a few years ago they are still more prone to spontaneous sequencing errors when compared with Illumina sequences. In practice, this error can be as little as 1% difference between platforms but when Nanopore sequencing is producing sequences which are thousands of nucleotides in length, a 1% error rate does result in significant divergence from the original sequence.

In addition to stochastic error, which all platforms suffer from, there are also instances of platform-specific errors in which particular nucleotide motifs exhibit a propensity to certain error types. Dedicated long-read mapping tools, including `minimap2`, are trained on these platform-specific error profiles and can achieve greater mapping accuracy by according for these errors when they are aware of the sequencing platform used to produce the input sequences.

It is therefore critical to note which sequencing method was used in producing your HTS data and perform reference mapping with parameters customised to the platform.

---
