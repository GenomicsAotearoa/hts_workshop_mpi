# 3.2 - Mapping Illumina sequences to a reference

!!! clock "time"

    * Teaching: 10 minutes
    * Exercises: 20 minutes

!!! circle-info "Objectives and Key points"

    #### Objectives
    
    * Use `bowtie2` to index a reference genome and map DNA paired-end reads against a reference genome.
    
    #### Keypoints
    
    * Understand how to index a reference sequence for mapping.
    * Understand how to apply `bowtie2` to map a set of DNA paired-end reads to the reference.

---

## Indexing the reference sequence

Before we can map our sequence data to the reference genome (or gene sequence) obtained in the [previous exercise](./31_coverage_mapping.md#obtaining-a-reference-sequence) we need to perform a step known as *indexing*. How this process works is well beyond the scope of this tutorial, but it is a process of performing a scan of the reference sequence and transforming it into an organised format amenable to the `bowtie2` rapid mapping algorithm.

To perform mapping, navigate to the `/nesi/project/nesi03181/phel/USERNAME/level2/mapping/` directory and perform the following commands:
 
```bash
$ module load Bowtie2/2.4.5-GCC-11.3.0
$ bowtie2-build references/Mbovis_87900.16S_rRNA.fna references/Mbovis_87900.16S_rRNA
```

This will create a series of files with the prefix `Mbovis_87900.16S_rRNA` and extensions `.bt2` in the `references/` folder. The number of files depends on the size of the reference sequence which was indexed, but collectively these files comprise the index for the reference. When we perform mapping we specify the target as the file path `reference/Mbovis_87900.16S_rRNA` and `bowtie2` will automatically find and make sense of the index files.

Once this has completed, also index the full genome file (`Mbovis_87900.genome.fna`) to practice mapping against a more realistic reference sequence.

---

## Mapping reads with `bowtie2`

Once you have an index produced, it is now time to align (map) the short sequences against the reference. The nature of the `bowtie2` mapping tool is that it can be run with one of several preset configurations depending on your requirements, or you can devote quite a bit of time ot fine-tuning the parameters to optimise your output. Have a quick look through the mapping options by loading the `bowtie2` module and running the help command. This can be done with a single command although we need to provide several pieces of information to the tool:

```bash
$ module load Bowtie2/2.4.5-GCC-11.3.0
$ bowtie2 -h
```

There's a lot to take in here, but really the most important parts to take away from the command are as follows:

|Parameter|Value|Purpose|
|:---:|:---:|:---|
|`--sensitive`||Use the 'sensitive' mapping parameters for end-to-end read mapping.<br />Mapping can be performed on a sliding scale changing sensitivity (thoroughness) for speed.|
|`-x`|`references/Mbovis_87900.16S_rRNA`|The path to the index file(s) corresponding to the reference sequence|
|`-1`|`reads/Mbovis_87900.miseq_R1.fq.gz`|The forward reads file to be mapped|
|`-2`|`reads/Mbovis_87900.miseq_R2.fq.gz`|The reverse reads file to be mapped|
|`-U`||The singleton reads file to be mapped, if necessary.<br />*We are not using singleton reads in this exercise.*|
|`-S`|`Mbovis_87900.16S_rRNA.bowtie2.sam`|The output `sam` file to which the results are to be written|

Assemble these values into a command and run:

```bash
$ bowtie2 --sensitive \
      -x references/Mbovis_87900.16S_rRNA \
      -1 reads/Mbovis_87900.miseq_R1.fq.gz \
      -2 reads/Mbovis_87900.miseq_R2.fq.gz \
      -S Mbovis_87900.16S_rRNA.bowtie2.sam
```

> ### Exercise
>
> Now perform mapping against your full reference genome. For this exercise, you will need to increase the number of computing threads to the maximum available in your `JupyterHub` session to get a reasonable run time. Use the help manual to find the parameter for adjusting the thread number.
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> $ module load Bowtie2/2.4.5-GCC-11.3.0
>
> # Where 'N' is the number of threads available in your JupyterHub session
> $ bowtie2 --sensitive \
>       --threads N \
>       -x references/Mbovis_87900.genome \
>       -1 reads/Mbovis_87900.miseq_R1.fq.gz \
>       -2 reads/Mbovis_87900.miseq_R2.fq.gz \
>       -S Mbovis_87900.genome.bowtie2.sam
> ```
>
> </details>

---
