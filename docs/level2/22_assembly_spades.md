# 2.2 - Assembling short-read Illumina data

## Overview

!!! clock "time"

    * Teaching: 15 minutes
    * Exercises: 45 minutes
    
!!! circle-info "Objectives and Key points"

    #### Objectives
    
    * Perform an assemble of a bacterial genome using the `SPAdes` assembly tool.
    * Use `QUAST` to assess the assembly status.
    
    #### Keypoints
    
    * ...

---

## Introduction to the SPAdes assembly tool

The `SPAdes` assembler is a very powerful and popular assembler which utilises de Bruijn graphs to assemble short read sequence data into larger contigs.

SPAdes is particularly good at:

* *De novo* assembly of short read sequences (i.e. Illumina or IonTorrent)
* Assembling small genomes (ideally <100 Mb, such as bacterial, viral, fungal, mitochondrial genomes)
  * Despite this, it is possible to apply `SPAdes` to larger genomes and obtain very good results.

Before beginning to work with `SPAdes` we need to ensure that our data is free from adapter sequences. The main reason for doing this is when we assemble sequences to form a genome, the assembler is looking for spans of nucleic acid sequence which are common to multiple reads, so that those reads can be joined together to create longer contigs. As the adapater sequence is an identical tag added to every read, these create regions of artificial homology between sequences which have no real connection to each other. Before we attempt assembly it is critical to remove these from our data.

In order to begin, we must first find the versions of `SPAdes` installed on NeSI and load the module of interest.

```bash
$ module load SPAdes/3.15.2-gimkl-2020a

$ spades.py -h
```

---

## Performing an assembly using SPAdes

The test data is a set of Illumina MiSeq sequencing reads from an *M. bovis* genome. To save time, the reads have already been quality filtered with `fastp` and the number of reads reduced to speed up analysis.

Navigate to the `assembly_illumina/` folder to begin.

!!! terminal "code"

    ```bash
    cd /nesi/project/nesi03181/phel/USERNAME/level2/assembly_illumina/
    ```

!!! question "Exercise"

    Use your knowledge of the module system on NeSI to find the most recent release of `SPAdes`.

    ??? circle-check "Solution"
 
        !!! terminal "code"
        
            ```bash
            module avail spades

            # OR...
            module spider spades
            ```

        ??? success "Output"

            ```
            ----------------------------------------------- /opt/nesi/CS400_centos7_bdw/modules/all -----------------------------------------------
              SPAdes/3.13.1-gimkl-2018b    SPAdes/3.15.0-gimkl-2020a    SPAdes/3.15.3-gimkl-2020a
              SPAdes/3.14.0-gimkl-2020a    SPAdes/3.15.2-gimkl-2020a    SPAdes/3.15.4-gimkl-2022a-Python-3.10.5 (D)
            ```

Create a `slurm` script with the following contents. Be sure to replace the `YOUR_EMAIL` and `USERNAME` values with your details.

!!! file-code "spades_asm.sl"

    ```bash
    #!/bin/bash -e
    #SBATCH --account       nesi03181
    #SBATCH --job-name      spades_asm
    #SBATCH --time          01:00:00
    #SBATCH --cpus-per-task 16
    #SBATCH --mem           20G
    #SBATCH --error         spades_asm.%j.err
    #SBATCH --output        spades_asm.%j.out
    #SBATCH --mail-type     END
    #SBATCH --mail-user     YOUR_EMAIL

    module purge
    module load SPAdes/3.15.4-gimkl-2022a-Python-3.10.5

    # Set the path from which the script will execute SPAdes
    cd /nesi/project/nesi03181/phel/USERNAME/assembly_illumina/

    # Execute SPAdes
    spades.py --isolate --threads 16 -1 reads/miseq_R1.qc.fq.gz -2 reads/miseq_R2.qc.fq.gz -o assembly/
    ```

Submit this job to `slurm`:

!!! terminal "code"

    ```bash
    sbatch spades_asm.sl
    ```

??? success "Output"

    ```
    Submitted batch job ########
    ```

When this job is complete, we will have a folder named `assembly/` which contains a lot of different files and pieces of information. A lot of the contents are files generated internally by `SPAdes` and we do not need to pay attention to them. The most important files for us are:

1. `contigs.fasta` - the assembled contigs, as a fasta file.
1. `scaffolds.fasta` - the scaffolded contigs, as a fasta file.
   1. This is similar to the contigs file, except it will contain sequences where contigs have been joined by an indeterminate piece of sequence.
   1. This occurs when `SPAdes` can tell from the pairing information in our library that two contigs belong adjacent to each other, but it has no information for filling the gap between them.
1. `spades.log` - the log file of the steps `SPAdes` performed and any warnings which occured during assembly.
1. `assembly_graph.fastg` - a map of how well the assembyl is resolved.
   1. This can be useful if we are trying to obtain a complete genome, as it shows us areas which have assembled cleaning and areas which were difficult to resolve.

!!! info "Different modes for running `SPAdes`"

    If you checked the help command for `SPAdes`, or were examining the contents of the `slurm` script carefully, you might have noted the flag `--isolate` which we passed to the command.

    `SPAdes` was originally created as a tool for a very specific use case - the assembly of single-cell amplified genomes, which were problematic to assembly due to the highly uneven coverage of the nucleic acid content in the input libraries ([Bankevich *et al.*, 2012](https://doi.org/10.1089/cmb.2012.0021)). It was also noted at the time that despite this fairly niche scope for use, the tool outperformed many traditional short read assemblers of the time and so over the years the team behind `SPAdes` have expanded the tools and added a number of alternate assembly modes into successive iterations of the tool.

    There are now specific subroutines in the `SPAdes` assembler for working with the following data types:

    1. metaSPAdes - assembly of metagenomic samples (sequence libraries from multiple organisms)
    1. plasmidSPAdes - selective assembly of plasmids from genomic libraries
    1. metaplasmidSPAdes - selective assembly of plasmids from metagenomic 
    1. rnaSPAdes - *de novo* assembly of transcriptome libraries
    1. biosyntheticSPAdes - selective assembly of biosynthetic gene clusters (BCG) from a library
    1. rnaviralSPAdes - selective *de novo* assembly of RNA viruses from transcriptome, metatranscriptome, and metavirome libraries

---

## Assessing the results of the assembly

Once assembly is complete, we have a complex process of determining the quality of the assembly. How 'good' a genome is can be difficult to measure, but as we are mostly working with well characterised pathogens a good starting place is to compare our assembled genome with previously characterised members of the same species to see how well the conserved genomic features have been reconstructed by our assembly tool.

You have been provided with a copy of the *Mycoplasmopsis bovis* PG45 genome in your `assembly_illumina/` folder, and this is the reference we will be comparing against, using a tool called [QUAST](https://github.com/ablab/quast).

Running `QUAST` is quite simple:

!!! terminal "code"

    ```bash
    module load QUAST/5.2.0-gimkl-2022a

    quast.py -r reference/GCF_000183385.1.fna --gene-finding -o quast/ assembly/contigs.fasta
    ```

??? success "Output"

    ```
    Version: 5.2.0

    System information:
      OS: Linux-3.10.0-693.2.2.el7.x86_64-x86_64-with-glibc2.17 (linux_64)
      Python version: 3.10.5
      CPUs number: 2

    Started: 2023-09-28 15:02:30

    # Text omitted...

    Finished: 2023-09-28 15:02:40
    Elapsed time: 0:00:10.422580
    NOTICEs: 4; WARNINGs: 1; non-fatal ERRORs: 0

    Thank you for using QUAST!
    ```

````

Open the resulting `quast/report.pdf` file in Jupyter using the file browser. How complete does the assembly appear to be, compared with the reference?

??? info "Visualising assemblies with `bandage`"

    We can also visualise the assembly by looking at how well the loops and fragments of the assembly graph were resolved. For this, we require only the `assembly_graph.fastg` file from the `SPAdes` output:

    !!! terminal "code"

    ```bash
    module load Bandage/0.8.1_Centos

    Bandage image assembly/assembly_graph.fastg assembly_bandage.svg
    ```

    You can then open the `assembly_bandage.svg` file in the Jupyter browser. Unfortunately, we cannot filter out the short contigs from this result. However, it should be clear that there is one long contig which has been assembled, and then a large number of short fragments.

---

## Concluding comments

...

---
