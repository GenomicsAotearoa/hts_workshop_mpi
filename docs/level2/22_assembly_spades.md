# 2.2 - Assembling short-read Illumina data

## Overview

!!! clock "time"

    * Teaching: 10 minutes
    * Exercises: 10 minutes
    
!!! circle-info "Objectives and Key points"

    #### Objectives
    
    * Perform an assemble of a bacterial genome using the `SPAdes` assembly tool.
    
    #### Keypoints
    
    * The `SPAdes` genome assembler is a powerful tool for assemnling genoms and contigs from a wide range of sampel types.

---

## Introduction to the SPAdes assembly tool

The `SPAdes` assembler is a very powerful and popular assembler which utilises de Bruijn graphs to assemble short read sequence data into larger contigs.

!!! note "SPAdes is particularly good at..."

    * *De novo* assembly of short read sequences (i.e. Illumina or IonTorrent)
    * Assembling small genomes (ideally <100 Mb, such as bacterial, viral, fungal, mitochondrial genomes)

    Despite this, it is possible to apply `SPAdes` to larger genomes and obtain very good results.

In order to begin, we must first find the versions of `SPAdes` installed on NeSI and load the module of interest.

```bash
$ module load SPAdes/3.15.2-gimkl-2020a

$ spades.py -h
```

!!! warning "Remember to filter your data prior to assembly!"

    Before beginning to work with `SPAdes` we need to ensure that our data is free from adapter sequences. The main reason for doing this is when we assemble sequences to form a genome, the assembler is looking for spans of nucleic acid sequence which are common to multiple reads, so that those reads can be joined together to create longer contigs. As the adapter sequence is an identical tag added to every read, these create regions of artificial homology between sequences which have no real connection to each other.

    Before we attempt assembly it is critical to remove these from our data but as this was covered in the level 1 training we will not be revisiting the process here.

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
    #SBATCH --time          00:40:00
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
    spades.py --isolate --threads ${SLURM_CPUS_PER_TASK} \
        -1 reads/Mbovis_87900.miseq_R1.fq.gz \
        -2 reads/Mbovis_87900.miseq_R2.fq.gz \
        -o assembly/
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
1. `spades.log` - the log file of the steps `SPAdes` performed and any warnings which occurred during assembly.
1. `assembly_graph.fastg` - a map of how well the assembly is resolved.
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

It will take a while for `SPAdes` to complete, so we will move to the next exercise and return to this folder later.

---
