# 2.3 - Evaluating the results of an assembly

## Overview

!!! clock "time"

    * Teaching: 20 minutes
    * Exercises: 15 minutes
    
!!! circle-info "Objectives and Key points"

    #### Objectives

    * Use `QUAST` to assess the assembly status.
    * (Optional) Use `Bandage` to view how well the assembly resolved.

    #### Keypoints

    * Tools like `QUAST` can be used to perform quick and easy comparisons between an assembly and a trusted reference genome.
    * It is important to make sure that your genome is sufficiently resolved to address your need, but we often do not need to go further than what an assembler provides.

---

## Assessing the results of an assembly

Once assembly is complete, we have a complex process of determining the quality of the assembly. How 'good' a genome is can be difficult to measure, but as we are mostly working with well characterised pathogens a good starting place is to compare our assembled genome with previously characterised members of the same species to see how well the conserved genomic features have been reconstructed by our assembly tool.

Navigate to the `assembly_evaluation/` folder and we will begin.

You have been provided with a copy of a reference *Mycoplasmopsis bovis* genome in the `reference/` folder, but we will need some draft assemblies to test as part of this module.

!!! question "Exercise"

     Create a new directory and copy in the `assembly_evaluation/` folder and copy in your `SPAdes` and `Flye` fasta and fastg assembly files.

    ??? circle-check "Solution"
 
        !!! terminal "code"
        
            ```bash
            mkdir assemblies/

            cp ../assembly_illumina/assembly/contigs.fasta assemblies/

            cp ../assembly_nanopore/assembly/assembly.fasta assemblies/
            ```

??? failure "Help, my assembly failed!"

    If your assembly did not complete, don't worry about it. There is a training set of assemblies we can provide for you if required.

Once you have a local copy of your assemblies, we will be comparing these to the reference genome using a tool called [QUAST](https://github.com/ablab/quast).

Running `QUAST` is quite simple:

!!! terminal "code"

    ```bash
    module load QUAST/5.2.0-gimkl-2022a

    quast.py -r reference/Mbovis_87900.genome.fna --gene-finding -o quast/ assemblies/*.fasta
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

Open the resulting `quast/report.pdf` file in Jupyter using the file browser. Take a look through the report and see if you can get a feel for how well your assemblies compare to the reference.

How do the Illumina and Nanopore assemblies differ, if at all?

---

## (Optional) Visualising assemblies with `Bandage`

We can also visualise the assemblies by looking at how well the loops and fragments of the assembly graph were resolved. For this, we require a different set of files from the assembly output folders.

!!! question "Exercise"

    Copy the `.fastg` (`SPAdes`) and `*.gfa` (`Flye`) files from your previous output folders into your current assembly directory, ready for analysis.

    ??? circle-check "Solution"
 
        !!! terminal "code"
        
            ```bash
            cp ../assembly_illumina/assembly/assembly_graph.fastg assemblies/
            cp ../assembly_nanopore/assembly/assembly_graph.gfa assemblies/
            ```

Running the tool is then a matter of:

!!! terminal "code"

    ```bash
    module load Bandage/0.8.1_Centos

    Bandage image assemblies/assembly_graph.fastg spades_bandage.svg
    ```

You can then open the `assembly_bandage.svg` file in the Jupyter browser. Unfortunately, we cannot filter out the short contigs from this result. However, it should be clear that there is one long contig which has been assembled, and then a large number of short fragments.

!!! question "Exercise"

    Repeat the `Bandage` command for your Nanopore assembly, then contrast the result from what you obtained with `SPAdes`. How do they differ?

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            Bandage image assemblies/assembly_graph.gfa flye_bandage.sv
            ```

        The data have assembled cleanly into a single contig, without bubbles, and there are no short fragments plotted.

---

## Concluding comments

As you can see from this exercise, getting a *pretty good* genome assembly is not particularly difficult with the right tools. However, the distance between a draft assembly, which we have produced, and a final completed genome is a very long process and involved multiple rounds of assembly refinement, scaffolding, and often requires the creation of custom primers to perform PCRs to close sequence gaps which were not covered in your HTS library.

It can be hard knowing when the assembly is good enough to move out of the assembly stage and into annotation. In research groups working with genomic data, the yardstick for working with these kinds of data is typically to ask whether the current assembly is sufficient to answer the research question which led to its sequencing in the first place.

We can copy this logic and ask, what was the purpose of sequencing this genome and can we achieve that with the current data. Typically, we are most likely looking to perform a species identification. If we find that the genome assembly contains the right marker genes or operons to perform the identification then, regardless of whether the genome is officially completed or not, it has served its purpose.

---
