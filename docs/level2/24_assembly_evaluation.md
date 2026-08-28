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

     Create a new directory and copy in the `assembly_evaluation/` folder and copy in your `SPAdes` and `raven` fasta and assembly graphs.

    ??? circle-check "Solution"
 
        !!! terminal "code"
        
            ```bash
            mkdir assemblies/

            cp ../assembly_illumina/assembly/contigs.fasta assemblies/spades.fna

            cp ../assembly_nanopore/raven.fna assemblies/
            ```

??? failure "Help, my assembly failed!"

    If your assembly did not complete, don't worry about it. There is a training set of assemblies we can provide for you if required.

Once you have a local copy of your assemblies, we will be comparing these to the reference genome using a tool called [QUAST](https://github.com/ablab/quast).

Running `QUAST` is quite simple:

!!! terminal "code"

    ```bash
    quast.py -r reference/Mbovis_87900.genome.fna --gene-finding -o quast/ assemblies/*.fna
    ```

??? success "Output"

    ```
    Version: 5.3.0

    System information:
    OS: Linux-5.14.0-503.40.1.el9_5.x86_64-x86_64-with-glibc2.35 (linux_64)
    Python version: 3.12.14
    CPUs number: 32

    # Text omitted...

    Elapsed time: 0:00:05.509653
    NOTICEs: 5; WARNINGs: 3; non-fatal ERRORs: 0

    Thank you for using QUAST!
    ```

Open the resulting `quast/report.pdf` file in Jupyter using the file browser, or download the `quast/report.html` file to view it locally. Take a look through the report and see if you can get a feel for how well your assemblies compare to the reference.

How do the Illumina and Nanopore assemblies differ, if at all?

---

## (Optional) Visualising assemblies with `Bandage`

We can also visualise the assemblies by looking at how well the loops and fragments of the assembly graph were resolved. For this, we require a different set of files from the assembly output folders.

!!! question "Exercise"

    Copy the `.fastg` (`SPAdes`) and `.gfa` (`raven`) files from your previous output folders into your current assembly directory, ready for analysis.

    ??? circle-check "Solution"
 
        !!! terminal "code"
        
            ```bash
            cp ../assembly_illumina/assembly/assembly_graph.fastg assemblies/spades.fastg

            cp ../assembly_nanopore/raven.gfa assemblies/
            ```

Running the tool is then a matter of:

!!! terminal "code"

    ```bash
    Bandage image assemblies/spades.fastg spades_bandage.svg
    ```

You can then open the `spades_bandage.svg` file in the Jupyter browser. Unfortunately, we cannot filter out the short contigs from this result. However, it should be clear that there is one long contig which has been assembled, and then a large number of short fragments.

!!! question "Exercise"

    Repeat the `Bandage` command for your Nanopore assembly, then contrast the result from what you obtained with `SPAdes`. How do they differ?

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            Bandage image assemblies/raven.gfa raven_bandage.svg
            ```

        The data have assembled cleanly into a single contig, without bubbles, and there are no short fragments plotted.

!!! info "Assemblies are not always this clean!"

    Although there are slight differences in the quality of these two assemblies, they are still extremely good attempts at recovering the *M. bovis* genome. This will not always be the case! With a different set of input data, significantly worse results can be seen.

    If you look inside the `examples/` directory, there are two examples of alternate forms of this genome produced from slightly different sets of the same input data.

---

## Concluding comments

As you can see from this exercise, getting a *pretty good* genome assembly is not particularly difficult with the right tools. However, the distance between a draft assembly, which we have produced, and a final completed genome is a very long process and involved multiple rounds of assembly refinement, scaffolding, and often requires the creation of custom primers to perform PCRs to close sequence gaps which were not covered in your HTS library.

It can be hard knowing when the assembly is good enough to move out of the assembly stage and into annotation. In research groups working with genomic data, the yardstick for working with these kinds of data is typically to ask whether the current assembly is sufficient to answer the research question which led to its sequencing in the first place.

We can copy this logic and ask, what was the purpose of sequencing this genome and can we achieve that with the current data. Typically, we are most likely looking to perform a species identification. If we find that the genome assembly contains the right marker genes or operons to perform the identification then, regardless of whether the genome is officially completed or not, it has served its purpose.

---
