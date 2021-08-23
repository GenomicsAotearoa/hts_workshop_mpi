# Level 3.2 revision

This revision material consists of three short exercises covering the main points of the second part of Level 3 training. For each of these exercises you will need to work from your `3_Assembly_mapping/` directory. **_However_** - keep a close eye on the names of your own folders during these exercises. Each person created their folders and therefore there are several variations in names that we have seen in the workshops, such as;

* `3_Assembly_mapping/`
* `3_Assembly-mapping/`
* `3_assembly_mapping/`

It does not matter which you use, but make sure that your commands are specific to your own folders and not those of the trainers. If you find that your `3_Assembly_mapping/` folder has become fragmented across different names, it might be worth consolidating your data into a single location before proceeding.

When you have completed all exercises, email the requested material and answers to the trainers.

> Remember, if you are struggling with any exercises:
>
> 1. Every exercise here was covered in the training material. Refer to this for hints.
> 1. It is not cheating to use Google to find help on how various commands work.
> 1. You are welcome to contact the trainers for advice on a particular exercise, but please attempt the first two options before resorting to this.

---

## Exercise 1 - Genome assembly with SPAdes

In the `/nesi/project/nesi03181/phel/module_3/2_Quality_filtered_data/` folder you will find quality filtered Illumina MiSeq reads for six different *M. bovis* isolates, different from the ones that we have used in the workshops so far:

1. Mb166
1. Mb182
1. Mb183
1. Mb194
1. Mb240
1. Mb267

These files ahve already been processed with `fastp` to remove low quality sequences and adapter sequence, producing familiar `Mb###_1.trim.fq.gz`, `Mb###_2.trim.fq.gz`, and `Mb###_s.trim.fq.gz` files. Select any one of these isolates, and write a `slurm` script to performing a genome assembly using `SPAdes`.

For this assembly, you must set the following parameters:

1. Assembly will take place in **isolate** mode.
1. Manually specify *k*-mer sizes for assembly. It is up to you which values to choose, but you must make the selection rather than let `SPAdes` pick automatically.
1. Set the number of threads to 20.
   1. *Note - you will have to account for this in your `slurm` batch file as well.*
1. Write the output files to your equivalent of the `3_Assembly-mapping/` folder.

You must also provide an appropriate amount of memory and run time for the command to finish. For most samples, selecting 10 - 20 GB of RAM and 2 hours of run time should be sufficient.

When the assembly is complete, copy and rename the `contigs.fasta` and `assembly_graph.fastg` files. These will be used in later exercises.

---

## Exercise 2 - Genome assembly with Flye

In the sample folder you will find a match set of quality filtered MinION sequences from the same six genomes, named as `Mb###.nanofilt.fq.gz`. Perform an assembly on the same isolate using `Flye`, or any other long-read assembly tool that you wish to trial.

For any assembly tool that requires an estimate of genome size, 1 Mb (1,000,000 bases) is sufficient for the *M. bovis* genome.

Write the output to your equivalent of the `3_Assembly-mapping/` folder. When assembly is complete, copy and rename the resulting contigs and assembly graph files. These will be named `assembly.fasta` and `assembly_graph.gfa` if you use `Flye`.

---

## Exercise 3 - Performing a hybrid assembly

Modify your `SPAdes` script to incorporate the MinION sequences into the assembly process. Submit the assembly job, then copy and rename the output contigs and assembly graph as for the previous exercises.

>**Note:** For this exercise you only need to add a single parameter to your existing script, and change the output folder name. All other parameters can be left the same, and the run time and memory limits used in the original assembly should be sufficient to complete the hybrid assembly.

---

## Exercise 4 - Evaluating the differences in assembly quality

For each of the three assembly graphs you have produced (short read, long read, hybrid) summarise the assembly outputs using `QUAST` with the *M. bovis* reference genome. Inspect the final output and identify features that favour each of the three methods.

You may also wish to produce `Bandage` plots of the assemblies for aid in your conclusions.

Email your report file and key findings to the trainers.

---
