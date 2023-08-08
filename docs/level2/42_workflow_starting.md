# Introduction to the Nextflow system

* Teaching: 60 minutes
* Exercises: 60 minutes

#### Objectives

* Build a small workflow in `Nextflow` to map sequences to a reference and produce a mapping summary.
* Know how to organise subunits of work into process blocks and how to link these using the `workflow{}` statement.

#### Keypoints

* A complex workflow is built from smaller sub-units of work, which are executed sequentially to produce a final outcome.
* The main sections of a `Nextflow` file are the `process{}` statements and the `workflow{}` block which governs the flow of data between the processes.

---

## Contents

1. [Getting started](#getting-started)
1. [Writing the outline for our workflow](#writing-the-outline-for-our-workflow)
1. [Executing our first workflow](#executing-our-first-workflow)
1. [Adding a second process to the workflow](#adding-a-second-process-to-the-workflow)
1. [Adding the final process and exporting the results](#adding-the-final-process-and-exporting-the-results)

---

## Getting started

For all of the work we do in this session, we require three modules to be loaded in NeSI. Load these now so that we do not need to think about them for the rest of the session.

```bash
$ module purge
$ module load minimap2/2.24-GCC-11.3.0
$ module load SAMtools/1.16.1-GCC-11.3.0
$ module load Nextflow/22.10.3
```

Once the modules are loaded navigate to the `/nesi/project/nesi03181/phel/USERNAME/level2/workflows/` directory for this exercise.

In this work, we are going to be using a 'learn by example' approach to teaching how a `Nextflow` workflow operates. If you wish for more technical information on how the interacting components are defined, interact, or how to extend upon the teachings of this tutorial and the next please consult the online `Nextflow` docmentation.

---

## Writing the outline for our workflow

A `Nextflow` workflow comprises two main features, **_processes_** and the **_workflow statement_**. Before we start to build our workflow, produce a small text file using your prefered text editor and enter the follow content into it. You can save this file under whatever name you wish - the traditional ending for a `Nextflow` rule is `.nf` but this is only a recommendation and not enforced by the software.

```
// Configuration
nextflow.enable.dsl=2

// Processes


// Workflow
workflow {

}
```

The lines which behin with `//` are comments - these are statements which are not read by the `Nextflow` engine and exist for our own benefit. Commenting in code is a useful way to break the lines up with notes to yourself or other users. You can use comments to explain what a particular region of a script is used for, or to explain the rationale for a particular implementation over another. In this instance we are using three comments to break our script up into three sections.

The first section is to do with the workflow configuration. We are inserting a single line here, which specifies to `Nextflow` that the workflow must be interpretted using the second release of its *domain-specific language*, which is the language style `Nextflow` uses to express the code blocks and relationships between them.

><details>
><summary>Methods of specifying the DSL version</summary>
>
> Although the second version of this language has been released for some time, the `Nextflow` engine still defaults to version 1 when interpretting commands so when producing a workflow we need to specify if we are deviating from this. There are three ways to tell `Nextflow` which version of the DSL we are using:
>
> 1. Add the `-dsl2` flag to the command line when executing the workflow.
> 1. Add the line `nextflow.enable.dsl=2` inside the workflow rule to overwrite the user-provided or default case for DSL interpretation.
> 1. Add the statement `nextflow.enable.dsl=2` to a configuration file which sits in the same directory as the `Nextflow` rule being executed (used here).
>    1. A configuration file must have the name `nextflow.config` and reside in the same directory as the `.nf` file to be recognised.
>    1. For large projects, or projects which can work across different server deployments this is the prefered way to address configuration settings but it is overly complex at this stage in our lesson.
></details>


The second comment is simply a placeholder for where we are going to write **_process_** statements. These are effectively the individual shell commands which we will write to perform steps of our analysis.

The third comment denotes the `workflow{}` statement. The **_workflow_** statement is the brains of the `Nextflow` rule. It assesses which input file(s) are to be analysed and organises the **_process_** statements to achieve the desired outcome.

---

## Writing the first process in the workflow

Now that we have the outline of our workflow we are going to take the `minimap2` reference mapping command from [the previous exercise on mapping](./33_nanopore_mapping.md) and transcribe it into `Nextflow`.

To begin, we are going to modify your basic document to include a single process, which maps a given input file against a pre-determined reference file.

>**Note:** If at any time you get lost, there are a set of completed `Nextflow` rules for each stage of this tutorial in a hidden folder named `.checkpoints/`

To begin, under your comment 'Processes', write the following commands:

```diff
// Processes
+process map_to_reference {

+   input:
+   path fq_file

+   output:
+   path "mapping.sam"

+   script:
+   """
+   minimap2 -ax map-ont ${launchDir}/references/Mbovis_87900.genome.mmi ${fq_file} > mapping.sam
+   """
+}
```

Each **_process_** is declared by the statement `process` followed by a unique name used to identify the process block. The process is defined as all file content between the opening and closing curly brace (`{`, and `}`). If you are familiar with programming in a C-like language you will recognise this as very similar to a function declaration.

Within the **_process_** block are three sections, `input`, `output` and `script`.

* `input:`
   * This block refers to the input data for the process.
    * This block provides two variables to us, the *value* `sample_id` and a list of *file paths* named `paired_reads`.
* `output:`
  * This block defines the name of any output files which must be produced upon *successful* completion of the **_process_**.
  * This *does not* have to capture every file produced as the **_process_** runs, these only the output file(s) which are to be made use of at the end of the **_process_**.
 * `script:`
   * This section is where we write the command to be executed.
   * By default this is read as the `bash` language you have been using in the tutorials leading up to this one so you should hopefully recognise the `${fasta_file}` style of accessing variables.
     * In the case of the `paired_reads` variable this is actually a list of values, where the first position corresponds to the forwards reads and the second position to the reverse reads. We access these using the square brackets in the command above (e.g. `${paired_reads[0]}`).
   * It is possible to replace this block with a different statement, `shell:` which you may see in some online examples.
     * This is mainly done when the code block needs to be written in a language other than `bash`, such as `perl`, `python`, or `R`.
     * Using the `shell:` statement changes the way in which variables are accessed so we will not be using it in these tutorials.

><details>
><summary>Variable types - path and tuple</summary>
>
> Don't worry about the keywords `path` and `tuple` in the above example. These are a keywords used by `Nextflow` to identify determine characteristics and expectations of the input and output files. This is not something we will be exploring in these tutorials.
></details>

With the **_process_** block completed, we now need to go to the **_workflow_** statement and tell `Nextflow` what to do with this **_process_**. Modify your **_workflow_** block to look like the following:

```diff
// Workflow
workflow {

+   input_files = Channel.fromPath("input_files/*.fq.gz")
+   map_to_reference(input_files)
}
```

This introduces two concepts at once. However, this is the worst part for understanding what is happening - for all subsequent processes the **_workflow_** requires only minor tweaking.

The first line of the command establishes a **_channel_** from a set of input files to be analysed by the workflow. We are ignoring the technical information on what a **_channel_** is for these tutorials, but see the documentation if you wish to understand them in detail.

For our purposes the important part is that the first line of the **_workflow_** is creating a list of files which match the shell expression `input_files/*R{1,2}.fq.gz`. Files detected this way are then organised in a manner such that they are grouped according to the sample they represent.

The second line takes the contents of the `input_files` list and applies the `map_to_reference` **_process_** to each entry in the list.

---

## Executing our first workflow

Once you have completed this workflow, you can execute it from the command line:

```bash
$ nextflow run nextflow_example.nf
```

This will produce some text on your console which will look something like the following:

```
N E X T F L O W  ~  version 22.10.3
Launching `nextflow_example.nf` [golden_allen] DSL2 - revision: 99157aa282
executor >  local (1)
[4f/87ba8d] process > map_to_reference (1) [100%] 1 of 1 ✔
```

If you check your current directory, you will see the although the `Nextflow` report claims that the workflow completed successfully, there is no output file produced. However there is a new folder present in your directory - `work/`. To find your output files, we need to make use of the text in the square brackets in the output. In this example the value is `4f/87ba8d`, you will have something different in your own output. Copy and paste this text into the command line to replicate the following command:

```bash
$ ls work/4f/87ba8d
```

Then use tab-completion to fill out the command:

```bash
$ ls work/4f/87ba8db01949863c26076aa52ff3b5/
```
```
mapping.sam
Mbovis_87900.nanopore.fq.gz
```

You will be able to match the fastq files to the contents of the `input_files/` folder. They are coloured in teal which is the shell's way of showing that these are shortcuts to the original files. The other file is our mapping output. The isolation of the input/outputs of each process execution into a randomly named subdirectory in the `work/` folder is a feature of `Nextflow` used to prevent multiple instances of the same rule from overwriting each others outputs. We will look at how to extract results from these randomly generated temporary directories later in this tutorial.

---

## Adding a second process to the workflow

Now that we have our first successful **_process_** lets add another one. We will this time take the sorting, filtering, and compression commands from the [mapping tutorial](./34_mapping_filters.md#sorting-and-compressing-sam-files). Open your `Nextflow` file and add the following process:


```diff
// Processes
process map_to_reference {
    ...
}

+process sort_and_filter {

+   input:
+   path sam_file

+   output:
+   path "mapped.bam"

+   script:
+   """
+   samtools view -bS ${sam_file} | samtools sort -o sorted.bam
+   samtools view -h -F 4 -b sorted.bam > mapped.bam
+   """
+}
```

The **_process_** should look familiar, as it is just taking commands we have already used and reproducing them within the `Nextflow` script block. You are hopefully also familiar enough with the form of the `input` and `output` blocks that you can see  that the input file is a variable which we must access in our `script` block and the output declares the name of the file we expect to see when the process is successfully completed.

The more interesting part is what will add to the **_workflow_** statement:

```diff
// Workflow
workflow {

    input_files = Channel.fromPath("input_files/*.fq.gz")
    map_to_reference(input_files)
+   sort_and_filter(map_to_reference.out)
}
```

We still declare our original input channel in the first line, then feed it into the `map_to_reference` process. But how do we get another channel to feed into the `sort_and_filter` process? This comes from the `map_to_reference.out` value. When each **_process_** completes, the file(s) produced from its `output` statement are collected in a fresh channel, which we can then feed into a new **_process_**. This is the fundamental building block of how we grow our workflow.

Run your new command, and ensure that it is working correctly. If so, you should see something similar to below:

```
N E X T F L O W  ~  version 22.10.3
Launching `nextflow_example.nf` [marvelous_sammet] DSL2 - revision: ab0e19254d
executor >  local (2)
[7a/d0d18c] process > map_to_reference (1) [100%] 1 of 1 ✔
[db/de7d4e] process > sort_and_filter (1)  [100%] 1 of 1 ✔
```

This is basically the same as before, but now we have a second **_process_** run, and a second set output directory created. Check the contents to make sure that `mapping.bam` is present before we continue. Note that there are two `bam` files in your output directory, but we have only marked the `mapped.bam` as a desired output. The `sorted.bam` file is ignored from this point forward.

---

## Adding the final process and exporting the results

We are now going to take a simple example from the [mapping statistics](./35_mapping_statistics.md) exercise and run the `samtools flagstats` command to produce our final output file. Add the following **_process_** to your workflow:

```diff
// Processes
process map_to_reference {
    ...
}

process sort_and_filter {
    ...
}

+process compute_flagstats {

+   publishDir "./", mode: "copy"

+   input:
+   path bam_file

+   output:
+   path "flagstats.txt"

+   script:
+   """
+   samtools flagstat ${bam_file} > flagstats.txt
+   """
+}
```

There is one new line in this statement that we haven't seen before - the `publishDir` line. As we know that this is the final stage of the workflow, we are telling `Nextflow` to make a copy of the output file in the executing directory. Hiding away the intermediate files in subdirectories is great as we usually do not need to see their contents. At some point we need to retrireve *something* from the workflow, or there is no point having run it.

You can add this tag to multiple processes, for example if you also wanted to copy out the `mapping.bam` file so that you had a mapping file and its summarise statistics produced in a single workflow.

> ### Exercise
>
> Since we're now getting quite proficient at writing `Nextflow` statements, have a go at modifying your **_workflow_** block to include the final **_process_**.
>
> <details>
> <summary>Solution</summary>
>
> ```diff
> workflow {
> 
>     input_files = Channel.fromPath("input_files/*.fq.gz")
>     map_to_reference(input_files)
>     sort_and_filter(map_to_reference.out)
> +   compute_flagstats(sort_and_filter.out)
> }
> ```
> </details>

If everything is successful, when you run your workflow for the last time we should see something like the following:

```
N E X T F L O W  ~  version 22.10.3
Launching `nextflow_example.nf` [pedantic_visvesvaraya] DSL2 - revision: dbbc67638e
executor >  local (3)
[c2/782a80] process > map_to_reference (1)  [100%] 1 of 1 ✔
[fd/1696c9] process > sort_and_filter (1)   [100%] 1 of 1 ✔
[95/275417] process > compute_flagstats (1) [100%] 1 of 1 ✔
```

Upon completion we should also have the `flagstats.txt` file in our directory with the following contents:

```
1733 + 0 in total (QC-passed reads + QC-failed reads)
1725 + 0 primary
0 + 0 secondary
8 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
1733 + 0 mapped (100.00% : N/A)
1725 + 0 primary mapped (100.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

---
