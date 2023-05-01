# Introduction to the Nextflow system

* Teaching: 30 minutes
* Exercises: 30 minutes

#### Objectives

* Build a small workflow in `Nextflow` to map sequences to a reference and produce reports from the execution of the workflow.
* Know how to produce summary reports and DAG outputs from a `Nextflow` workflow.

#### Keypoints

* Know that a complex workflow is built from smaller sub-units of work, which are executed sequentially to produce a final outcome.
* Understand the main sections of the `Nextflow` workflow and how the `workflow{}` block governs the individual processes.

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

```bash
// Configuration
nextflow.enable.dsl=2

// Processes


// Workflow
workflow {

}
```

The lines which behin with `//` are comments - these are statements which are not read by the `Nextflow` engine and exist for our own benefit. Commenting in code is a useful way to break the lines up with notes to yourself or other users. You can use comments to explain what a particular region of a script is used for, or to explain the rationale for a particular implementation over another. In this instance we are using three comments to break our script up into three sections.

The first section is to do with the workflow configuration. We are inserting a single line here, which specifies to `Nextflow` that the workflow must be interpretted using the second release of its *domain-specific language*, which is the language style `Nextflow` uses to express the code blocks and relationships between them. Although the second version of this language has been released for some time, the `Nextflow` engine still defaults to version 1 when interpretting commands so when producing a workflow we need to specify if we are deviating from this. There are three ways to tell `Nextflow` which version of the DSL we are using:

1. Add the `-dsl2` flag to the command line when executing the workflow.
1. Add the line `nextflow.enable.dsl=2` inside the workflow rule to overwrite the user-provided or default case for DSL interpretation.
1. Add the statement `nextflow.enable.dsl=2` to a configuration file which sits in the same directory as the `Nextflow` rule being executed.
   1. A configuration file must have the name `nextflow.config` and reside in the same directory as the `.nf` file to be recognised.
   1. For large projects, or projects which can work across different server deployments this is the prefered way to address configuration settings but it is overly complex at this stage in our lesson.

The second comment is simply a placeholder for where we are going to write **_process_** statements. These are effectively the individual shell commands which we will write to perform steps of our analysis.

The third comment denotes the `workflow{}` statement. The **_workflow_** statement is the brains of the `Nextflow` rule. It assesses which input file(s) are to be analysed and then assembles the individual **_process_** statements in the correct order to achieve their outcome.

---

## Writing the first process in the workflow

Now that we have the outline of our workflow we are going to take the `minimap2` reference mapping command from [the previous exercise of Nanopore mapping](./33_nanopore_mapping.md#mapping-reads-with-minimap2) and transcribe it into `Nextflow`.

To begin, we are going to modify your basic document to include a single process, which maps a given input file against a pre-determined reference file.

>**Note:** If at any time you get lost, there are a set of completed `Nextflow` rules for each stage of this tutorial in a hidden folder named `.checkpoints/`

To begin, under your comment 'Processes', write the following commands:

```bash
// Processes
process map_to_reference {

    input:
    path fasta_file

    output:
    path "mapping.sam"

    script:
    """
    minimap2 -ax map-ont ${launchDir}/references/Mbovis_87900.genome.mmi ${fasta_file} > mapping.sam
    """
}
```

Each **_process_** is declared by the statement `process` followed by a unique name used to identify the process block. The process is defined as all file content between the opening and closing curly brace (`{`, and `}`). If you are familiar with programming in a C-like language you will recognise this as very similar to a function declaration.

Within the **_process_** block are three headers, `input`, `output` and `script`. What these refer to should be intuitive by their naming:

* `input:`
   * This block refers to the input file(s) for the process.
   * Any entries under the `input` block are provided as variables as we do not know in advance what file names are going to be provided when the workflow is executed.
* `output:`
  * This block defines the name of any output files which must be produced upon *successful* completion of the **_process_**.
  * We have three choices for setting the name of any output files produced.
    1. Hardcode the names as string values (as done here).
    1. Infer the output names from the input file (which we will do in the next tutorial)
    1. Pull the names from user-provided variables.
 * `script:`
   * This section is where we write the actual command to be executed.
   * By default this is read as the `bash` language you have been using in the tutorials leading up to this one so you should hopefully recognise the `${fasta_file}` style of accessing variables.
   * It is possible to replace this block with a different statement, `shell:` which you may see in some online examples.
     * This is mainly done when the code block needs to be written in a language other than `bash`, such as `perl`, `python`, or `R`.
     * Using the `shell:` statement changes the way in which variables are accessed so we will not be using it in these tutorials.

>**:Note:** Don't worry about the keyword `path` in the above example. This is a keyword used by `Nextflow` to identify particular characteristics and expectations of the input and output file names. While there are other keywords that can be used, we will not be exploring them in these tutorials.

With the **_process_** block completed, we now need to go to the **_workflow_** statement and tell `Nextflow` what to do with this **_process_**. Modify your **_workflow_** block to look like the following:

```bash
// Workflow
workflow {

    input_files = Channel.fromPath('input_files/*.fna')
    map_to_reference(input_files)
}
```

This introduces two concepts at once. However, this is the worst part for understanding what is happening - for all subsequent processes the **_workflow_** requires only minor tweaking.

The first line of the command establishes a **_channel_** from a set of input files to be analysed by the workflow. We are ignoring the technical information on what a **_channel_** is (see the documentation if you wish to understand them in detail) but for our purposes all that we need to understand is that the first line of the **_workflow_** is creating a list of files which match the pattern `input_files/*.fna` (i.e. files in the `input_files/` directory which end with the `.fna` extension).

The second line takes the contents of the `input_files` list and applies the `map_to_reference` **_process_** to each entry in the list.

---

## Executing our first workflow

Once you have completed this workflow, you can execute it from the command line:

```bash
$ nextflow run your_command.nf
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
temp.fna
```

You will be able to match the `temp.fna` file to the contents of the `input_files/` folder and the other file is our mapping file. The isolation of the input/outputs of each process execution into a randomly named subdirectory in the `work/` folder is a feature of `Nextflow` used to prevent multiple instances of the same rule from overwriting each others outputs and is also used in resuming aborted runs. We will look at how to extract out results from these randomly generated temporary directories later in this tutorial.

---

## Adding a second process to the workflow

Now that we have our first successful **_process_** lets add another one. We will this time take the sorting, filtering, and compression commands from the [mapping tutorial](./34_mapping_filters.md#sorting-and-compressing-sam-files). Open your `Nextflow` file and add the following process:


```bash
process sort_and_filter {

    input:
    path sam_file

    output:
    path "mapping.bam"

    script:
    """
    samtools view -bS ${sam_file} | samtools sort -o sorted.bam
    samtools view -h -F 4 -b sorted.bam > mapping.bam
    """
}
```

One you are done, also add the following line to the end of your **_workflow_**. Make sure that this is entered before the closing `}` character of the **_workflow_** block, but after the `map_to_reference(input_files)` line.

```bash
sort_and_filter(map_to_reference.out)
```

The **_process_** should look familiar, as it is just taking commands we have already used and reproducing them within the `Nextflow` script block. You are hopefully also familiar enough with the form of the `input` and `output` blocks that you can see  that the input file is a variable which we must access in our `script` block and the output declares the name of the file we expect to see when the process is successfully completed.

The more interesting part is what we have added to the **_workflow_** statement. In it's entirety, the **_workflow_** block now looks like:

```bash
// Workflow
workflow {

    input_files = Channel.fromPath('input_files/*.fna')
    map_to_reference(input_files)
    sort_and_filter(map_to_reference.out)
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

This is basically the same as before, but now we have a second **_process_** run, and a second set output directory created. Check the contents to make sure that `mapping.bam` is present before we continue.

---

## Adding the final process and exporting the results

We are now going to take a simple example from the [mapping statistics](./35_mapping_statistics.md) exercise and run the `samtools flagstats` command to produce our final output file. Add the following **_process_** to your workflow:

```bash
process compute_flagstats {

    publishDir "./", mode: "copy"

    input:
    path bam_file

    output:
    path "flagstats.txt"

    script:
    """
    samtools flagstat ${bam_file} > flagstats.txt
    """
}
```

There is one new line in this statement that we haven't seen before - the `publishDir` line. As we know that this is the final stage of the workflow, we are telling `Nextflow` to make a copy of the output file in the executing directory. Hiding away the intermediate files in subdirectories is great as in this context we don't really need to see the contents of the mapping and filtering commands. However since we are interested in the final mapping statistics it is not ideal that this file is also hidden aawy in a cryptic `work/` directory.



> ### Exercise
>
> Since we're now getting quite proficient at writing `Nextflow` statements, have a go at modifying your **_workflow_** block to include the final **_process_**.
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> workflow {
> 
>     input_files = Channel.fromPath('input_files/*.fna')
>     map_to_reference(input_files)
>     sort_and_filter(map_to_reference.out)
>     compute_flagstats(sort_and_filter.out)
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

Upon completion we should also have the `flagstats.txt` file in our directory.

<details>
<summary>Contents of `flagstats.txt`</summary>
```
39 + 0 in total (QC-passed reads + QC-failed reads)
39 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
39 + 0 mapped (100.00% : N/A)
39 + 0 primary mapped (100.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
</details>

---
