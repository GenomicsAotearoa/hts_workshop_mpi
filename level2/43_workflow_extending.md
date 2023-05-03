# Introduction to the Nextflow system

* Teaching: 60 minutes
* Exercises: 60 minutes

#### Objectives

* ...

#### Keypoints

* ...

---

## Contents

1. [Getting started](#getting-started)
1. [Internalising the module load commands](#internalising-the-module-load-commands)
1. [Streamlining the workflow block](#streamlining-the-workflow-block)
1. [Allowing dynamic input files](#allowing-dynamic-input-files)
1. [Generating output reports](#generating-output-reports)
1. [Other useful modifications](#other-useful-modifications)

---

## Getting started

In the previous tutorial we produced a simple `Nextflow` workflow to perform read mapping and export the results of a `samtools flagstats` command. This was a good first effort for producing a workflow but there were a few parts of the rule file which could be improved.

In this tutorial we are going to be iteratively improving our `Nextflow` workflow by adding a number of changes to the base file. Since it can be tricky to debug errors in `Nextflow` files, there are a series of completed examples given in the `.checkpoints/` folder which you can use to reinstate your current work if necessary.

Navigate to the `/nesi/project/nesi03181/phel/USERNAME/level2/workflows/` directory to begin this tutorial.

---

## Internalising the module load commands

When we first wrote out `Nextflow` command, the `bowtie2` and `samtools` modules were loaded outside of the workflow file, alongside the `Nextflow` module itself. While this approach does work, it produces a gap in our reproducibility as the workflow file cannot tell us which version of the mapping tools were used to execute the **_process_** blocks. If the workflow was successfully run we know that *some* version of those tools must have been present but it is anyone's guess which version they were.

This could potentially also result in a failure of the `Nextflow` workflow to run at all, if the available version of `bowtie2` or `samtools` version is not compatible with the commands we provide.

To remedy this, we are going to migrate our module load statements for the mapping tools inside the `Nextflow` workflow, so that we only need to load `Nextflow` itself for our workflow to run. To do this, make the following changes to your script.

```diff
process map_to_reference {

+    module "Bowtie2/2.4.5-GCC-11.3.0"

    input:
    tuple val(sample_id), path(paired_reads)

    output:
    path "mapping.sam"

    script:
    """
    bowtie2 --sensitive -x ${launchDir}/references/Mbovis_87900.genome -1 ${paired_reads[0]} -2 ${paired_reads[1]} -S mapping.sam
    """
}

process sort_and_filter {

+    module "SAMtools/1.16.1-GCC-11.3.0"

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

process compute_flagstats {

+    module "SAMtools/1.16.1-GCC-11.3.0"
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

We will still need to load the `Nextflow` module itself, as we cannot invoke `Nextflow` without `Nextflow` being present in our shell environment, but now we only need two commands to get the workflow up and running:

```bash
$ module load Nextflow/22.10.3
$ nextflow run nextflow_example.nf
```

Not only does this simplify the steps needed to get the workflow up and running, it also provides a form of quality control on the workflow itself, since the exact tools used in its execution are recorded in the workflow and these workflow file can be subject to version control.

>**Note:** This is actually the approach used in the Virology & Phytoplasmology team - the main `BLASTn` commands used for diagnostic testing are wrapped inside a `Nextflow` workflow, which is recorded in the MPI GitLab server.

---

## Streamlining the workflow block

The manner in which we wrote the **_workflow_** block in the previous tutorial is functional but not particularly elegant to read. While it does explicitly detail how the input channel is created, and the flow of data from the first process into the second and third it's quite dense.

`Nextflow` supports an alternate way of orchestrating the flow of data between processes, which is basically identical to redirection in the shell. Modify your **_workflow_** statement to the following:

```diff
workflow {

-    input_files = Channel.fromFilePairs("input_files/*R{1,2}.fq.gz")
-    map_to_reference(input_files)
-    sort_and_filter(map_to_reference.out)
-    compute_flagstats(sort_and_filter.out)
+    Channel.fromFilePairs("input_files/*R{1,2}.fq.gz") | map_to_reference | sort_and_filter | compute_flagstats
}
```

Once that is done, run your workflow again! This provides the same flow of data as the original version of the **_workflow_** block, but in a more succinct form. You may or may not prefer this implementation. If not, change it back.

---

## Allowing dynamic input files

One feature of the workflow we are currently working with is that the rule has a hardcoded statement which identifies input file(s) for the mapping processes. This means that if we want to run this workflow, we must create a folder of the name `input_files/` in the same directory as the rule and then copy our sequence files into the folder for the rule to pick them up.

This is limiting, as we may not want to be copying our data all around our file system when we need to perform mapping. We also probably don't want to be copying the workflow file into many different locations either. To avoid these workarounds, we are going to modify the workflow to allow for a user-provided input path. This requires that we change the workflow in two stages.

The first change we need to make is that we need to replace the `Channel.fromFilePairs('input_files/*R{1,2}.fq.gz')` to something which is influenced by user-specified input rather than set with a predetermined value. There is also a second set of changes we need to make.

> **Exercise**
>
> Take a moment to read through your current workflow and see if you can spot a problem in the current implementation. The current version is limited in its current form, and this will become a major problem in an updated version, particularly if we write the workflow in a manner that allows multiple input files.
>
> <details>
> <summary>Solution</summary>
>
> The issue is in the `compute_flagstats` process which declares a hardcoded output file name of `flagstats.txt`. When we are mapping with a single, known, input file this will not matter as we always know which file was used to produce the output. However, if we are changing the input code to allow variable input files we lose this traceability from input to output file.
>
> This is potentially even worse, as if we write the workflow in a way to allow multiple input files (which we are going to) then each output file will be copied to the path `./flagstats.txt` resulting in the output of the first completed mapping process to be overwritten by the second, then the second by the third etc.
>
> ```diff
> process compute_flagstats {
>
>    module "SAMtools/1.16.1-GCC-11.3.0"
>    publishDir "./", mode: "copy"
>
>    input:
>    path bam_file
>
>    output:
>-    path "flagstats.txt"
>
>    script:
>    """
>-    samtools flagstat ${bam_file} > flagstats.txt
>    """
> }
> ```
> </details>

Once we are done allowing for flexible inputs, we must account for this issue in the updated workflow.

To begin, we are going to ignore the file naming issue and just focus on changing the input code. This requires us to make use of a new concept, which is the [`params` scope](https://www.nextflow.io/docs/latest/config.html#scope-params). This is a feature through which users can provide dynamically written inputs to the workflow in a manner similar to the flags we provide to tools like `grep` and `cut`. We need to make two changes to the `Nextflow` rule to do this:

```diff
// Configuration
nextflow.enable.dsl=2
+params.input = "input_files/*R{1,2}.fq.gz"
```

```diff
workflow {

-    Channel.fromFilePairs("input_files/*R{1,2}.fq.gz") | map_to_reference | sort_and_filter | compute_flagstats
+    Channel.fromFilePairs(params.input) | map_to_reference | sort_and_filter | compute_flagstats
}
```

With these changes we have declared an optional user input flag called `input`. If the user provides a custom value for this flag then that will be used, otherwise the rule will default to the original `input_files/*.fq.gz` path. 

Now for the bigger issue. If there was a single file in the `input_files/` folder, and we only ever read inputs from that folder, then we would always know which output file corresponded to the input. We would probably need to manually rename the output after each run so that we knew which input file the output corresponded to, but that's manageable with a simple `mv` command. But we're aiming for an all-in-one solution so will not be doing that!

We are now going to solve this issue by updating the `Nextflow` code to dynamically name its output files, using the input name to provide some uniquely identifying text for each input file. Unfortunately, to do this we need to make a slight modification to each **_process_** block. Since we now have a fairly long rule, we have to make changes in multiple places. Lesson for the future, always build as you intend to go on!

```diff
process map_to_reference {

    module "Bowtie2/2.4.5-GCC-11.3.0"

    input:
    tuple val(sample_id), path(paired_reads)

    output:
+    val sample_id
    path "mapping.sam"

    script:
    """
    bowtie2 --sensitive -x ${launchDir}/references/Mbovis_87900.genome -1 ${paired_reads[0]} -2 ${paired_reads[1]} -S mapping.sam
    """
}
```

```diff
process sort_and_filter {

    module "SAMtools/1.16.1-GCC-11.3.0"

    input:
+    val sample_id
    path sam_file

    output:
+    val sample_id
    path "mapping.bam"

    script:
    """
    samtools view -bS ${sam_file} | samtools sort -o sorted.bam
    samtools view -h -F 4 -b sorted.bam > mapping.bam
    """
}
```

```diff
process compute_flagstats {

    module "SAMtools/1.16.1-GCC-11.3.0"
    publishDir "./", mode: "copy"

    input:
+    val sample_id
    path bam_file

    output:
-    path "flagstats.txt"
+    path flagstats_file

    script:
+    flagstats_file = "${sample_id}.flagstats.txt"
    """
-    samtools flagstat ${bam_file} > flagstats.txt
+    samtools flagstat ${bam_file} > ${flagstats_file}
    """
}
```

What we are really doing here is making use of the unique identifying `sample_id` value which was extracted in the original channel creation. We *could* make use of it to name each intermediate `sam` and `bam` file but this isn't really necessary for our purposes today. What we are doing instead is taking that original `sample_id` value in the `map_to_reference` **_process_** and passing it foward as an output value to the next **_process_**, `sort_and_filter`. In turn, this **_process_** passes the value forward to the final **_process_** `compute_flagstats` which then incorporates this variable into the output file name. We are retaining the uniquely identifying part of the file pairs and using these to name the final output file.

>**Note:** In practice it is often helpful to use variables to name the output file of each process as we did in the `compute_flagstats` command above. This makes it easier to trace errors and avoids transcription mistakes which can happen when updating the `Nextflow` code. We are avoiding this today to minimise the changes necessary to get a working example but remember this is writing your own workflows.

If everything is successful, we should now be able to just do the following:

```bash
$ nextflow run nextflow_example.nf
$ ls *.txt
```

```
Mbovis_87900.miseq.flagstats.txt
```

Running that `ls` command you should now see that rather an `flagstats.txt`, the output file is called `Mbovis_87900.miseq.flagstats.txt`. This gets even better though because now, if we make use of the `input` parameter we created, we can point the command towards a directory with multiple sequence files and get an output for each one with a single command:

```bash
$ nextflow run nextflow_example.nf --input 'other_input_files/*R{1,2}.fq.gz'
$ ls *.txt
```

```
Mbovis_87900.miseq.flagstats.txt
Mbovis.SRR20305625.flagstats.txt
Mbovis.SRR20305628.flagstats.txt
Mbovis.SRR20305631.flagstats.txt
```

Note that to do this, we wrap the search pattern with quotation marks. This is so that the input is passed directly into `Nextflow` as an expression to be evaluated. Without the quotes, the shell would expand the command into a list of files, which do not get handled the same way.

---

## Generating output reports

>TO DO - Easily achieved on the command line or with a config file, but looking for a way to internalise this into the workflow itself

```bash
$ nextflow run nextflow_example.nf -with-timeline other_input_files.timeline.html -with-reportother_input_files.report.html -with-dag other_input_files.dag.png
```

```
config:
timeline {
  enabled = true
  file = "$params.outdir/timeline.html"
}

report {
  enabled = true
  file = "$params.outdir/report.html"
}
```

---

## Other useful modifications

1. Generalising the input expression further
1. Adding adjustable CPU counts to the command
1. Using variables to name output files in processes rather than hardcoding until the final step
1. Adding multiple reference files for mapping?

---
