# Revision

* Teaching: 30 minutes
* Exercises: 30 minutes

#### Objectives

* Revisit the basics of working in a HPC environment
* Refamiliarise with the project directory
* Practice with simple command line tools
* Submit simple `slurm` scripts to NeSI

---

## Contents

1. [Returning to NeSI](#returning-to-nesi)
1. [Shell exercises](#shell-exercises)
1. [xxx](#xxx)

---

## Returning to NeSI

It has been quite a few months since the last training workshop, so the main point of the training today is to make sure everyone is up to speed on the main considerations when working with NeSI and high-performance computing. To begin, log into NeSI through the [JupyterHub web portal](https://jupyter.nesi.org.nz/) and start a session. Make sure that you use the training project (nesi03181) rather than the usual project.

Once you have logged in, open a terminal and navigate to **your** directory.

```bash
$ cd /nesi/project/nesi03181/phel/YOUR_USER_NAME/
```

We are going to do a few quick exercises to refresh on the basics of command line operations.

You should still have a folder named `shell_data/` in your directory. This is where we did our initial practice introduction to the command line and common operations. Ignoring any folders you created during your own working, there were three folders provided to you at the start of this training workshop.

> ### Exercise
>
> There are three folders in your `shell_data/` directory. One of these is hidden, the other two are not. List these folders, and report back with their names and how long ago they were created.
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> $ cd shell_data/
> # ls -a -l
> drwxrws---+  2 dsen018        nesi03181 4096 Jun  2  2021 .hidden
> drwxrws---+  3 dsen018        nesi03181 4096 Jul 11  2021 sra_metadata
> drwxrws---+  3 dsen018        nesi03181 4096 Jun 23  2021 untrimmed_fastq
> ```
> </details>

From here, we will perform a few more exercises to refresh ourselves on some of the common command line operations.

---

## Shell exercises

There are no answers given for this next section - you can either work through this section independently, or with the group.

**Exercise 1.**

How many lines are in the files `SRR097977.fastq` and `SRR098026.fastq`? Based on this information, how many sequences are in each file?

**Exercise 2.**

Sometimes we need to peak inside large files to get a view of the content. We don't neccessarily want to open the full file, but getting a glimpse of the first few lines can be very useful.

Write a command to extract the first 8 lines of `SRR097977.fastq` and write them into a new file, named `SRR097977.sample.fastq`

**Exercise 3.**

We previously wrote a script called `bad-reads-script.sh`, which reads the contents of each fastq file in a location and copies out sequences with a run of 10 `N` characters. The results are then redirected to a new file, called `scripted_bad_reads.txt`.

Perform the following tasks:

1. Create a copy of the `bad-reads-script.sh` script, using whatever name you like.
1. In the copied script, modify the first `grep` command in the script to extract the sequence name and nucleotide sequence only (no quality information).
1. In the copied script, change the destination file to redirect into a file with the extension fasta.

Once you have completed this exercise, take a look and the contents of the output file using the `less` or `cat` command. Is the output a valid fasta file? If not, what would you need to change to make it one?

---
