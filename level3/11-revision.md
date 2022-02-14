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
1. [slurm exercises](#slurm-exercises)

---

## Returning to NeSI

It has been quite a few months since the last training workshop, so the main point of the training today is to make sure everyone is up to speed on the main considerations when working with NeSI and high-performance computing. To begin, log into NeSI through the [JupyterHub web portal](https://jupyter.nesi.org.nz/) and start a session. Make sure that you use the training project (nesi03181) rather than the usual project.

Once you have logged in, open a terminal and navigate to **your** directory.

```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/
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

**Exercise 4.**

Finally, let's create a new `bash` script. The `bad-reads-script.sh` is not particularly useful, as it dumps all of the sequences into a snigle output file. For real work, we would probably prefer to keep a per-sample record of which sequences failed our filtering criteria.

Create a copy of the script produced in **Exercise 3**, then make the following changes:

1. Use a `for loop` to iterate over each `.fastq` file in your directory, instead of a wildcard
1. For each iteration of the loop, create a new output file to store the reads (i.e. each `.fastq` file should have a separate output file for bad reads rather than a single output file for all samples).
   1. Remember that you can use the `basename` command to extract a unique identifier from the loop variable

---

## slurm exercises

Now that we're back up to speed with `bash` and the basics of command line navigation and operation, let's take a look at `slurm`.

Navigate over to your folder named `1_Raw_data`. We are going to write  a quick `slurm` script to perform quality filtering over one of the MiSeq samples in this directory, and one of the MinION samples.

**Exercise 1.**

Firstly, we need to identify appropriate tools for performing filtering on these data sets and ensure that they are available on NeSI.

Which tools can we use for filtering Illumina and Oxford data? Why can we not use the same tool for both types of data?

**Exercise 2.**

We now know the names of the tools for performing the quality filtering, but we do not yet have access to them from the command line. How do we identify what software is available in NeSI, and load it into our environment?

**Exercise 3.**

Write a `slurm` script to execute your MiSeq filtering tool on one of the samples in your `1_Raw_data` folder. You can use the following `slurm` template to begin your file:

```bash
#!/bin/bash -e
#SBATCH --account       nesi03181
#SBATCH --job-name      revision_filtering
#SBATCH --time          00:30:00
#SBATCH --cpus-per-task 10
#SBATCH --mem           10G
#SBATCH --error         revision_filtering.%j.err
#SBATCH --output        revision_filtering.%j.out
#SBATCH --mail-type     END
#SBATCH --mail-user     YOUR_EMAIL

module purge
module load ...

cd /nesi/project/nesi03181/phel/USERNAME/1_Raw_data/

# Execute your command
# ...
```

>**Note:** You can keep all filtering parameters at the defaul settings for this exercise if you wish, we are more interested in getting the `slurm` script working for this set of exercises.

---

