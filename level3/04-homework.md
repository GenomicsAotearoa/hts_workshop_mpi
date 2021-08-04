# Level 3.1 revision

This revision material consists of three short exercises covering the main points of the first part of Level 3 training. For each of these exercises you will be performing read mapping against the reference genome downloaded in the Level 2 homework. You can either copy the reference genome (`GCF_000696015.1_ASM69601v1_genomic.fna`) into your `3_Assembly_mapping/` directory for these exercises, or index and map against the current location.

Two of these exercises require you to write `slurm` scripts to perform read mapping. Each exercise should require around 30 minutes of run time, but you are welcome to use more time as different nodes on NeSI have different processing speeds and it is hard to accurately predict the minimum time needed for a job.

When you have completed all exercises, email the requested material and answers to the trainers.

**_All participants must complete the homework to facilitate the smooth running of the sessions._**

> Remember, if you are struggling with any exercises:
>
> 1. Every exercise here was covered in the training material. Refer to this for hints.
> 1. It is not cheating to use Google to find help on how various commands work.
> 1. You are welcome to contact the trainers for advice on a particular exercise, but please attempt the first two options before resorting to this.

---

## Exercise 1 - Mapping to a genome using `bowtie2`

Create a `slurm` batch file to perform mapping of the quality filtered sequences corresponding to the `Mb168` isolate against the reference genome and convert the output to a sorted `bam` file.

For the sake of writing your `slurm` file, `bowtie2` should need about 30 minutes to perform read mapping with 10 threads and to then sort and compress the `bam` file. You can search the `bowtie2` manual to find the command for changing the thread limit.

Once you have completed this task, download both the reference genome and your `bam` file and import them into `Geneious`. Extract the consensus sequence from the mapped output and send this to the trainers along with your `slurm` script as the results of this exercise.
          
---

## Exercise 2 - Mapping to a genome using `bwa`

Create a `slurm` batch file to perform mapping of the quality filtered sequences corresponding to the `Mb168` isolate against the reference genome and convert the output to a sorted `bam` file.

>Note: Remember that `bwa mem` cannot map paired and unpaired reads in a single command, so you will need to perform two `bwa mem` steps and then merge the outputs into the final file.

Similar to with `bowtie2`, you can increase the number of threads for `bwa mem`. This will require approximately 25 minutes to map both data sets, perform the merge, and then compress the results as a `bam` file.

Once you have completed this task, download the `bam` file and import it into `Geneious`. Extract the consensus sequence from the mapped output and send this to the trainers along with your `slurm` script as the results of this exercise.

---

## Exercise 3 - Mapping to a genome using `minimap2`

Map the quality filtered sequences corresponding to the `Mb168` isolate against the reference genome and convert the output to a sorted `bam` file. For this job, the resource requirement is very low so you can perform mapping from the command line rather than submit a `slurm` job.

Once you have obtained your `bam` file, use your knowledge of the `samtools` software to split the file into two additional `bam` files, separating the mapped and unmapped reads.

Finally, you are going to convert the `bam` files into `fastq` files in order to determine how many reads were successfully mapped, and how many failed to map. In this case the command to export a `bam` file to `fastq` format is:

```bash
$ samtools fastq Mb168.minimap2.bam > Mb168.minimap2.fastq
```

>Note: This in a very simple use of the `samtools fastq` command and it will not work for Illumina data where we have paired and unpaired reads to sort through.

Perform this for both of your `bam` files then ascertain the number of sequences in each `fastq` file in whatever way seems appropriate to you. Report the number of sequences in each file to the trainers as your results for this exercise.

---
