# 3.4 - Filtering and compressing sam files

!!! clock "time"

    * Teaching: 20 minutes
    * Exercises: 20 minutes

!!! circle-info "Objectives and Key points"

    #### Objectives
    
    * Use `samtools` to sort and compress a raw `sam` file into the `bam` format.
    * Use `samtools` to filter a `bam` file into either the successfully mapped, or unmapped reads.
    * Use `samtools` to recover reads in `fastq` format from a `bam` file.
    
    #### Keypoints
    
    * Understand the reasons for sorting and compressing files in the `sam` and `bam` formats.
    * Understand the situations in which you may wish to filtering a `sam`/`bam` file and what the downstream applications of the output would be.

---

## Contents

1. [Why we need to compress and filter `sam` files](#why-we-need-to-compress-and-filter-sam-files)
1. [Sorting and compressing `sam` files](#sorting-and-compressing-sam-files)
1. [Splitting `bam` files to separate mapped and unmapped reads](#splitting-bam-files-to-separate-mapped-and-unmapped-reads)

---

## Why we need to compress and filter `sam` files

Now that we have created some basic mapping data in the `sam` files from the last exercises([here](./32_illumina_mapping.md) and [here](./33_nanopore_mapping.md)), let's take a look at the size of these files. Navigate to the `/nesi/project/nesi03181/phel/USERNAME/level2/mapping/` directory and run the following command:

```bash
ls -sh *.sam
```

You should see an output similar to the following:

```
1.3G Mbovis_87900.16S_rRNA.bowtie2.sam
1.5G Mbovis_87900.genome.bowtie2.sam
52M Mbovis_87900.genome.nanopore.sam
```

These are quite large files - especially in the case of the 16S rRNA mapping file where the mapping success rate was only around 0.4% of the Illumina sequences. Why is this file so far? There are two reasons.

1. The `sam` file format is a human-readable file format.
   * You can view the contents of this file with the `less` or `head` command (**do not** use `cat`), and see that although it's a bit cryptic, you can read recognise some values and patterns in the data such as file names and sequence content.
   * This is valuable when we need to read the file directly but is inefficient in terms of the hard drive space used.
1. When mapping is performed, *all* of the reads in the input file are recorded in the output `sam` file, along with a marker denoting whether or not they were successfully mapped to the reference.

This second point is the major cause of the file size. Although there are only about 8,700 reads actually mapped to the reference sequence (we'll see how to calculate that number in the next session), the `sam` file still contains all 2 million reads from the Illumina library.

This is a deliberate design feature of the mapping tools as depending on our context we may prefer to have the mapped or unmapped reads. For example, if we are trying to perform a reference-based genome assembly or determine the sequencing coverage of the genome (or a particular region) having the mapped reads is critical. However, if we are trying to subtract the host reads from a set of reads, such as when we are trying to find a pathogen, the reads which do not map to the host genome are actually the ones we are interested in examining in further detail.

---

## Sorting and compressing `sam` files

The first step of reducing the file size is to efficienctly compress the contents of the `sam` file. We could technically do this with a compression tool like `gzip` (see [exercise 1.1](./11_shell_manipulation.md#compressing-and-uncompressing-large-files)) but this would then make it difficult to read and manipulate the contents of the file. Fortunately there is a built-in solution for this - the `sam` file specification also has a binary-encoded equivalent which records the exact same information, in a much more efficieny format. This is the Binary Alignment/Map (`bam`) format.

When performing this compression from `sam` to `bam` we also use this opportunity to sort the mapped reads in terms of their starting position in the reference sequence. This sorting is important as it increases the speed of many of the operations we need to perform using a `sam` file, particular when producing coverage statistics. Because the sorting and compression can be performed in a single command line operation we tend to do these things together once and never worry about it again.

>In the situations where you need your mapped reads sorted, operations will fail if the reads are not sorted. In situations where you do not need them sorted, operations will succeed in either case so it is generally just easier to sort as soon as mapping completes then never worry about it again.

the main tool used for handling `sam` and `bam` files is the `samtools` package (). Load the module and execute the file below. While it is running, refresh yourself on what the `|` operator in the command is doing.

```bash
$ module load SAMtools/1.16.1-GCC-11.3.0

$ samtools view -bS Mbovis_87900.16S_rRNA.bowtie2.sam | samtools sort -o Mbovis_87900.16S_rRNA.bowtie2.bam
```

The reason we need to redirect the data from one `samtools` command to another is due to the behaviour of the `samtools sort` subcommon. If you examine the manual for this command, you will see that while the *output* of the command can be written in `sam` or `bam` format the *input* must be in `bam`. We therefore need to use the `samtools view` subcommand to first convert the `sam` file into the `bam` format. Rather than write the results to NeSI's hard drive then perform `samtools sort` as a second command we can redirect between the commands to save on hard drive space and speed up the operation.

Compare the file sizes between the `sam` and `bam` files:

```bash
ls -sh Mbovis_87900.16S_rRNA.bowtie2.sam Mbovis_87900.16S_rRNA.bowtie2.bam
```

```
1.3G Mbovis_87900.16S_rRNA.bowtie2.sam
470M Mbovis_87900.16S_rRNA.bowtie2.bam
```

> **Note:** There is no point in retaining the original `sam` file at this point, as the information it contains is more efficiently encoded within the `bam` format. If working with real data, this is the point you should delete your `sam` file. In this training exercise, it does not matter whether you delete the `sam` file or not.

By compressing the data of the `sam` file into the `bam` file we have already reduced its size to about one third of the original. Not only will this save us space on NeSI, it also makes downloading the data much quicker.

> ### Exercise
>
> Sort and compress any other `sam` files you have produced during the previous mapping tutorials.
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> $ for i in bowtie2 nanopore;
> > do
> >     samtools view -bS Mbovis_87900.genome.${i}.sam | samtools sort -o Mbovis_87900.genome.${i}.bam
> > done
> ```
> </details>

---

## Splitting `bam` files to separate mapped and unmapped reads

Most of the time when we are mapping against a reference sequence, we are interested in the sequences which successfully mapped to the target. In these cases, having a `bam` file which contains all of the unmapped reads is not particularly useful so we will now practice filtering `bam` files according to the mapping state of the reads in the file.

This is simple to do in terms of the command which needs to be run, but it the meaning of the command can be a bit confusing. Within the `sam` and `bam` file format specification is a numeric flag which may be assign to an unmapped read, denoting its status as such. Because it is a *positive* marker of unmapped state, if we want to filter for mapped reads we need to filter for sequences which do not have the marker.

When running the `samtools view` command there is a pair of optional paramters which allow us to filter reads by their mapping flags. We use the `-f` to *keep* reads with the flag(s) we specify, or `-F` to *reject* reads with the flag(s) we specify. As the unmapped flag is only applied to reads which fail to map to the reference we use the following logic:

* `-f 4` = Include reads with the unmapped flag = Keep unmapped reads
* `-F 4` = Reject reads with the unmapped flag = Keep mapped reads

>**Note** The value `4` is the numeric flag for ummapped reads.

We will run `samtools` to filter the `Mbovis_87900.16S_rRNA.bowtie2.bam` file to only keep mapped reads:

```bash
$ samtools view -h -F 4 -b Mbovis_87900.16S_rRNA.bowtie2.bam > Mbovis_87900.16S_rRNA.bowtie2.mapped.bam
```

In this command, the difference between keeping mapped or unmapped reads is the case of the `-f` flag. It is really easy to get confused about these values so be very careful when applying this command. Of course, in this case we are expecting only a handful of reads to have actually mapped to the reference sequence so a good acid test for whether we used the right commnad will be to check the size of the output file. If it is roughly the same size as the input then we have probably kept the unmapped reads. If it is a lot smaller, we have most likely kept only the mapped reads.

```bash
$ ls -sh Mbovis_87900.16S_rRNA.bowtie2.bam Mbovis_87900.16S_rRNA.bowtie2.mapped.bam
```

```
470M Mbovis_87900.16S_rRNA.bowtie2.bam
1.8M Mbovis_87900.16S_rRNA.bowtie2.mapped.bam
```

> ### Exercise
>
> Write a loop to filter the genome mapping files, creating an output file for both the mapped and unmapped reads for each case.
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> $ for i in bowtie2 nanopore;
> > do
> >     samtools view -h -F 4 -b Mbovis_87900.genome.${i}.bam > Mbovis_87900.genome.${i}.mapped.bam
> >     samtools view -h -f 4 -b Mbovis_87900.genome.${i}.bam > Mbovis_87900.genome.${i}.unmapped.bam
> > done
> ```
> </details>

---
