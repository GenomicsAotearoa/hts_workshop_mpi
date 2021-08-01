# Mapping Oxford Nanopore sequences to a reference

* Teaching: 20 minutes
* Exercises: 10 minutes

#### Objectives

* Use `minimap2` to map reads against a reference genome or sequence.
* Use `samtools` to sort and compress the alignment file.
* Import and visualise the results using `Geneious`.

#### Keypoints

* There are specific tools for aligning reads against a reference sequence when working with Nanopore sequences.

---

## Contents

1. [Mapping reads with `minimap2`](#mapping-reads-with-minimap2)
1. [Filtering unmapped reads](#filtering-unmapped-reads)

---

## Mapping reads with `minimap2`

Once more, navigate to the working directory before beginning the exercise.

```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/3_Assembly-mapping/
```

As in the previous module, we will create a persistent index of the reference sequence, to accelerate the repeated mapping attempts.

```bash
$ module load minimap2/2.20-GCC-9.2.0

$ mkdir minimap_index/
$ minimap2 -d minimap_index/NZ_CP005933_16S_rRNA.mmi minimap_index/NZ_CP005933_16S_rRNA.fasta
```

Once the index is created we will map the Nanopore reads using `minimap2` with the Oxford Nanopore Technologies settings. These are toggled using the `-x` parameter.

```bash
$ mkdir minimap_results/

$ minimap2 -ax map-ont minimap_index/NZ_CP005933_16S_rRNA.mmi ../2_Quality_filtered_data/Mb152_trimmed.minion.fastq > minimap_results/Mb152.16S_M_bovis.ont.sam
```

We can then sort and compress the results in exactly the same way as for Illumina data.

```bash
$ module load SAMtools/1.12-GCC-9.2.0

$ samtools view -bS minimap_results/Mb152.16S_M_bovis.ont.sam | samtools sort -o minimap_results/Mb152.16S_M_bovis.ont.bam
```

The results can now be imported into `Geneious` just as before, but we will not do this just yet.

> ### Exercise
>
> Write a loop to perform the `minimap2` and `samtools` operations over the other two samples.
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> $ module load minimap2/2.20-GCC-9.2.0
> $ module load SAMtools/1.12-GCC-9.2.0
>
> $ for i in Mb1 Mb168;
> do
>     minimap2 -ax map-ont minimap_index/NZ_CP005933_16S_rRNA.mmi ../2_Quality_filtered_data/${i}_trimmed.minion.fastq > minimap_results/${i}.16S_M_bovis.ont.sam
>     samtools view -bS minimap_results/${i}.16S_M_bovis.ont.sam | samtools sort -o minimap_results/${i}.16S_M_bovis.ont.bam
> done
> ```
> </details>

---

## Filtering unmapped reads

As you will have seen by now when importing `bam` files into `Geneious`, the mapping file contains not only the coordinates of all mapped sequences but also a list of the reads which failed to map to the reference.

As we are performing alignment against only the 16S rRNA gene of *M. bovis* the majority of our sequences are not mapped. This produces files which are much larger than necessary so in the interests of speed we can remove unmapped reads from the `bam` file.

```bash
$ module load SAMtools/1.12-GCC-9.2.0

$ samtools view -h -F 4 -b minimap_results/Mb152.16S_M_bovis.ont.bam > minimap_results/Mb152.16S_M_bovis.ont_mapped.bam
$ samtools view -h -f 4 -b minimap_results/Mb152.16S_M_bovis.ont.bam > minimap_results/Mb152.16S_M_bovis.ont_unmapped.bam
```

These two commands look very similar and they only differ by a single parameter - the `-f` versus the `-F`. If you recall the table of `sam` file contents which was presented at the start of this module, the second column of the `sam` and `bam` file is the `FLAG` column which contains numeric values representing different information about the mapping. The value `4` tells us that read is **ummapped**.

The `-f` and `-F` are filters for whether to retain or reject sequences that possess the `FLAG` value that follows. Therefore:

* `-f 4` = Include reads **with** the unmapped flag = Keep unmapped reads
* `-F 4` = Include reads **without** the unmapped flag = Keep mapped reads

It is really easy to get confused about these values...

> ### Exercise
>
> Write a loop to split the remaining `bam` files into mapped and unmapped files.
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> $ module load SAMtools/1.12-GCC-9.2.0
>
> $ for i in Mb1 Mb168;
> do
>     samtools view -h -F 4 -b minimap_results/${i}.16S_M_bovis.ont.bam > minimap_results/${i}.16S_M_bovis.ont_mapped.bam
>     samtools view -h -f 4 -b minimap_results/${i}.16S_M_bovis.ont.bam > minimap_results/${i}.16S_M_bovis.ont_unmapped.bam
> done
> ```
> </details>

Once done, we can see the difference in file size with the following command:

```bash
$ ls -sh minimap_results/*.bam
```

---
