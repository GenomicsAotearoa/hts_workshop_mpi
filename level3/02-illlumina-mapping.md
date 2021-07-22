# Mapping Illumina sequences to a reference

* Teaching: 20 minutes
* Exercises: 10 minutes

#### Objectives

* Use `bowtie2` or `bwa` to map paired-end reads against a reference genome or sequence.
* Use `samtools` to sort and compress the alignment file.
* Import and visualise the results using `Geneious`.

#### Keypoints

* There are multiple viable tools for aligning reads against a reference sequence.

---

## Contents

1. [Mapping reads with `bowtie2`](#mapping-reads-with-bowtie2)
1. [Mapping reads with `bwa`](#mapping-reads-with-bwa)

---

## Mapping reads with `bowtie2`

<Text to come>
  
```bash
# Navigate to the correct location
$ cd /nesi/project/nesi03181/phel/USERNAME/3_Assembly-mapping/

# Create an index of the reference sequence
$ module load Bowtie2/2.4.1-GCC-9.2.0
$ mkdir bowtie2_index/
$ bowtie2-build NZ_CP005933_16S_rRNA.fasta bowtie2_index/NZ_CP005933_16S_rRNA
```

We are now ready to align (map) the sequences against the reference.

```bash
$ mkdir bowtie2_results/
$ bowtie2 --sensitive -x bowtie2_index/NZ_CP005933_16S_rRNA \
          -1 ../2_Quality_filtered_data/Mb152_1_trimmed.miseq.fastq \
          -2 ../2_Quality_filtered_data/Mb152_2_trimmed.miseq.fastq \
          -U ../2_Quality_filtered_data/Mb152_unpaired_trimmed.miseq.fastq \
          -S bowtie2_results/Mb152.16S_M_bovis.sam
```

We will now sort and compress the contents of the `sam` file. This will not require us to perform a `module purge` as `bowtie2` and `samtools` were built with the same environment, but there is no harm in performing a purge if you wish.

```bash
$ module load SAMtools/1.12-GCC-9.2.0
$ samtools view -bS bowtie2_results/Mb152.16S_M_bovis.sam | samtools sort -o bowtie2_results/Mb152.16S_M_bovis.bam
```

> **Note:** There is no point in retaining the original `sam` file at this point, as the information it contains is more efficiently encoded within the `bam` format.

We can now download the resulting `bam` file and import it into `Geneious`.

> ### Exercise
>
> You have also obtained paired-end reads from two other *M. bovis* isolates. Write a loop to perform the `bowtie2` and `samtools` operations over these data.
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> $ module load Bowtie2/2.4.1-GCC-9.2.0
> $ module load SAMtools/1.12-GCC-9.2.0
>
> $ for i in Mb1 Mb168;
> do
>      bowtie2 --sensitive -x bowtie2_index/NZ_CP005933_16S_rRNA \
>              -1 ../2_Quality_filtered_data/${i}_trimmed.miseq.fastq \
>              -2 ../2_Quality_filtered_data/${i}_2_trimmed.miseq.fastq \
>              -U ../2_Quality_filtered_data/${i}_unpaired_trimmed.miseq.fastq \
>              -S bowtie2_results/${i}.16S_M_bovis.sam
>     samtools view -bS bowtie2_results/${i}.16S_M_bovis.sam | samtools sort -o bowtie2_results/${i}.16S_M_bovis.bam
> done
> ```
> </details>

---

## Mapping reads with `bwa`

<Text to come>

```bash
# Create an index of the reference sequence
$ module load BWA/0.7.17-gimkl-2017a
$ mkdir bwa_index/
$ bwa index -p bwa_index/NZ_CP005933_16S_rRNA NZ_CP005933_16S_rRNA.fasta
```

We are now ready to align (map) the sequences against the reference. Unforetunatey, `bwa` does not allow us to map paired and unpaired in the same operation. We can still achieve the same outcome as when mapping with `bowtie2` but must perform a few more steps.

```bash
$ mkdir bwa_results/
$ bwa mem bwa_index/NZ_CP005933_16S_rRNA \
      ../2_Quality_filtered_data/Mb152_1_trimmed.miseq.fastq \
      ../2_Quality_filtered_data/Mb152_2_trimmed.miseq.fastq > bwa_results/Mb152.16S_M_bovis.paired.sam
$ bwa mem bwa_index/NZ_CP005933_16S_rRNA \
      ../2_Quality_filtered_data/Mb152_unpaired_trimmed.miseq.fastq > bwa_results/Mb152.16S_M_bovis.unpaired.sam
```

We will now sort and compress the contents of the `sam` file. This time, however, we will need to account for the fact that there are two `sam` files that need to be combined.

```bash
$ module load SAMtools/1.12-GCC-9.2.0

$ samtools view -bS bwa_results/Mb152.16S_M_bovis.paired.sam | samtools sort -o bwa_results/Mb152.16S_M_bovis.paired.bam
$ samtools view -bS bwa_results/Mb152.16S_M_bovis.unpaired.sam | samtools sort -o bwa_results/Mb152.16S_M_bovis.unpaired.bam

$ samtools merge bwa_results/Mb152.16S_M_bovis.bam bwa_results/Mb152.16S_M_bovis.paired.bam bwa_results/Mb152.16S_M_bovis.unpaired.bam
```

We can now download the final `bam` file and import it into `Geneious`.

> ### Exercise
>
> You have also obtained paired-end reads from two other *M. bovis* isolates. Write a loop to perform the `bowtie2` and `samtools` operations over these data.
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> $ module load Bowtie2/2.4.1-GCC-9.2.0
> $ module load SAMtools/1.12-GCC-9.2.0
>
> $ for i in Mb1 Mb168;
> do
>     # Map
>     bwa mem bwa_index/NZ_CP005933_16S_rRNA \
>             ../2_Quality_filtered_data/${i}_1_trimmed.miseq.fastq \
>             ../2_Quality_filtered_data/${i}_2_trimmed.miseq.fastq > bwa_results/${i}.16S_M_bovis.paired.sam
>     bwa mem bwa_index/NZ_CP005933_16S_rRNA ../2_Quality_filtered_data/${i}_unpaired_trimmed.miseq.fastq > bwa_results/${i}.16S_M_bovis.unpaired.sam
>
>     # Sort and merge
>     samtools view -bS bwa_results/${i}.16S_M_bovis.paired.sam | samtools sort -o bwa_results/${i}.16S_M_bovis.paired.bam
>     samtools view -bS bwa_results/${i}.16S_M_bovis.unpaired.sam | samtools sort -o bwa_results/${i}.16S_M_bovis.unpaired.bam
>     samtools merge bwa_results/${i}.16S_M_bovis.bam bwa_results/${i}.16S_M_bovis.paired.bam bwa_results/${i}.16S_M_bovis.unpaired.bam
>
>     # Tidy up
>     rm bwa_results/${i}.16S_M_bovis.{paired,unpaired}.{bam,sam}
> done
> ```
> </details>

---

[Next lesson](03-ont-mapping.md)
