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

---

## Mapping reads with `minimap2`

---

<Text to come>

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

The results can now be imported into `Geneious` just as before.

---
