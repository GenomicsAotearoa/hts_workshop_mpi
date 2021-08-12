# *De novo* assembly of Oxford Nanopore sequences

* Teaching: 20 minutes
* Exercises: 20 minutes

#### Objectives

* Use `Canu` to assemble quality filtered Nanopore sequences
* Use `racon` and `medaka` to polish the assembly

#### Keypoints

* Assembling Oxford Nanopore sequences is a different process to assembling Illumina data
* Rounds of refinement are required to reduce error and improve the initial assembly

---

## Contents

1. [Performing *de novo* assembly with `Canu`](#performing-de-novo-assembly-with-Canu)
1. [Preparing to polish our assembly](#preparing-to-polish-our-assembly)
1. [Performing an initial tidy up with `racon`](#performing-an-initial-tidy-up-with-racon)
1. [Performing additional polishing with `medaka`](#performing-additional-polishing-with-medaka)

---

## Performing *de novo* assembly with `Canu`

As ONT data are fundamentally more error prone than the sequences we obtained through Illumina sequencing, a considerable amount of an ONT assembly is spent identifying and correcting errors to produce high-quality contigs from a comparably low-quality set of reads. If you think back to the original presentation, recall this figure of the error rates of our Illumina (green) and MinION (purple) sequences:

![](../img/03_ont_vs_illumina_quality.png)

The median sequence quality for the MinION data sits around Q20 for most of the sequence. This corresponds to 99% accuracy which might sound good but by definition half of the sequences have lower quality than this. At the low end of this plot the sequences are slightly above Q10, which denotes 90% accuracy. Finding consensus regions between pairs of reads, when one of them might differ by up to 10% of it's composition **_just due to sequencing error alone_** makes assembly a complicated process and assembly tools which are aware of the error profiles of our long read data are essential.

The complete workflow of `Canu` is published ([Koren *et al.*, 2017](http://www.genome.org/cgi/doi/10.1101/gr.215087.116)) but it is quite complicated process which involved multiple rounds of aligning reads and attempting to resolve errors.

For our purposes, the main points to understand are:

1. Read correction
   1. Reads are mapped against each other to find regions of alignment common to multiple sequences.
   1. High quality sequences are used to correct ambiguous regions of lower quality sequences, but the number of times a sequence can be used to correct others is restricted.
1. Overlap-based trimming and error adjustment
   1. Corrected reads are trimmed to only high coverage regions, and short and low quality sequences are removed.
   1. Overlaps are re-computed and overlap error rates are computed
1. Graph construction
   1. An assembly graph is constructed from the remaining reads. *Dovetail* overlaps, where a pair of reads overlap at their ends, are identified and the 'best' overlaps for each read are identified.
   1. Filtering is performed to remove spurious overlaps, and the resulting set of reads are used to construct an initial set of contigs.
1. Contig consensus
   1. For each initial contig, a template sequence is created by splicing together the individual reads which comprise the contig.
   1. Reads are aligned to the template sequence and a final pass at indel correction is attempted.

Generally speaking, we do not need to know the specifics of how these steps are performed but they are important to know as `Canu` uses a set of default values for each of these processes. When assemblies fail, or do not produce a good result, tuning these values may drastically improve the final set of contigs as different alignment sensitivities or overlap sizes are considered.

To run `Canu`, navigate to your `3_Assembly-mapping/` directory, and we will load the latest module from NeSI:

```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/3_Assembly-mapping/

$ module load Canu/2.1.1-GCC-9.2.0
```

You can explore the `Canu` documentation [online](https://canu.readthedocs.io/en/latest/index.html) or through the help menu, but at a minimum to make an assembly work, you can use the following parameters:

```bash
$ canu -d Mb152_canu/ -p Mb152 \
       genomeSize=1m useGrid=false \
       -nanopore ../2_Quality_filtered_data/Mb152.trimmed.minion.fastq
```

There are several parameters here which we are using, so we will go through them briefly:

|Parameter|Meaning|
|:---|:---|
|`-d`|Output directory for assembly files|
|`-p`|Output file prefix (i.e. file names)|
|`genomeSize`|The approximate genome size for our organism. This is used to calculate the initial read coverage.|
|`useGrid`|By default, `Canu` will detect that it is running on a compute cluster with `slurm` installed and will dispatch itself in a `slurm` script. We don't want this behaviour.|
|`-nanopore`|Use pre-computed error rates representative of ONT data|


---

## Preparing to polish our assembly

Once we have our draft assembly, we want to revisit it and attempt to improve regions by aligning the original reads against the assembled contigs and trying to correct or resolve regions which might have been modelled incorrectly by `Canu`. We are going to use two tools for performing this polishing process, the first is a tool called `racon` and the second is `medaka`.

Whether or not out assembly will benefit from polishing is hard to predict. When working through this process it is a good idea to make copies of your data as you perform different correction procedures (or combinations of procedures) and to evaluate each outcome.

In practice, I (David) have had good experience performing a workflow of `Canu` -> `racon` -> `medaka` but this is may not be the case for your data.

Before proceding, we will create a directory for all polishing attempts and make a copy of the initial assembly which will form the basis of our polishing steps.

```bash
$ mkdir ont_polishing/
$ cp asdasdasd/asdasdasd ont_polishing/adasdasd
```

---

## Performing an initial tidy up with `racon`

Strictly speaking, `racon` ([source](https://github.com/isovic/racon)) is designed for polishing assemblies which have been obtained form tools that do not perform extensive error correction themselves. However, we are going to apply it to the `Canu` assembly as we have limited time in this workshop. Anecdotally, I have generally not found any harm in applying `racon` to a `Canu` assembly although this is not a guarantee.

Before running `racon` we must produce a mapping file of the quality filtered sequences against the assembly. We can do this with `minimap2`:

```bash
$ module load minimap2/2.17-GCC-9.2.0

$ minimap2 -ax map-ont ${TARGET} ../2_Quality_filtered_data/Mb152.trimmed.minion.fastq > ${TARGET}.sam
```

We can then use this mapping file as the input for `racon`:

```bash
$ module load Racon/1.4.13-GCC-9.2.0

$ racon ../2_Quality_filtered_data/Mb152.trimmed.minion.fastq ${TARGET}.sam ${TARGET} > ont_polishing/${TARGET}.racon.fna
```

**Note:** It is possible to perform the `racon` process iteratively, remapping the reads to each output and then re-running `racon` on the new alignment. There is some data in the `medaka` documentation ([link here](https://nanoporetech.github.io/medaka/draft_origin.html#discussion)) which suggests that up to four rounds of `racon` polishing, in conjunction with `medaka`, produces better quality output than running tools individually. However there are costs associated with this approach both in terms of time invested and over-zealous correction to repeat regions. Whether or not improvement with multiple rounds will be seen in your is unclear, and ultimately it is your decision whether or not to perform this approach.

---

## Performing additional polishing with `medaka`

The `medaka` software is developed by Oxford Nanopore Technologies and claims to provide drastically improved assemblies when used in conjunction with `racon`. The [documentation](https://nanoporetech.github.io/medaka/) for `medaka` recommends using the `Flye` assembler for creating draft assemblies but they also report data showing that `Canu` produces the highest quality assembly.

When working with `medaka` it is important to note that the basecalling model used by `guppy` during sequencing must be known for polishing to perform correctly. Unfortunately, as our mock data were obtained from NCBI we do not know for certain which basecaling model was used. According to the project description ([PRJEB38523](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB38523)) for these data, parallel basecalling attempts were made with `guppy` (the data we have) and `bonito v0.1.3`. This version of `bonito` was released in [April 2020](https://pypi.org/project/ont-bonito/0.1.3/#history) so with a bit of detective work I *think* that the correct basecalling model for our data will be either `r941_min_high_g351` or `r941_min_high_g360`.

For the purposes of this exercise we will assume XXXXX, but this may yield sub-optimal results. As `medaka` is not installed on NeSI, we have created a small software environment which contains `medaka` and its dependencies. To access this environment perform the following steps:

```bash
$ module purge
$ module load Miniconda3/4.10.3

$ conda activate -p /nesi/project/nesi03181/phel/module_3/envs/medaka
$ medaka -h
```

While this environment is active we can access `medaka` as if it were a native piece of software.
