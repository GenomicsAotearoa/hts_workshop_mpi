# *De novo* assembly of Oxford Nanopore sequences

* Teaching: 20 minutes
* Exercises: 40 minutes

#### Objectives

* Use `Flye` to assemble quality filtered Nanopore sequences
* Use `racon` and `medaka` to polish the assembly
* Compare the outputs of several assemblies using `QUAST`

#### Keypoints

* Assembling Oxford Nanopore sequences is a different process to assembling Illumina data
* Rounds of refinement are required to reduce error and improve the initial assembly

---

## Contents

1. [A word on assembly tools](#a-word-on-assembly-tools)
1. [Performing *de novo* assembly with `Flye`](#performing-de-novo-assembly-with-flye)
1. [Preparing to polish our assembly](#preparing-to-polish-our-assembly)
1. [Performing an initial tidy up with `racon`](#performing-an-initial-tidy-up-with-racon)
1. [Performing additional polishing with `medaka`](#performing-additional-polishing-with-medaka)
1. [Comparing outputs with `QUAST`](#comparing-outputs-with-QUAST)

---

## A word on assembly tools

For the exercise today we will be using the `Canu` assembly tool to work with one of the *M. bovis* genomes. Like with other areas of genomics, there are many good options for assembly tools and our usage of `Canu` today is in no way an endorsement that we consider this tool to be the 'best' long read assembler. `Canu` is a very good tool an will give us good results with the data we process today, but when working with real data there are many other good options to try, including:

1. `MaSuRCA` ([Zimin *et al.*, 2013](https://doi.org/10.1093/bioinformatics/btt476)) - [https://github.com/alekseyzimin/masurca](https://github.com/alekseyzimin/masurca)
1. `UniCycler` (and `TriCycler`) ([Wick *et al.*, 2017](https://doi.org/10.1371/journal.pcbi.1005595)) - [https://github.com/rrwick/Unicycler](https://github.com/rrwick/Unicycler)
1. `Canu` ([Koren *et al.*, 2017](http://www.genome.org/cgi/doi/10.1101/gr.215087.116))

A recent comparison of assembly tools was published by [Wick & Holt (2021)](https://doi.org/10.12688/f1000research.21782.4) which tests some of the options listed above along with several other tools.

In practice, there are sometimes particular cases where a tool will not be compatible with your data, so it is helpful to be aware of several tools so that you have options is assembly proves problematic for a particular sample.

---

## Performing *de novo* assembly with `Flye`

As ONT data are fundamentally more error prone than the sequences we obtained through Illumina sequencing, a considerable amount of an ONT assembly is spent identifying and correcting errors to produce high-quality contigs from a comparably low-quality set of reads. If you think back to the original presentation, recall this figure of the error rates of our Illumina (green) and MinION (purple) sequences:

![](../img/03_ont_vs_illumina_quality.png)

The median sequence quality for the MinION data sits around Q20 for most of the sequence. This corresponds to 99% accuracy which might sound good but by definition half of the sequences have lower quality than this. At the low end of this plot the sequences are slightly above Q10, which denotes 90% accuracy. Finding consensus regions between pairs of reads, when one of them might differ by up to 10% of it's composition **_just due to sequencing error alone_** makes assembly a complicated process and assembly tools which are aware of the error profiles of our long read data are essential.

The complete workflow of `Flye` is published ([Kolmogorov *et al.*, 2019](https://doi.org/10.1038/s41587-019-0072-8)) but it is quite complicated process. The novel aspect of assembly with `Flye` when compared with other asssembly tools was the developers observation that when working with noisy reads (as mentioned above) mapping sequences against each other is confounded by highly similar repeat regions, which create hotspots of local alignment between reads with very different flanking sites.

For our purposes, the main points of the assembly process to understand are:

1. Repetitive regions of the genome are identified
1. Repeat regions are clustered together to form deliberately misassembled **_disjointigs_**
1. Disjointigs are concatenated and a repeat graph (similar to a de Bruijn graph) is created
1. The input reads are mapped to the repeat graph
1. Using the reads that map to the repeat region and it's 5' and 3' flanking regions, the loops in the assembly graph are unwound

To run `Flye`, navigate to your `3_Assembly-mapping/` directory, and we will load the latest module from NeSI:

```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/3_Assembly-mapping/

$ module load Flye/2.8.3-gimkl-2020a-Python-3.8.2
```

You can explore the `Flye` documentation [online](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md) or through the help menu, but at a minimum to make an assembly work, you can use the following parameters:

```bash
$ flye --threads 10 --genome-size 1m \
       --nano-raw ../2_Quality_filtered_data/Mb1.trimmed.minion.fastq \
       --out-dir Mb1_flye/
```

---

## Preparing to polish our assembly

Once we have our draft assembly, we want to revisit it and attempt to improve regions by aligning the original reads against the assembled contigs and trying to improve regions which might have been modelled incorrectly by `Canu`. We are going to use two tools for performing this polishing process, the first is a tool called `racon` and the second is `medaka`.

Whether or not out assembly will benefit from polishing is hard to predict. When working through this process it is a good idea to make copies of your data as you perform different correction procedures (or combinations of procedures) and to evaluate each outcome.

In practice, I (David) have obtained good assemblies applying a workflow of `Canu` -> `racon` -> `medaka`, but this is may not be the case for your data.

Before proceding, we will create a directory for all polishing attempts and make a copy of the initial assembly which will form the basis of our polishing steps.

```bash
$ mkdir ont_assemblies/

$ cp Mb1_flye/Mb1.contigs.fasta ont_assemblies/Mb1.flye.fna
```

---

## Performing an initial tidy up with `racon`

Strictly speaking, `racon` ([source](https://github.com/isovic/racon)) is designed for polishing assemblies which have been obtained form tools that do not perform extensive error correction themselves. However, we are going to apply it to the `Canu` assembly as we have limited time in this workshop. Anecdotally, I have not found any harm in applying `racon` to a `Canu` assembly, although this is not a guarantee.

Before running `racon` we must produce a mapping file of the quality filtered sequences against the assembly. We can do this with `minimap2`:

```bash
$ module load minimap2/2.17-GCC-9.2.0

$ minimap2 -t 10 -ax map-ont ont_assemblies/Mb1.flye.fna ../2_Quality_filtered_data/Mb1.trimmed.minion.fastq > Mb1.sam
```

We can then use this mapping file as the input for `racon`:

```bash
$ module load Racon/1.4.13-GCC-9.2.0

$ racon -t 10 ../2_Quality_filtered_data/Mb1.trimmed.minion.fastq Mb1.sam ont_assemblies/Mb1.flye.fna > ont_assemblies/Mb1.racon.fna
```

>**Note:** It is possible to perform the `racon` process iteratively, remapping reads to the output and then running the polishing cycle again. There is some data ([link here](https://nanoporetech.github.io/medaka/draft_origin.html#discussion)) which suggests that up to four rounds of `racon` polishing, in conjunction with `medaka`, produces better quality output than running a single polishing step. However there are costs associated with this approach both in terms of time invested and over-zealous correction to repeat regions. Whether or not improvement with multiple rounds will be seen in your data is unclear, and ultimately it is your decision whether or not to perform this approach.

---

## Performing additional polishing with `medaka`

The `medaka` software is developed by Oxford Nanopore Technologies and claims to provide drastically improved assemblies when used in conjunction with `racon`. The [documentation](https://nanoporetech.github.io/medaka/) for `medaka` recommends using the `Flye` assembler for creating draft assemblies but they also report data showing that `Canu` produces the highest quality assembly.

When working with `medaka` it is important to note that the basecalling model used by `guppy` during sequencing must be provided as a parameter. Unfortunately, as our mock data were obtained from NCBI we do not know for certain which basecaling model was used. According to the project description ([PRJEB38523](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB38523)) for these data, parallel basecalling attempts were made with `guppy` (the data we have) and version 0.1.3 of the `bonito` basecaller. This version was released in [April 2020](https://pypi.org/project/ont-bonito/0.1.3/#history) so with a bit of detective work I *think* that the correct basecalling model for our data should be either **r941_min_high_g351** or **r941_min_high_g360**.

In practice we will always know which model was used during basecalling as we provide this parameter to `guppy` as part of basecalling. In addition, if we perform basecalling live, then `MinKNOW` tells us the model used in the run report. This is only an uncertainty in this training exercise as we do not have access to the raw `fast5` sequencing files.

> <details>
> <summary>Understanding `medaka` model names</summary>
>
> The model names might look cryptic, but they are actually extremely informative. The `medaka` [github page](https://github.com/nanoporetech/medaka) explains the naming convention used, which takes the form
> ```
> [PORE TYPE]_[SEQUENCING DEVICE]_[CALLING METHOD]_[VERSION]
> ```
> The pore type and sequencing device are dictated by the sequencing device we use - generally this would be MinION with either the R9.4.1 or R10.3 pore chemistry. If you are working on the sequencing computers look up the relationship between device, kit, and model using the command:
>
> ```bash
> $ guppy_basecaller --print_workflows
> ```
> This will show you which model corresponds to which sequencing kit. The version number is taken from the version of `guppy` used during basecalling.
>
> The only detail which is really in the users control is whether we used the fast, high-accuracy, or super-accuracy (`guppy` v5 only) models for basecalling. This is something that you should know, and will be in the run report if live basecalling was performed by `MinKNOW`.
> </details>

For the purposes of this exercise we will assume that the correct model is **r941_min_high_g360**. As `medaka` is not installed on NeSI, we have obtained a small software container which contains `medaka` and its dependencies. To access this container perform the following steps:

```bash
$ module purge
$ module load Miniconda3/4.10.3

$ conda activate /nesi/project/nesi03181/phel/module_3/envs/medaka
$ medaka -h
```

While this environment is active we can access `medaka` as if it were a native piece of software.

```bash
$ medaka_consensus -t 10 -m r941_min_high_g360 \
                   -i ../2_Quality_filtered_data/Mb1.trimmed.minion.fastq \
                   -d ont_assemblies/Mb1.racon.fna \
                   -o Mb1_medaka/

$ cp Mb1_medaka/consensus.fasta ont_assemblies/Mb1.medaka.fna
```

One important piece on information to note with the polishing process is that `racon` and `medaka` **_do not_** change the names of contigs during polishing. This is helpful, as it allows us to easily compare contigs between different polishing steps but it also means that you have to be careful when importing the data into `Geneious` as it might become hard to track which step of the analysis your contig comes from.

As an easy solution to this is to rename your contigs using `seqmagick` to append some versioning information to each sequence name:

```bash
$ module load seqmagick/0.7.0-gimkl-2018b-Python-3.7.3

$ seqmagick mogrify --name-suffix _flye ont_assemblies/Mb1.flye.fna
$ seqmagick mogrify --name-suffix _racon ont_assemblies/Mb1.racon.fna
$ seqmagick mogrify --name-suffix _medaka ont_assemblies/Mb1.medaka.fna
```

>**Note:** When running `seqmagick` you will get a warning which looks like:
>
>```bash
>/opt/nesi/CS400_centos7_bdw/Python/3.7.3-gimkl-2018b/lib/python3.7/importlib/_bootstrap.py:219: RuntimeWarning: This module has been deprecated. We encourage users to switch to alternative libraries implementing a trie data structure, for example pygtrie.
>    return f(*args, **kwds)
>/opt/nesi/CS400_centos7_bdw/Python/3.7.3-gimkl-2018b/lib/python3.7/site-packages/Bio/triefind.py:34: BiopythonDeprecationWarning: This module has been deprecated. We encourage users to switch to alternative libraries implementing a trie data structure, for example pygtrie.
>    "for example pygtrie.", BiopythonDeprecationWarning)
>```
>
>Don't worry about this. It is a warning to developers and does not affect our work.

---

## Comparing outputs with `QUAST`

To finish this exercise, we will compare out assembly produced using `Flye` to comparable assemblies produced with `Canu` and `Unicycler`.

You will be able to find a copy of these assemblies in the directory `/nesi/project/nesi03181/phel/module_3/3_Assembly-mapping/`. Copy the references to your `3_Assembly-mapping/ont_assemblies/` folder and run `QUAST` then compress the output with the following commands


```bash
$ module purge
$ module load QUAST/5.0.2-gimkl-2018b

$ quast.py -r GCF_000696015.1.fna -o quast/ --gene-finding \
           ont_assemblies/Mb1.flye.fna ont_assemblies/Mb1.racon.fna ont_assemblies/Mb1.unicycler.fna ont_assemblies/Mb1.canu.fna
```

>**Note:** You may need to change the path to your reference genome `GCF_000696015.1.fna`, depending on where you have downloaded it to.

We can now open the `quast/report.pdf` file in the Jupyter browser to inspect the results.

---

[Next lesson](08-hybrid-assembly.md)
