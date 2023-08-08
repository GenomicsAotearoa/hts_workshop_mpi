# Polishing of Oxford Nanopore assemblies

* Teaching: 20 minutes
* Exercises: 40 minutes

#### Objectives

* Use `racon` and `medaka` to polish an assembly
* Compare the outputs of several assemblies using `QUAST`

#### Keypoints

* Rounds of refinement are sometimes required to reduce error and can improve the initial assembly
* Sometimes the above sentence does not apply and polishing does not really improve the assembly. The only way to know is trying!

---

## Preparing to polish our assembly

Once we have our draft assembly, we want to revisit it and attempt to improve regions by aligning the original reads against the assembled contigs and trying to improve regions which might have been modelled incorrectly by ther assembler. We are going to use two tools for performing this polishing process, the first is a tool called `racon` and the second is `medaka`.

Whether or not our assembly will benefit from polishing is hard to predict. When working through this process it is a good idea to make copies of your data as you perform different correction procedures (or combinations of procedures) and to evaluate each outcome.

In theory, good assemblies can be obtained by applying a workflow of `Canu` -> `racon` -> `medaka`, but this is may not be the case for your data.

Before proceding, we will create a directory for all polishing attempts and make a copy of the initial assembly which will form the basis of our polishing steps.

```bash
$ mkdir ont_assemblies/

$ cp Mb1_flye/assembly.fasta ont_assemblies/Mb1.flye.fna
```

---

## Performing an initial tidy up with `racon`

Strictly speaking, `racon` ([source](https://github.com/isovic/racon)) is designed for polishing assemblies which have been obtained from tools that do not perform extensive error correction themselves. However, we are going to apply it to the `Flye` assembly as we have limited time in this workshop. Anecdotally, I have not found any harm in applying `racon` to a `Flye` assembly, although this is not a guarantee.

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

When working with `medaka` it is important to note that the basecalling model used by `guppy` during sequencing must be provided as a parameter. Unfortunately, as our mock data were obtained from NCBI we do not know for certain which basecalling model was used. According to the project description ([PRJEB38523](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB38523)) for these data, parallel basecalling attempts were made with `guppy` (the data we have) and version 0.1.3 of the `bonito` basecaller. This version was released in [April 2020](https://pypi.org/project/ont-bonito/0.1.3/#history) so with a bit of detective work I *think* that the correct basecalling model for our data should be either **r941_min_high_g351** or **r941_min_high_g360**.

In practice we will always know which model was used during basecalling as we provide this parameter to `guppy` as part of basecalling. In addition, if we perform basecalling live, then `MinKNOW` tells us the model used in the run report. This is only an uncertainty in this training exercise as we do not have access to the raw sequencing files.

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

For the purposes of this exercise we will assume that the correct model is **r941_min_high_g360**. We can load `medaka` and execute it with the following commands:

```bash
$ module purge
$ module load medaka/1.4.3-Miniconda3-4.10.3

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

You will be able to find a copy of these assemblies in the directory `/nesi/project/nesi03181/phel/module_3/3_Assembly-mapping/`. Copy the references to your `ont_assemblies/` folder and run `QUAST` with the following command:


```bash
$ module purge
$ module load QUAST/5.0.2-gimkl-2018b

$ quast.py -r GCF_000696015.1.fna -o quast/ --gene-finding \
           ont_assemblies/Mb1.flye.fna \
           ont_assemblies/Mb1.racon.fna \
           ont_assemblies/Mb1.unicycler.fna \
           ont_assemblies/Mb1.canu.fna
```

>**Note:** You may need to change the path to your reference genome `GCF_000696015.1.fna`, depending on where you have downloaded it to.

We can now open the `quast/report.pdf` file in the Jupyter browser to inspect the results.

---
