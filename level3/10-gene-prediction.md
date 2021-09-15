# Prediction of protein coding genes in an assembly

* Teaching: 30 minutes
* Exercises: 30 minutes

#### Objectives

* Understand the limitations of *de novo* gene prediction tools
* Learn to use the `prodigal` tool for predicting protein coding sequences in a bacterial genome
* Learn some basic use of a tool designed for eukaryotic gene prediction

#### Keypoints

* *De novo* protein sequence prediction is a good starting point, but proper annotation will require significant manual curation
* Not all tools are able to predict gene boundaries over splice junctions - be careful when interpreting predictions
* Protein sequence prediction does not predict all informative genetic elements - additional tools are required for a complete annotation
  * For example, the tools used in this workshop cannot identify rRNA, tRNA, or ncRNA elements

---

## Contents

1. [Complexities of prediction protein coding sequences from a *de novo* assembly](#complexities-of-prediction-protein-coding-sequences-from-a-de-novo-assembly)
1. [Prediciting prokaryote coding regions with `prodigal`](#prediciting-prokaryote-coding-regions-with-prodigal)
   1. [`prodigal` prediction parameters](#prodigal-prediction-parameters)
   1. [Running `prodigal`](#running-prodigal)

---

## Complexities of prediction protein coding sequences from a *de novo* assembly

Prediction of genes in a genome assembly is a complicated process - there are many tools which can perform good initial predictions from assembled contigs, but there are often many biological features which confound the prediction process and make it more complicated than simply finding start and stop codons within a sequence. Depending on the work that you are doing, you may be most interested in identifying protein coding regions or untranslated genetic elements such as ribosomal and transfer RNA sequences (rRNA and tRNAs, respectively). For features such as these, which are functional but not not translated a different set of tools is required for prediction and we will cover these in the next workshop.

At the most basic level, searching for proteins is simply looking for open reading frames (ORF) within a contig, but in practice there are many features (some biological, some as a consequence of having a draft assembly) which confound the process. At the biological level, features such as pseudogenes, the presence of alternate stop codons (the [amber, umber, and ochre codons](https://en.wikipedia.org/wiki/Stop_codon#Nomenclature) and stop codon read-through or intron splicing make obtaining the protein coding sequence more complicated than simply translating the nucleotide sequence between a start/stop pairing. In addition, if our genome assembly is not complete we run the risk of encountering partial coding sequeences in which either the 5' or 3' region of the sequence were not assembled. In these cases, a simple search for ORFs will fail to detect the partial sequence.

The prediction of protein coding sequences is therefore done using more complicated techniques. Similar to assembly, we can perform gene prediction in either a reference-guided manner, or through the use of *ab initio* prediction tools. We will not be covering the reference-guided approach here, as it is quite simple to perform in `Geneious`, but it is not to be underestimated as a technique, particularly when working with viruses or other organisms with complex read-through or splicing properties. *Ab initio* pprediction is akin to *de novo* assembly - the tool is created with some internal models for what coding regions 'look' like, which are then applied to query sequences to find putative coding regions. Depending on the intended use of the tool, each prediction tool may be better tuned for partiular assumptions of the data. We are going to use two different tools today, one designed for prediction of prokaryotic coding sequences (which generally lack introns) and eukaryotic sequences.

---

## Prediciting prokaryote coding regions with `prodigal`

Navigate to your working directory on NeSI. Create a folder `4_Gene_predictions/` and copy the contents of the `/nesi/project/nesi03181/phel/module_3/4_Gene_predictions/` folder into your folder.

> <details>
> <summary>Solution</summary>
> 
> ```bash
> $ cd /nesi/nobackup/nesi03181/phel/USERNAME/
>
> $ mkdir 4_Gene_predictions/
> $ cp ../../module_3/4_Gene_predictions/*.fna ./4_Gene_predictions/
> ```
> </details>

We are going to be predicting protein coding sequences for the *M. bovis* genome using the tool `prodigal`, produced by [Hyatt *et al.*](https://doi.org/10.1186/1471-2105-11-119). This is a powerful prediction tool which is quick to run, and flexible enough for most projects.

Use `module spider` or `module avail` to find the `prodigal` modules on NeSI and then load one. As the modules are both the same version of the tool it will not matter which one you load for this exercise. Once you have loaded `prodigal`, run it with the `-h` parameter to view the help menu.

```bash
$ module purge
$ module spider prodigal

$ module load prodigal/2.6.3-GCCcore-7.4.0

$ prodigal -h
```

### `prodigal` prediction parameters

There are quite a few options given here, which we can split into input and output parameters. From the input parameters, the ones we need to be most aware of are the following.

|Parameter|Function|
|:---:|:---|
|`-i`|Specify FASTA/Genbank input file (default reads from stdin).|
|`-p`|Select procedure (single or meta).  Default is single.|
|`-g`|Specify a translation table to use (default 11).|
|`-c`|Closed ends.  Do not allow genes to run off edges.|
|`-q`|Run quietly (suppress normal `stderr` output).|

The first parameter here should be obvious, so will not be discussed. The `-p` parameter is an important one to pay attention to, as it determines whether or not we are predicting sequences in an single genome (i.e. a genome obtained from a pure culture) or a metagenomic mix of sequences. In `single` mode prediction goes through a round of training against the input sequences before producing a prediction output tailored to your contigs. The `meta` mode uses pre-calculated profiles of gene features. Strictly speaking it is possible to run `meta` mode on an isolate genome, or `single` mode on a metagenome but the results do differ.

Selecting the translation table (`-g`) is something we usually do not need to wory about. `prodigal` supports nubmers 1 through 25 of the [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). By default, it will start with code 11 (bacteria, archaea, and plastids) but shift to code 4 if the predictions are too short. This is important to us because genetic code 4 corresponds to the *Mycoplasmataceae* - the bacterial family that contains the genera *Mycoplasma* and *Ureaplasma*, and the family *Spiroplasmataceae*. These lineages have repurposed `UGA` from a stop codon to a tryptophan codon.

The `-c` parameter is one we generally do not need to use. As we will almost always be working with draft genomes, there is no problem with allowing genes to run over the edges of contigs, as we know that the contigs are not completely assembled.

The quiet (`-q`) parameter is simply a quality of life feature, which reduces the feedback provided by `prodigal` during prediction.

The output parameters are more simple:

|Parameter|Function|
|:---:|:---|
|`-a`|Write protein translations to the selected file.|
|`-d`|Write nucleotide sequences of genes to the selected file.|
|`-o`|Specify output file (default writes to `stdout`).|
|`-f`|Select output format. Default is gbk.|

Three of these capture the prediction information in commonly used formats. The output of `-d` and `-a` are simply fasta format, and the `-o` (or `stdout`) output are in GenBank format. At the very least, we should be capturing either `-d` **and** `-a`, or `-o` to get the full prediction.

There are also some features for producing customised prediction models for a lineage of interest. The full description of all parameters can be found [here](https://github.com/hyattpd/prodigal/wiki) if you want more information. Note that the latest version of `prodigal` is not available on NeSI so some of the parameter names differ between the version we are using and the latest.

### Running `prodigal`

One of the nice features of `prodigal` is that it does not take a lot of resources to run, so we can easily run it without resorting to `slurm` for a single genome.

```bash
$ prodigal -i M_bovis.NZ_CP005933.fna -p single -g 4 \
           -d M_bovis.NZ_CP005933.prod.fna -a M_bovis.NZ_CP005933.prod.faa -o M_bovis.NZ_CP005933.prod.gbk
```

This will only take a few seconds to complete. Once done, we can quickly see how many protein coding sequences were predicted by counting the number of sequences in either of the output fasta files.

> <details>
> <summary>Solution</summary>
> 
> ```bash
> $ grep -c ">" M_bovis.NZ_CP005933.prod.faa
> ```
> </details>

Since this is a chromosome downloaded from the [NCBI RefSeq database](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP005933) we can look to the official annotation and see how many proteins we should expect to find. This number (782) is slightly fewer than what `prodigal` predicts (792) but is very close.

If we inspect the output of the fasta files, you will notice that there is a lot of metadata stored in each line. For example:

```bash
$ head -n1 M_bovis.NZ_CP005933.prod.faa
# >NZ_CP005933.1_1 # 1 # 1401 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.293
```

There is quite a lot to unpack here - some is quite important to know and other parts are just descriptive. We can break down the results like so:

```
>NZ_CP005933.1_1 # 1 # 1401 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.293
 |             |   |   |      |   |      |
 |             |   |   |      |   |      Are the gene boundaries complete or not
 |             |   |   |      |   Unique gene ID
 |             |   |   |      Orientation
 |             |   |   Stop position
 |             |   Start position
 |             Gene suffix
 Contig name
```

As `prodigal` is only concerned with *predicting* sequences, not annotating them, there is no functional imformation which can be used to guide the name of each prediction. For simplicity, protein predictions are simply named as `[CONTIG NAME]_[PREDICTION]`. This has some nice implications when we want to study gene synteny but for the most part numbering off the predicitons in the order they are made is the simpliest way to generate names.

The next three numbers provide the nucleotide coordinates of the coding sequence and the orientation (1 for forward, -1 for reverse). Again, these are very useful when tracking down the genomic context of a sequence and can sometimes provide a quick-and-dirty means for spotting rearrangements.

The unique gene identifier is used to link the entries in the fasta file to the output obtained from the `-o` (or `stdout`) channel. Simiarly to the prediction number, it is simply derived from the order of contigs and sequence of predictions.

Finally, the last piece of information we are going to inspect is whether or not the prediction is complete of not.








---
