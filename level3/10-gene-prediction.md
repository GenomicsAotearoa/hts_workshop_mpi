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

1. [Complexities of predicting protein coding sequences from a *de novo* assembly](#complexities-of-predicting-protein-coding-sequences-from-a-de-novo-assembly)
1. [Predicting prokaryote coding regions with `prodigal`](#predicting-prokaryote-coding-regions-with-prodigal)
   1. [`prodigal` prediction parameters](#prodigal-prediction-parameters)
   1. [Running `prodigal`](#running-prodigal)
1. [Predicting eukaryote coding regions with `AUGUSTUS`](#predicting-eukaryote-coding-regions-with-augustus)

---

## Complexities of predicting protein coding sequences from a *de novo* assembly

Prediction of genes in a genome assembly is a complicated process - there are many tools which can perform good initial predictions from assembled contigs, but there are often many biological features which confound the prediction process and make it more complicated than simply finding start and stop codons within a sequence. Depending on the work that you are doing, you may be most interested in identifying protein coding regions or untranslated genetic elements such as ribosomal and transfer RNA sequences (rRNA and tRNAs, respectively). For features such as these, which are functional but not not translated, a different set of tools is required for prediction and we will cover these in the next workshop.

At the most basic level, searching for proteins is simply looking for open reading frames (ORF) within a contig, but in practice there are many features (some biological, some as a consequence of having a draft assembly) which confound the process. At the biological level a number of features complicate the process of predicting protein coding regions:

1. Alternate coding schemes, including the [amber, umber, and ochre stop codons](https://en.wikipedia.org/wiki/Stop_codon#Nomenclature)
1. Stop codon read-through
1. Splicing of intronic regions

As a consequences, simply translating the nucleotide sequence between a start/stop pairing of not sufficient to correctly identify the protein complement of the genome. In addition, if our genome assembly is not complete we run the risk of encountering partial coding sequeences in which either the 5' or 3' region of the sequence were not assembled. In these cases, a simple search for ORFs will fail to detect the partial sequence. The prediction of protein coding sequences is therefore done using more complicated techniques than a simple `grep` search for the start and stop codon(s).

Similar to assembly, we can perform gene prediction in either a reference-guided manner, or through the use of *ab initio* prediction tools. We will not be covering the reference-guided approach here, as it is quite simple to perform in `Geneious`, but it is not to be underestimated as a technique, particularly when working with viruses or other organisms with complex read-through or splicing properties.

*Ab initio* prediction is akin to *de novo* assembly - the tool is created with some internal models for what coding regions 'look' like, which are then applied to query sequences to find putative coding regions. Depending on the intended use of the tool, each prediction tool may be better tuned for partiular assumptions of the data. We are going to use two different tools today, one designed for prediction of prokaryotic coding sequences (which generally lack introns) and one designed primarily for eukaryotic sequences, where splicing is common.

---

## Predicting prokaryote coding regions with `prodigal`

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

We are going to be predicting protein coding sequences for the *M. bovis* genome using the tool `prodigal` ([Hyatt *et al.*, 2010](https://doi.org/10.1186/1471-2105-11-119)). This is a powerful prediction tool which is quick to run, and flexible enough for most projects.

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

Selecting the translation table (`-g`) is something we usually do not need to wory about. `prodigal` supports numbers 1 through 25 of the [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). By default, it will start with code 11 (bacteria, archaea, and plastids) but shift to code 4 if the predictions are too short. This is important to us because genetic code 4 corresponds to the *Mycoplasmataceae* - the bacterial family that contains the genera *Mycoplasma* and *Ureaplasma*, and the family *Spiroplasmataceae*. These lineages have repurposed `UGA` from a stop codon to a tryptophan codon.

The `-c` parameter is one we generally do not need to use. As we will almost always be working with draft genomes, there is no problem with allowing genes to run over the edges of contigs, as we know that the contigs are not completely assembled.

The quiet (`-q`) parameter is simply a quality of life feature, which reduces the feedback provided by `prodigal` during prediction.

The output parameters are more simple:

|Parameter|Function|
|:---:|:---|
|`-a`|Write protein translations to the selected file.|
|`-d`|Write nucleotide sequences of genes to the selected file.|
|`-o`|Specify output file (default writes to `stdout`).|
|`-f`|Select output format. Default is gbk.|

Three of these capture the prediction information in commonly used formats. The output of `-d` and `-a` are simply fasta format, and the `-o` (or `stdout`) output are in GenBank format. At the very least, we should be capturing either `-d` **_and_** `-a`, or `-o` to get the full prediction. It is helpful to have both the nucleotide and protein sequence of a coding region, as the different annotation techniques can be applied to each one, and certain methods of phylogenetic analysis favour one sequence type over the other.

>**Note:** One nice feature of `prodigal`, in addition to it running quickly, is that it is a fully deterministic algorithm. If you run it multiple times on the same input with the same settings, you will always get the exact same output. This means you can always go back and reproduce the analysis with additional output files if required,

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

After these parameters are done, there are a series of keyword-linked pieces of information about the sequence. The unique gene identifier is used to link the entries in the fasta file to the output obtained from the `-o` (or `stdout`) channel. Similar to the prediction number, it is simply derived from the order of contigs and sequence of predictions.

The main piece of information we are going to inspect is whether or not the prediction is complete of not. The 'partial' keyword provides a two digit code that reports the status of the prediction, the first digit corresponds to the start of the sequence and the second to the end of the sequence. The values these can take are either 0 (complete) or 1 (incomplete). A fully complete prediction will therefore have a code of `00`, and a partial predictions can be:

* `01` - Started, but no end found - typically when a prediction runs off the end of a contig
* `10` - No start identified, but a complete end found - often a sequence which occurs at the start of the contig
* `11` - No ends found - likely due to predicting for a very short contig

You can also inspect the sequences to see if they appear complete. Typically protein predictions begin with a methionine (M) amino acid residue, as this is the translation of the `ATG` codon. There is no residue which corresponds to stop (as stop is a gap in translation) so `prodigal` reports the stop position with an asterisk (`*`).

---

## Predicting eukaryote coding regions with `AUGUSTUS`

Unlike prokaryotic genomes, the genes of eukaryotes carry intronic sequences, which need to be spliced out of the gene sequence before undergoing translation. The detection of splicing boundaries is a difficult task, as there are many organism-specific patterns used to mark splice sites. Protein prediction tools for this purpose typically come with  number of pre-trained models for finding protein domains within contigs, but if there is no model for your organism, or a closely related lineage, then results may not be ideal.

A recently published article ([Scalzitti *et al.*, 2020](https://doi.org/10.1186/s12864-020-6707-9)) which profiled a number of these tools found `AUGUSTUS` to be one of the best performing tools for gene prediction in eukaryotic organisms, so this is what we will use today.

>**Note:** `AUGUSUTUS` does require training against a closely related model organisms to generate accurate predictions, which we do not have for this workshop. We will instead be performing predictions with a few different models and seeing how the outputs differ.

To begin, find the latest version of `AUGUSTUS` on NeSI and load it.

> <details>
> <summary>Solution</summary>
> 
> ```bash
> $ module spider augustus
> $ module load AUGUSTUS/3.3.3-gimkl-2020a
> ```
> </details>

Once loaded, you can view the available models for gene prediction using the following command:

```bash
$ augustus --species=help
```

You will see a number of results, including some from the insecta, but nothing particularly closely related to the *Pentatomidae*. Instead, we will use two different models for an initial round of prediction on the BMSB mitochondrial genome - one insect and one bacterial species.

```bash
$ augustus --protein=on --codingseq=on --species=honeybee1 H_halys.ML746646.fna > H_halys.ML746646.aug_hb.gff
$ getAnnoFasta.pl H_halys.ML746646.aug_hb.gff

$ augustus --protein=on --codingseq=on --species=s_aureus H_halys.ML746646.fna > H_halys.ML746646.aug_sa.gff
$ getAnnoFasta.pl H_halys.ML746646.aug_sa.gff
```

If you count the number of coding sequences in each prediction, you will notice that they are different. This is because neither of these gene models are accurate for the organism we are trying to characterise. It is possible to create a custom species profile for an organism to get a more accurate prediction. Creating a new model is a slow process so we will not cover it here, although you can see the steps required in the text below.

> <details>
> <summary>Creating a custom gene profile</summary>
> 
> There are a few steps we need to perform in advance of the new model training. The first is to do with file permissions - the location of the prediction databases that `AUGUSTUS` uses for gene prediction are not writable to us, so we cannot add new data into them. We must create our own copy of the configuration information and point `AUGUSTUS` towards this new location in order to create new models
>
> ```bash
> $ ECHO $AUGUSTUS_CONFIG_PATH
> #/opt/nesi/CS400_centos7_bdw/AUGUSTUS/3.3.3-gimkl-2020a/config
>
> # Create a local copy
> $ cp -r /opt/nesi/CS400_centos7_bdw/AUGUSTUS/3.3.3-gimkl-2020a/config/ /nesi/project/nesi03181/phel/USERNAME/config/
>
> # Update the `AUGUSTUS_CONFIG_PATH` variable to point to the new location
> $ AUGUSTUS_CONFIG_PATH="/nesi/project/nesi03181/phel/USERNAME/config"
> ```
>
> Now we just need to obtain a reference genome to train against. For *H. halys*, this can found on the [NCBI website](https://www.ncbi.nlm.nih.gov/assembly/GCA_000696795.3) and downloaded from the command line:
>
> ```bash
> $ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/696/795/GCA_000696795.3_Hhal_1.1/GCA_000696795.3_Hhal_1.1_genomic.gbff.gz
> $ gunzip GCA_000696795.3_Hhal_1.1_genomic.gbff.gz
> ```
>
> Once these steps are completed, training is a single command:
>
> ```bash
> $ autoAugTrain.pl --trainingset=GCA_000696795.3_Hhal_1.1_genomic.gbff --species=hhalys
> ```
>
> We would then perform prediction using the new model as usual. Note that if this was performed in a new session, we would need to set the `AUGUSTUS_CONFIG_PATH` variable again.
>
> ```bash
> $ AUGUSTUS_CONFIG_PATH="/nesi/project/nesi03181/phel/module_3/4_Gene_predictions/config"
>
> $ augustus --protein=on --codingseq=on --species=hhalys H_halys.ML746646.fna > H_halys.ML746646.aug_hh.gff
> $ getAnnoFasta.pl H_halys.ML746646.aug_hh.gff
> ```
> </details>

We will take a quick look at the outputs in `Geneious`. Download the output fasta files, ready to analyse.

---

[Next lesson](12-annotation.md)
