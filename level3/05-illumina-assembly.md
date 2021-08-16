# *De novo* genome assembly of short sequence reads

* Teaching: 30 minutes
* Exercises: 30 minutes

---

## Contents

1. [Genome assembly and *de Bruijn* graphs](#genome-assembly-and-de-bruijn-graphs)
1. [Introduction to the `SPAdes` assembly tool](#introduction-to-the-spades-assembly-tool)
1. [Performing an assembly using `SPAdes`](#performing-an-assembly-using-spades)
1. [Assessing the results of the assembly](#assessing-the-results-of-the-assembly)

---

## Genome assembly and *de Bruijn* graphs

Each sequence read can be disseted into *k*-mers of the length *k*. For example, for the sequence `ATGCAT`, where *k*=3, the sequence (string) can be dissected into four overlapping substrings of *3*-mers: 
<br>`ATG` 
<br>`TGC` 
<br>`GCA` 
<br>`CAT`

Each *k*-mer can be further split into two substrings: **prefix** and **suffix**. Both substrings have the length of *k-1*. For example, the *3*-mer `ATG` can be split into **prefix** `AT` and **suffix** `TG`. *k*-mers that shared a common **prefix** and **suffix** can be connected. For example, the **sufix** of `ATG` and the **prefix** of `TGC` are identical, therefore they can be connected and form a new sequence `ATGC`. In this way, more overlapping *k*-mers can be connected and forming a longer new sequence. The linkage of *k*-mers along these connections is known as a [directed acyclic graph](https://en.wikipedia.org/wiki/Directed_acyclic_graph). When many *k*-mers are connected in this way, they form paths which represent the original sequence of the input nucleic acid. 

However, due to the fact that short nucleic acid *k*-mers may not be unique in a genome, some *k*-mers will form multiple paths which may form loops in the graph - this structure is called a [de Bruijn graph](https://en.wikipedia.org/wiki/De_Bruijn_graph) and these structures in the graph must be resolved if we wish to correctly resolve our genome.

![](../img/03_debruijn_graph.png)

**Fig. 1.** Standard and multisized *de Bruijn graph*. A circular Genome CATCAGATAGGA is covered by a set Reads consisting of nine 4-mers, {ACAT, CATC, ATCA, TCAG, CAGA, AGAT, GATA, TAGG, GGAC}. Three out of 12 possible 4-mers from Genome are missing from Reads (namely {ATAG,AGGA,GACA}), but all 3-mers from Genome are present in Reads. (A) The outside circle shows a separate black edge for each 3-mer from Reads. Dotted red lines indicate vertices that will be glued. The inner circle shows the result of applying some of the glues. (B) The graph DB(Reads, 3) resulting from all the glues is tangled. The three h-paths of length 2 in this graph (shown in blue) correspond to h-reads ATAG, AGGA, and GACA. Thus Reads3,4 contains all 4-mers from Genome. (C) The outside circle shows a separate edge for each of the nine 4-mer reads. The next inner circle shows the graph DB(Reads, 4), and the innermost circle represents the Genome. The graph DB(Reads, 4) is fragmented into 3 connected components. (D) The multisized *de Bruijn* graph DB(Reads, 3, 4). From [Bankevich et. al 2012](https://dx.doi.org/10.1089%2Fcmb.2012.0021).
 
---

## Introduction to the `SPAdes` assembly tool

The `SPAdes` assembler is a very powerful tool for assembling genomes from short read sequence data using de Bruijn graphs to join *k*-mers into longer contiguous sequences (contigs). The tool was originally developed as an assembler for processing single-cell genomic data but it has since expanded to cover many more areas of biology (see the next module for more information). The `SPAdes` assembler is designed for the following cases:

* *De novo* assembly of short read sequences (i.e. Illumina or IonTorrent)
* *De novo* assembly of short reads and long reads (i.e. PacBio or Oxford Nanopore Technologies)
* Small size genomes (ideally <100 Mb, such as bacterial, viral, fungal, mitochondrial genomes)
  * Despite this, it is possible to apply `SPAdes` to larger genomes and obtain very good results.

Before beginning to work with `SPAdes` we need to ensure that our data is free from adapter sequences. The main reason for doing this is when we assemble seuqences to form a genome, the assembler is looking for spans of nucleic acid sequence which are common to multiple reads, so that those reads can be joined together to create longer contigs. As the adapater sequence is an identical tag added to every read, these create regions of artificial homology between sequences which have no real connection to each other. Before we attempt assembly it is critical to remove these from our data.

In order to begin, we must first find the versions of `SPAdes` installed on NeSI and load the module of interest.

```bash
$ module avail spades
$ module load SPAdes/3.15.2-gimkl-2020a

$ spades.py -h
```

---

## Performing an assembly using `SPAdes`

The test data is a set of Illumina MiSeq sequencing reads from a BMSB sample. Morphological identification indicated that the specimen was a brown marmorated stink bug (BMSB).

To save time, the reads have already been quality filtered with `fastp` and the number of reads reduced to speed up analysis.

```bash
# Copy the BMSB reads to your SPAdes working directory
$ cd /nesi/project/nesi03181/phel/USERNAME/2_Quality_filtered_data/
$ cp /nesi/project/nesi03181/phel/module_3/2_Quality_filtered_data/SRR13090255_*.fq.gz ./

$ cd ../
```

Create a `slurm` script with the following contents. Be sure to replace the `YOUR_EMAIL` and `USERNAME` values with your details.

```bash
#!/bin/bash -e
#SBATCH --account       nesi03181
#SBATCH --job-name      bmsb_spades
#SBATCH --time          00:10:00
#SBATCH --cpus-per-task 10
#SBATCH --mem           4G
#SBATCH --error         bmsb_spades.%j.err
#SBATCH --output        bmsb_spades.%j.out
#SBATCH --mail-type     END
#SBATCH --mail-user     YOUR_EMAIL

module purge
module load  SPAdes/3.15.2-gimkl-2020a

# Set the path from which the script will execute SPAdes
cd /nesi/project/nesi03181/phel/USERNAME/

# Execute SPAdes
spades.py --thread 10 \
          -1 2_Quality_filtered_data/SRR13090255_1.fq.gz \
          -2 2_Quality_filtered_data/SRR13090255_2.fq.gz \
          -s 2_Quality_filtered_data/SRR13090255_s.fq.gz \
          -o 3_Assembly-mapping/bmsb_spades/ 
```
Submit this job to `slurm`:

```bash
$ sbatch bmsb_spades.sl
```

---

## Assessing the results of the assembly

We will use the tool `QUAST` [source](http://bioinf.spbau.ru/quast) to obtain assembly statistics for the assembly. `QUAST` works better with a reference genome for comparison, so we will download one from the NCBI website using the `efetch` tool on NeSI.

We will also perform some filtering on the `SPAdes` assembly to remove some of the very short contigs.

```bash
$ module load entrez-direct/13.3
$ module load seqmagick/0.7.0-gimkl-2018b-Python-3.7.3

# Obtain the reference sequence for comparison
$ efetch -format fasta -db sequences -id NC_013272.1 > NC_013272.fna

# Filter out contigs shorter than 1 kbp in length
$ seqmagick convert --min-length 1000 bmsb_spades/contigs.fasta BMSB_mitochondria.fna
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

We can see how much data was filtered from the assembly file with a quick `grep` command:

```bash
$ grep -c ">" bmsb_spades/contigs.fasta BMSB_mitochondria.fna 
bmsb_spades/contigs.fasta:490
BMSB_mitochondria.fna:5
```

We can now use `QUAST` to compare the results of our filtered assembly with the BMSB reference mitochondria genome.

```bash
$ module purge
$ module load QUAST/5.0.2-gimkl-2018b

$ quast.py -r NC_013272.fna --gene-finding -o bmsb_quast/ BMSB_mitochondria.fna
```

Open the resulting `bmsb_quast/report.pdf` file in Jupyter using the file browser. How complete does the assembly appear to be, compared with the reference?

We can also visualise the assembly by looking at how well the loops and fragments of the assembly graph were resolved using the tool `Bandage` ([source](https://rrwick.github.io/Bandage/)). It is a fairly easy process to apply this to our assembly - the tool only requires a single output file from the `SPAdes` assembly as its input

```bash
$ module load Bandage/0.8.1_Centos
$ module load Qt5/5.13.2-GCCcore-9.2.0

$ Bandage image bmsb_spades/assembly_graph.fastg bmsb_bandage.svg 2> /dev/null
```

You can then open the `bmsb_bandage.svg` file in the Jupyter browser. Unfortunately, we cannot filter out the short contigs from this result in the same way we can filter the assembled contigs file. However, it should be clear that there is one long contig which has been assembled, and then a large number of short fragments.

---

[Next lesson](06-assembly-choices.md)
