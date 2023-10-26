# 4.1 - Prediction of protein coding sequences

!!! clock "time"

    * Teaching: 10 minutes
    * Exercises: 30 minutes

!!! circle-info "Objectives and Key points"

    #### Objectives
    
    * Understand the limitations of protein coding predictions when working with annotations

    #### Keypoints

    * *De novo* protein sequence prediction is a good starting point, but proper annotation will require significant manual curation
    * Not all tools are able to predict gene boundaries over splice junctions - be careful when interpreting predictions
    * Protein sequence prediction does not predict all informative genetic elements - additional tools are required for a complete annotation
      * For example, the tools used in this workshop cannot identify rRNA, tRNA, or ncRNA elements

---

## Complexities of predicting protein coding sequences from a *de novo* assembly

Prediction of genes in a genome assembly is a complicated process - there are many tools which can perform good initial predictions from assembled contigs, but there are often many biological features which confound the prediction process and make it more complicated than simply finding start and stop codons within a sequence.

At the most basic level, searching for proteins is simply looking for open reading frames (ORF) within a contig, but in practice there are many factors which confound the process. At the biological level a number of features complicate the process of predicting protein coding regions:

1. Alternate coding schemes, including the [amber, umber, and ochre stop codons](https://en.wikipedia.org/wiki/Stop_codon#Nomenclature)
1. Stop codon read-through
1. Splicing of intronic regions

Simply translating the nucleotide sequence between a start/stop pairing is not sufficient to correctly identify the complete protein complement of the genome.

Furthermore, if our genome assembly is not complete we run the risk of encountering partial coding sequeences in which either the 5' or 3' region of the sequence were not assembled. In these cases, a simple search for ORFs will fail to detect the partial sequence. The prediction of protein coding sequences must be achieved using more complicated techniques than a simple `grep` search for the start and stop codon(s).

!!! info "Annotating non-coding sequences"

    Depending on the work that you are doing, you may be most interested in identifying protein coding regions or untranslated genetic elements such as ribosomal and transfer RNA sequences (rRNA and tRNAs, respectively).

    For features such as these, which are functional but not not translated, a different set of tools is required for prediction. These will be mentioned later in this section but we will not be working with them today.

Similar to assembly we can perform gene prediction in either a reference-guided manner or through the use of *ab initio* prediction tools. We will not be covering the reference-guided approach, as it is quite simple to perform in `Geneious`, but it is not to be underestimated as a technique - particularly when working with viruses or other organisms with complex read-through or splicing properties.

*Ab initio* prediction is akin to *de novo* assembly - the tool is created with some internal models for what coding regions look like, which are then applied to query sequences to find putative coding regions.

Depending on the intended use of the tool, each prediction tool may be better tuned for partiular assumptions of the data. We are going to use two different tools today, one designed for prediction of prokaryotic coding sequences (which generally lack introns) and one designed primarily for eukaryotic sequences, where splicing is common.

---

## Limitations of protein coding predictions

As you will see from both examples above, protein coding prediction is at best a good starting point for identifying genes. Careful validation of each sequence needs to be performed if you are trying to produce a comprehensive annotation.

In addition, these tools are **_only_** for prediction of protein coding sequences. There are many other genomic features you may need to find in order to find your markers of interest. Some additional tools to examine if you are looking for a complete annotation:

1. [Metaxa2](https://microbiology.se/software/metaxa2/) ([Bengtsson-Palme *et al*, 2015](http://dx.doi.org/10.1111/1755-0998.12399)) - Prediction of small and large subunit ribosomal RNA sequences
1. [Barrnap](https://github.com/tseemann/barrnap) - Prediction of small and large subunit ribosomal RNA sequences
1. [ARAGORN](https://www.ansikte.se/ARAGORN/) ([Laslett & Canback, 2004](https://doi.org/10.1093%2Fnar%2Fgkh152)) - Prediction of tRNA and tm RNA sequences
1. [Infernal](http://eddylab.org/infernal/) ([Nawrocki & Eddy *et al*, 2013](https://doi.org/10.1093%2Fbioinformatics%2Fbtt509)) - Prediction of non-coding RNA sequences, requires [rfam database](https://rfam.org/)

---
