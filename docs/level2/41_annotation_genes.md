# 4.1 - Prediction of protein coding sequences

!!! clock "time"

    * Teaching: 30 minutes
    * Exercises: 30 minutes

!!! circle-info "Objectives and Key points"

    #### Objectives
    
    * Understand the limitations of protein coding predictions when working with annotations
    * Understand how to run `prodigal` for predicting protein coding regions in prokaryotic genomes
    * Understand how to run `AUGUSTUS` for predicting protein coding regions in genomes

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

*Ab initio* prediction is akin to *de novo* assembly - the tool is created with some internal models for what coding regions look like, which are then applied to query sequences to find putative coding regions. Depending on the intended use of the tool, each prediction tool may be better tuned for partiular assumptions of the data. We are going to use two different tools today, one designed for prediction of prokaryotic coding sequences (which generally lack introns) and one designed primarily for eukaryotic sequences, where splicing is common.

---

## Predicting protein coding regions with `prodigal`

We are going to be predicting protein coding sequences in a reference *M. bovis* genome using the tool `prodigal` ([Hyatt *et al.*, 2010](https://doi.org/10.1186/1471-2105-11-119)). This is a powerful prediction tool which is quick to run, and flexible enough for most projects.

Navigate to the `/nesi/project/nesi03181/phel/USERNAME/level2/annotation_prodigal/` directory and prepare to run `prodigal`.

!!! question "Exercise"

    Use your knowledge of `slurm` to find and load the more current version of `prodigal` available on NeSI. Once loaded, run with the help paramter (`-h`) to confirm that the module has successfully loaded.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            # One of the following:
            module spider prodigal

            module avail prodigal

            # Load and confirm
            module purge
            module load prodigal/2.6.3-GCCcore-7.4.0

            prodigal -h
            ```

            ??? success "Output"

                ```
                Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
                                [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
                                [-p mode] [-q] [-s start_file] [-t training_file] [-v]

                        -a:  Write protein translations to the selected file.
                        -c:  Closed ends.  Do not allow genes to run off edges.
                        -d:  Write nucleotide sequences of genes to the selected file.
                        -f:  Select output format (gbk, gff, or sco).  Default is gbk.
                        -g:  Specify a translation table to use (default 11).
                        -h:  Print help menu and exit.
                        -i:  Specify FASTA/Genbank input file (default reads from stdin).
                        -m:  Treat runs of N as masked sequence; don't build genes across them.
                        -n:  Bypass Shine-Dalgarno trainer and force a full motif scan.
                        -o:  Specify output file (default writes to stdout).
                        -p:  Select procedure (single or meta).  Default is single.
                        -q:  Run quietly (suppress normal stderr output).
                        -s:  Write all potential genes (with scores) to the selected file.
                        -t:  Write a training file (if none exists); otherwise, read and use
                            the specified training file.
                        -v:  Print version number and exit.
                ```

There are quite a few options given here, which we can split into input and output parameters.

!!! info "Input parameters"

    |Parameter|Function|
    |:---:|:---|
    |`-i`|Specify FASTA/Genbank input file (default reads from stdin).|
    |`-p`|Select procedure (single or meta).  Default is single.|
    |`-g`|Specify a translation table to use (default 11).|

    The first parameter here should be obvious, so will not be discussed.

    The `-p` parameter is an important one to pay attention to, as it determines whether or not we are predicting sequences in an single genome (i.e. a genome obtained from a pure culture) or a metagenomic mix of sequences. In `single` mode prediction goes through a round of training against the input sequences before producing a prediction output tailored to your contigs. The `meta` mode uses pre-calculated profiles of gene features. Strictly speaking it is possible to run `meta` mode on an isolate genome, or `single` mode on a metagenome but the results do differ.

    Selecting the translation table (`-g`) is something we usually do not need to wory about. `prodigal` supports numbers 1 through 25 of the [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). By default, it will start with code 11 (bacteria, archaea, and plastids) but shift to code 4 if the predictions are too short. This is important to us because genetic code 4 corresponds to the *Mycoplasmataceae* - the bacterial family that contains the genera *Mycoplasma* and *Ureaplasma*, and the family *Spiroplasmataceae*. These lineages have repurposed `UGA` from a stop codon to a tryptophan codon.

!!! info "Output parameters"

    |Parameter|Function|
    |:---:|:---|
    |`-a`|Write protein translations to the selected file.|
    |`-d`|Write nucleotide sequences of genes to the selected file.|
    |`-o`|Specify output file (default writes to `stdout`).|
    |`-f`|Select output format. Default is gbk.|

    Three of these capture the prediction information in commonly used formats. The output of `-d` and `-a` are simply fasta format, and the `-o` (or `stdout`) output are in GenBank format. This is really the same information reported in multiple different ways, but the different files have different uses so it's recommended to take a few of them.

    At the very least, we should be capturing either `-d` **_and_** `-a`, or `-o` to get the full prediction. It is helpful to have both the nucleotide and protein sequence of a coding region, as the different annotation techniques can be applied to each one, and certain methods of phylogenetic analysis favour one sequence type over the other.

One of the nice features of `prodigal` is that it does not take a lot of resources to run, so we can easily run it without resorting to `slurm` for a single genome.

!!! terminal "code"

    ```bash
    prodigal -p single -g 4 \
        -i input/M_bovis.NZ_CP005933.fna \
        -d outputs/M_bovis.NZ_CP005933.prod.fna \
        -a outputs/M_bovis.NZ_CP005933.prod.faa \
        -o outputs/M_bovis.NZ_CP005933.prod.gbk
    ```

    ??? success "Output"

    ```
    -------------------------------------
    PRODIGAL v2.6.3 [February, 2016]         
    Univ of Tenn / Oak Ridge National Lab
    Doug Hyatt, Loren Hauser, et al.     
    -------------------------------------
    Request:  Single Genome, Phase:  Training
    Reading in the sequence(s) to train...948516 bp seq created, 29.29 pct GC
    Locating all potential starts and stops...30752 nodes
    Looking for GC bias in different frames...frame bias scores: 2.38 0.43 0.19
    Building initial set of genes to train from...done!
    Creating coding model and scoring nodes...done!
    Examining upstream regions and training starts...done!
    -------------------------------------
    Request:  Single Genome, Phase:  Gene Finding
    Finding genes in sequence #1 (948516 bp)...done!
    ```

This will only take a few seconds to complete. Once done, we can quickly see how many protein coding sequences were predicted by counting the number of sequences in either of the output fasta files.

!!! terminal "code"

    ```bash
    grep -c ">" outputs/M_bovis.NZ_CP005933.prod.fna outputs/M_bovis.NZ_CP005933.prod.faa
    ```

    ??? success "Output"

    ```
    outputs/M_bovis.NZ_CP005933.prod.fna:792
    outputs/M_bovis.NZ_CP005933.prod.faa:792
    ```

Since this is a chromosome downloaded from the [NCBI RefSeq database](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP005933) we can look to the official annotation and see how many proteins we should expect to find. A copy of the protein coding sequences from this genome have been provided in your `outputs/` folder.

!!! question "Exercise"

    Use `grep` to count the number of protein coding sequences in the official annotation for this *M. bovis* genome.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            grep -c ">" outputs/M_bovis.canonical.faa
            ```

            ??? success "Output"

                ```
                765
                ```

                This number is slightly fewer than what `prodigal` predicts for us, but the numbers are very close. We would need to examine the predictions in detail to understand what is causing the difference.

If we inspect the output of the fasta files, you will notice that there is a lot of metadata stored in each line. For example:

!!! terminal "code"

    ```bash
    head -n1 outputs/M_bovis.NZ_CP005933.prod.faa
    ```

    ??? success "Output"

    ```
    NZ_CP005933.1_1 # 1 # 1401 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.293
    ```

There is quite a lot to unpack here - some is quite important to know and other parts are just descriptive. We can break down the results like so:

!!! file-code "prodigal metadata"

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

!!! info "Interpreting the prodigal header"

    As `prodigal` is only concerned with *predicting* sequences, not annotating them, there is no functional imformation which can be used to guide the name of each prediction. For simplicity, protein predictions are simply named as `[CONTIG NAME]_[PREDICTION]`. This has some nice implications when we want to study gene synteny but for the most part numbering off the predicitons in the order they are made is the simpliest way to generate names.

    The next three numbers provide the nucleotide coordinates of the coding sequence and the orientation (1 for forward, -1 for reverse). Again, these are very useful when tracking down the genomic context of a sequence and can sometimes provide a quick-and-dirty means for spotting rearrangements.

    After these parameters are done, there are a series of keyword-linked pieces of information about the sequence. The unique gene identifier is used to link the entries in the fasta file to the output obtained from the `-o` (or `stdout`) channel. Similar to the prediction number, it is simply derived from the order of contigs and sequence of predictions.

    The main piece of information we are going to inspect is whether or not the prediction is complete of not. The 'partial' keyword provides a two digit code that reports the status of the prediction, the first digit corresponds to the start of the sequence and the second to the end of the sequence. The values these can take are either 0 (complete) or 1 (incomplete). A fully complete prediction will therefore have a code of `00`, and a partial predictions can be:

    * `01` - Started, but no end found - typically when a prediction runs off the end of a contig
    * `10` - No start identified, but a complete end found - often a sequence which occurs at the start of the contig
    * `11` - No ends found - likely due to predicting for a very short contig

    You can also inspect the sequences to see if they appear complete. Typically protein predictions begin with a methionine (M) amino acid residue, as this is the translation of the `ATG` codon. There is no residue which corresponds to stop (as stop is a gap in translation) so `prodigal` reports the stop position with an asterisk (`*`).

---

## Predicting protein coding regions with `AUGUSTUS`

Unlike prokaryotic genomes, the genes of eukaryotes carry intronic sequences which need to be spliced out of the gene sequence before undergoing translation. The detection of splicing boundaries is a difficult task, as there are many organism-specific patterns used to mark splice sites. Protein prediction tools for this purpose typically come with  number of pre-trained models for finding protein domains within contigs, but if there is no model for your organism, or a closely related lineage, then results may not be ideal.

A recently published article ([Scalzitti *et al.*, 2020](https://doi.org/10.1186/s12864-020-6707-9)) which profiled a number of these tools found `AUGUSTUS` to be one of the best performing tools for gene prediction in eukaryotic organisms, so this is what we will use today.

>**Note:** `AUGUSUTUS` does require training against a closely related model organisms to generate accurate predictions, which we do not have for this workshop. We will instead be performing predictions with a few different models and seeing how the outputs differ.

!!! question "Exercise"

    Find the latest version of `AUGUSTUS` on NeSI and load it.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            module spider augustus

            module purge
            module load AUGUSTUS/3.5.0-gimkl-2022a

            augustus
            ```

            ??? success "Output"

                ```
                AUGUSTUS (3.5.0) is a gene prediction tool.
                Sources and documentation at https://github.com/Gaius-Augustus/Augustus

                usage:
                augustus [parameters] --species=SPECIES queryfilename

                'queryfilename' is the filename (including relative path) to the file containing the query sequence(s)
                in fasta format.

                SPECIES is an identifier for the species. Use --species=help to see a list.

                parameters:
                ...

                For a complete list of parameters, type "augustus --paramlist". A description of the important ones can be found in the file RUNNING-AUGUSTUS.md.
                ```

The sequence(s) we have to work with today are from the brown marmorated stink bug (*Halyomorpha halys*) and if we use the prompt from the exercise above to view the available species models you will see that there is no good model for this organism, or even a closely related species from the *Pentatomidae*.

We will use two different models for an initial round of prediction on the *Halyomorpha halys* sequences - one insect and one bacterial species.

!!! terminal "code"

    ```bash
    augustus --species=help
    ```

    ??? success "Output"

        ```
        usage:
        augustus [parameters] --species=SPECIES queryfilename

        where SPECIES is one of the following identifiers

        identifier                               | species
        -----------------------------------------|----------------------
        pea_aphid                                | Acyrthosiphon pisum
        aedes                                    | Aedes aegypti
        amphimedon                               | Amphimedon queenslandica
        ancylostoma_ceylanicum                   | Ancylostoma ceylanicum
        adorsata                                 | Apis dorsata
        honeybee1                                | Apis mellifera
        arabidopsis                              | Arabidopsis thaliana
        ...
        (maize5)                                 | Zea mays
        ```

We will run `AUGUSTUS` on a randomly selected insect model, and see how well the results compare to the expected number of proteins in the sequences. Pick whatever model seems most likely to you, and use the following command to run `AUGUSTUS`.

!!! terminal "code"

    ```bash
    augustus --genemodel=partial --protein=on --codingseq=on --species=honeybee1 input/NW_020110202.fna > outputs/NW_020110202.aug_hb1.gff
    ```

This will take about 15 minutes to run, so while it is running set up the follow exercise:

!!! question "Exercise"

    Select a second model organism, preferably a non-insect model, and begin a second round of prediction using `AUGUSTUS`.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            augustus --genemodel=partial --protein=on --codingseq=on --species=E_coli_K12 input/NW_020110202.fna > outputs/NW_020110202.aug_ecoli.gff
            ```

Once your jobs have finished (about 15 minutes per attempt), we need to run a helper scfript that bundles with `AUGUSTUS` to extract the gene and coding sequence predicitons from the output file.

!!! terminal "code"

    ```bash
    getAnnoFasta.pl outputs/NW_020110202.aug_hb1.gff
    ```

This will return no output to the console, but creates two new files, with names dervied from the file name above.

1. `*.codingseq` - Nucleotide sequences for each predicted coding region
1. `*.aa` - Amino acid residue translation from each predicted coding region

!!! question "Exercise"

    Use `grep` to count the number of protein coding sequences in your predicted files (both models) and the official annotation for this *H. halys* genome.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            grep -c ">" outputs/NW_020110202.aug_hb1.aa outputs/NW_020110202.aug_ecoli.aa outputs/NW_020110202.canonical.faa 
            ```

            ??? success "Output"

                ```
                outputs/NW_020110202.aug_hb1.aa:170
                outputs/NW_020110202.aug_ecoli.aa:17
                outputs/NW_020110202.canonical.faa:110
                ```

Similar to when running predictions with `prodigal`, you can see that neither of the models identified the correct number of protein coding sequences. This is because neither of these gene models are accurate for the organism we are trying to characterise. However, the more closely related model was much closer in its prediction that the non-insect model.

??? info "Creating a custom species profile for gene prediciton"

    Creating a new model is a slow process so we will not cover it here, although you can use the steps below if this is something you want to try.

    **Creating a custom gene profile**

    There are a few steps we need to perform in advance of the new model training. The first is to do with file permissions - the location of the prediction databases that `AUGUSTUS` uses for gene prediction are not writable to us, so we cannot add new data into them. We must create our own copy of the configuration information and point `AUGUSTUS` towards this new location in order to create new models

    !!! terminal "code"

        ```bash
        cp -r ${AUGUSTUS_CONFIG_PATH} /your/location/to/store/data/

        # Update the `AUGUSTUS_CONFIG_PATH` variable to point to the new location
        AUGUSTUS_CONFIG_PATH="/your/location/to/store/data"

    Now we just need to obtain a reference genome to train against. For *H. halys*, this can found on the [NCBI website](https://www.ncbi.nlm.nih.gov/assembly/GCA_000696795.3) and downloaded from the command line:

    !!! terminal "code"

        ```bash
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/696/795/GCA_000696795.3_Hhal_1.1/GCA_000696795.3_Hhal_1.1_genomic.gbff.gz
        gunzip GCA_000696795.3_Hhal_1.1_genomic.gbff.gz
        ```

    Once these steps are completed, training is a single command:

    !!! terminal "code"

        ```bash
        autoAugTrain.pl --trainingset=GCA_000696795.3_Hhal_1.1_genomic.gbff --species=hhalys
        ```

    We could then perform prediction using the new model as usual. Note that if this was performed in a new session, we would need to set the `AUGUSTUS_CONFIG_PATH` variable again.

    !!! terminal "code"

        ```bash
        AUGUSTUS_CONFIG_PATH=""

        augustus --AUGUSTUS_CONFIG_PATH=/your/location/to/store/data --genemodel=partial --protein=on --codingseq=on --species=hhalys ...
        ```

---

## Limitations of protein coding predictions

As you will see from both examples above, protein coding prediction is at best a good starting point for identifying genes. Careful validation of each sequence needs to be performed if you are trying to produce a comprehensive annotation.

In addition, these tools are **_only_** for prediction of protein coding sequences. There are many other genomic features you may need to find in order to find your markers of interest. Some additional tools to examine if you are looking for a complete annotation:

1. [Metaxa2](https://microbiology.se/software/metaxa2/) ([Bengtsson-Palme *et al*, 2015](http://dx.doi.org/10.1111/1755-0998.12399)) - Prediction of small and large subunit ribosomal RNA sequences
1. [Barrnap](https://github.com/tseemann/barrnap) - Prediction of small and large subunit ribosomal RNA sequences
1. [ARAGORN](https://www.ansikte.se/ARAGORN/) ([Laslett & Canback, 2004](https://doi.org/10.1093%2Fnar%2Fgkh152)) - Prediction of tRNA and tm RNA sequences
1. [Infernal](http://eddylab.org/infernal/) ([Nawrocki & Eddy *et al*, 2013](https://doi.org/10.1093%2Fbioinformatics%2Fbtt509)) - Prediction of non-coding RNA sequences, requires [rfam database](https://rfam.org/)

---
