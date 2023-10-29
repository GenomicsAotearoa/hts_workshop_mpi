# 4.2 - Prediction with `prodigal`

!!! clock "time"

    * Teaching: 10 minutes
    * Exercises: 15 minutes

!!! circle-info "Objectives and Key points"

    #### Objectives
    
    * Understand how to run `prodigal` for predicting protein coding regions in prokaryotic genomes

    #### Keypoints

    * Selecting the correct prediction mode (genome or metagenome), and translation table is key to getting good outputs.
    * `prodigal` reports output in several formats - you do not need all of the output files, but know which ones you want to keep and why.

---

## Getting started and loading the tool

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

---

## Predicting protein coding regions

One of the nice features of `prodigal` is that it does not take a lot of resources to run, so we can easily run it without resorting to `slurm` for a single genome.

!!! terminal "code"

    ```bash
    prodigal -p single -g 4 -i input/M_bovis.NZ_CP005933.fna \
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

!!! question "Which value of `-g` are we using?"

    As mentioned above, for most cases the standard genetic code is the correct starting place, but we are working with an *M. bovis* sequence, which uses [translation table 4](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG4).

This will only take a few seconds to complete. Once done, we can quickly see how many protein coding sequences were predicted by counting the number of sequences in either of the output fasta files.

---

## Interpretting the output format

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
    |              |   |   |      |   |      |
    |              |   |   |      |   |      Gene completeness
    |              |   |   |      |   Unique gene ID
    |              |   |   |      Orientation
    |              |   |   Stop position
    |              |   Start position
    |              Gene suffix
    Contig name
    ```

!!! info "Interpreting the prodigal header"

    **Contig name** and **Gene suffix**

    As `prodigal` is only concerned with *predicting* sequences, not annotating them, there is no functional information which can be used to guide the name of each prediction. For simplicity, protein predictions are simply named as `[CONTIG NAME]_[PREDICTION]`. This has some nice implications when we want to study gene synteny but for the most part numbering off the predicitons in the order they are made is the simpliest way to generate names.

    ---

    **Start**, **Stop** and **Orientation**

    The next three numbers provide the nucleotide coordinates of the coding sequence and the orientation (1 for forward, -1 for reverse). Again, these are very useful when tracking down the genomic context of a sequence and can sometimes provide a quick-and-dirty means for spotting rearrangements.

    ---

    **Unique gene ID**

    After these parameters are done, there are a series of keyword-linked pieces of information about the sequence. The unique gene identifier is used to link the entries in the fasta file to the output obtained from the `-o` (or `stdout`) channel. Similar to the prediction number, it is simply derived from the order of contigs and sequence of predictions.

    ---

    **Gene completeness**

    The main piece of information we are going to inspect is whether or not the prediction is complete of not. The 'partial' keyword provides a two digit code that reports the status of the prediction, the first digit corresponds to the start of the sequence and the second to the end of the sequence. The values these can take are either 0 (complete) or 1 (incomplete). A fully complete prediction will therefore have a code of `00`, and a partial predictions can be:

    * `01` - Started, but no end found - typically when a prediction runs off the end of a contig
    * `10` - No start identified, but a complete end found - often a sequence which occurs at the start of the contig
    * `11` - No ends found - likely due to predicting for a very short contig

    You can also inspect the sequences to see if they appear complete. Typically protein predictions begin with a methionine (M) amino acid residue, as this is the translation of the `ATG` codon. There is no residue which corresponds to stop (as stop is a gap in translation) so `prodigal` reports the stop position with an asterisk (`*`).

---
