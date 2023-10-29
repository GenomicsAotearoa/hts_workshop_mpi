# 4.5 - Classification of sequences with *k*-mer composition profiles

!!! clock "time"

    * Teaching: 15 minutes
    * Teaching: 15 minutes

!!! circle-info "Objectives and Key points"

    #### Objectives

    * Understand when we may need to use something other than `BLAST` for getting quick classification results.
    * Know how to run `kraken2` using a standard database, and interpret the output files.

    #### Keypoints
    
    * *k*-mer based classification tools provide a rapid means of assigning taxonomy to sequence data.
    * The quality of the classification is tool and database dependent, results are also harder to interpret than a BLAST output
      * For example, there are no identity/coverage statistics from which judgement calls can be made.

---

## Why not `BLAST`?

As we saw in the level 1 training, running `BLASTn` to classify sequences is a good way to get annotation information about a sequence, and from this we can infer taxonomic origin of the sequence. However, obtaining classification information is not the primary purpose of `BLASTn` so we need to perform some manual inspection of the output to make a call on the classification of each sequence.

When we do not care about functional information, and *only* want classification of sequences, there are a suite of tools designed specifically with this purpose in mind. Rather than perform pairwise sequence alignment between query and target sequences, these tools process their reference database into a compressed form containing sequence fingerprints.

Query sequences are fingerprinted and the results compared to the database providing the most likely, **_lowest common ancestor_** for each query sequence.

!!! note "Options for *k*-mer classification"

    How this fingerprinting is done is tool specific, and goes outside the realm of biology. There are a number of tools developed for this kind of work, including;

    1. [kraken2](https://github.com/DerrickWood/kraken2/) ([Wood *et al*., 2019](https://doi.org/10.1186/s13059-019-1891-0))
    1. [CLARK](http://clark.cs.ucr.edu/) ([Ounit *et al*., 2015](https://doi.org/10.1186/s12864-015-1419-2))
    1. [MetaCache](https://github.com/muellan/metacache) ([Müller *et al*., 2017](https://doi.org/10.1093/bioinformatics/btx520))
    1. [ganon](https://gitlab.com/rki_bioinformatics/ganon) ([Piro *et al*., 2020](https://doi.org/10.1093/bioinformatics/btaa458))

We will be using `kraken2` today, as it is a tool which is used often within our operations.

---

## Getting started with `kraken2`

Navigate to `/nesi/project/nesi03181/phel/USERNAME/level2/annotation_kraken2/` and prepare to run `kraken2`.

!!! question "Exercise"

    Find and load the most current version of `kraken2`. Once loaded, run with the help paramter (`-h`) to confirm that the module has successfully loaded.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            module spider kraken2

            module purge
            module load Kraken2/2.1.3-GCC-11.3.0

            kraken2 -h
            ```

            ??? success "Output"

                ```
                Usage: kraken2 [options] <filename(s)>

                Options:
                --db NAME               Name for Kraken 2 DB
                                        (default: "/opt/nesi/db/Kraken2/standard-2018-09")
                --threads NUM           Number of threads (default: 1)
                --quick                 Quick operation (use first hit or hits)
                --unclassified-out FILENAME
                                        Print unclassified sequences to filename
                --classified-out FILENAME
                                        Print classified sequences to filename
                --output FILENAME       Print output to filename (default: stdout); "-" will
                                        suppress normal output
                --confidence FLOAT      Confidence score threshold (default: 0.0); must be
                                        in [0, 1].
                --minimum-base-quality NUM
                                        Minimum base quality used in classification (def: 0,
                                        only effective with FASTQ input).
                --report FILENAME       Print a report with aggregrate counts/clade to file
                --use-mpa-style         With --report, format report output like Kraken 1's
                                        kraken-mpa-report
                --report-zero-counts    With --report, report counts for ALL taxa, even if
                                        counts are zero
                --report-minimizer-data With --report, report minimizer and distinct minimizer
                                        count information in addition to normal Kraken report
                --memory-mapping        Avoids loading database into RAM
                --paired                The filenames provided have paired-end reads
                --use-names             Print scientific names instead of just taxids
                --gzip-compressed       Input files are compressed with gzip
                --bzip2-compressed      Input files are compressed with bzip2
                --minimum-hit-groups NUM
                                        Minimum number of hit groups (overlapping k-mers
                                        sharing the same minimizer) needed to make a call
                                        (default: 2)
                --help                  Print this message
                --version               Print version information

                If none of the *-compressed flags are specified, and the filename provided
                is a regular file, automatic format detection is attempted.
                ```

One of the nice perks of working with `kraken2` is that the maintainers of the tool provide pre-computed, up-to-date [databases for classification](https://benlangmead.github.io/aws-indexes/k2). These are released in many different variants - some for use of HPCs and some for local machines, and databases which cover a range of target organisms.

!!! info "We use `PlusPFP` database, which contains the following lineages:"

    1. Archaea
    1. Bacteria
    1. Fungi
    1. Human
    1. Plant
    1. Protozoa
    1. Virus
    1. Plasmid sequences
    1. Sequencing and transformation vectors (UniVec)

    ---

    For training, we are using a reduced version of this database (`PlusPFP-16`) so that our jobs do not need to queue as long, but the same lineages are represented in this smaller database.

This makes it a comprehensive database for most purposes, but if you are tring to annotate insect sequences this is not the right database for you. For training purposes, a copy of this database can be found at `/nesi/project/nesi03181/phel/databases/` and we also have a copy in our diagnostic workspace.

!!! warning "Make sure you have the right database for your work"

    It is critical to know what is in your classification database, as a failure to include the correct reference material in the database will result in sequences being either misclassified, or not classified at all.

---

## Performing basic classification with `kraken2`

To begin, we are going to write a `slurm` script to run `kraken2` with a fairly minimal set in input parameters just to see the tool in action. Although `kraken2` is very fast to run, it requires a lot of memory (RAM) to run so we cannot run it from the command line in most cases.

Complete the following script, then submit your `slurm` job.

!!! file-code "kraken2_basic.sl"

    ```bash
    #!/bin/bash -e
    #SBATCH --account       nesi03181
    #SBATCH --job-name      kraken2_basic
    #SBATCH --time          00:01:00
    #SBATCH --cpus-per-task 16
    #SBATCH --mem           20G
    #SBATCH --error         kraken2_basic.%j.err
    #SBATCH --output        kraken2_basic.%j.out

    module purge
    module load Kraken2/2.1.3-GCC-11.3.0

    cd /nesi/project/nesi03181/phel/USERNAME/level2/annotation_kraken2/

    PFP_DB="/nesi/project/nesi03181/phel/databases/k2_pluspfp_16gb_20231009"
    kraken2 --db ${PFP_DB} --threads ${SLURM_CPUS_PER_TASK} --use-names \
        --output outputs/input_seqs.out \
        input/input_seqs.fna
    ```

!!! note "No email block in the `slurm` script"

    This job will take less than 30 seconds to complete, so we're not going to bother with the email notifications for job status. The only reason we need to push this job through `slurm` at all is because of the memory required to run the tool.

!!! terminal "code"

    ```bash
    sbatch kraken2_basic.sl
    ```

    ??? success "Output"

        ```
        Submitted batch job XXXXXXXX
        ```

This will not take long to run.

---

## The `kraken2` output format

Once the job is complete, use `less` or `head` to take a look at the output file. The output file is a tab-delimited table with the following columns:

!!! file-code "kraken2_basic.sl"

    |Column|Content|
    |:---:|:---|
    |1|A single letter (C or U) designating whether the sequence was classified or not|
    |2|Sequence name|
    |3|The [NCBI taxonomy ID](https://www.ncbi.nlm.nih.gov/taxonomy) applied to the sequence (0 if unclassified)<br />Because we used the `--use-names` flag, this takes the form `TAXONOMY NAME (TAXONOMY ID)`<br />By default this is just the numeric value|
    |4|Sequence length<br />For paired-end data, this will be in the form `FORWARD|REVERSE`|
    |5|A list of how many *k*-mers mapped to each taxonomy ID<br />Results are reported as colon-delimited pairings of `TAXID:NUMBER`|

It is also possible to generate a report file that summarises the number of sequences classified to each taxonomy/rank, using the modified command:

!!! terminal "code"

    ```bash
    kraken2 --db ${PFP_DB} --threads ${SLURM_CPUS_PER_TASK} --use-names \
        --output outputs/input_seqs.out \
        --report outputs/input_seqs.report.txt \
        input/input_seqs.fna
    ```

Depending on your situation, this may be more useful.

!!! info "Classification output versus report"

    **When would you need the per-sequence output**

    If you are screening a sequence set for the specific reads obtained from a particular organism, or are trying to remove contamination from a sample, having the per-sequence results are invaluable.

    The content of this file can be easily parsed with the `cut` and `grep` tools to extract lists of sequences match or failing to match different search criteria.

    ---
    
    **When would you need the high-level report**

    If you trust that your data is uncontaminated, you might be trying ot achive one of the following:
    
    1. Detecting the presence of a particular organism in a mixed population (i.e. eDNA surveillance).
    1. Confirming the identity of a pure culture.
       1. In this case you might be comfortable taking a majority-rules approach to your data - if most of the sequences are classified as a single species, this is likely the identity of your culture.

---

## Refining the `kraken2` classification

By default `kraken2` does not apply any filtering to the classification data. Reads are classified based on the most commonly reported taxonomy for a sequence, with no regard for how good the evidence for this classification can be. There is a built-in parameter to increase the stringency of classification. The filter is applied during classification and cannot be run after the fact.

This means that if you want to change your filtering criteria you must perform classification again. This differs from tools like `BLAST`, where a more stringent filter can be applied after running classification.

>For example, if we run a `BLAST` search looking to report all hits with greater than 30% identity, we can later filter the table to keep results with greater than 50% identity by filtering the appropriate column.

!!! warning "Limitations with the `kraken2` filtering"

    It is also difficult to interpret the numeric values of filters.

    According to the `kraken2` manual, the confidence scoring does not have a simple interpretation.

    For example, with a `BLAST` identitiy of 90%, we can infer that 90% of positions are identical between query and target match, or a bootstrap support of 90% would mean that 90% of the resamplings agree with the full result.

    However, there is no such interpretation of what a `0.9` confidence value in `kraken2` means, other than it is more strict that a value of `0.5`.

We generally don't need to apply much of a cutoff to `kraken2` - the tool is already good at identifying sequences which can't be reliably classified.

!!! file-code "kraken2_refined.sl"

    ```bash
    #!/bin/bash -e
    #SBATCH --account       nesi03181
    #SBATCH --job-name      kraken2_refined
    #SBATCH --time          00:01:00
    #SBATCH --cpus-per-task 16
    #SBATCH --mem           20G
    #SBATCH --error         kraken2_refined.%j.err
    #SBATCH --output        kraken2_refined.%j.out

    module purge
    module load Kraken2/2.1.3-GCC-11.3.0

    cd /nesi/project/nesi03181/phel/USERNAME/level2/annotation_kraken2/

    PFP_DB="/nesi/project/nesi03181/phel/databases/k2_pluspfp_16gb_20231009"
    kraken2 --db ${PFP_DB} --threads ${SLURM_CPUS_PER_TASK} --confidence 0.1 --use-names \
        --output outputs/input_seqs.conf_0.1.out \
        input/input_seqs.fna
    ```

!!! terminal "code"

    ```bash
    sbatch kraken2_refined.sl
    ```

    ??? success "Output"

        ```
        Submitted batch job XXXXXXXX
        ```

---

## Comparing the outputs

Without going into too much detail, we can briefly compare the contents of both output files. For starters, we can use some of the basic `bash` command line tools to work out how many sequences were in our input and output files.

!!! terminal "code"

    ```bash
    grep -c ">" input/input_seqs.fna 
    ```

    ??? success "Output"

        ```
        5614
        ```

!!! terminal "code"

    ```bash
    wc -l outputs/input_seqs.* 
    ```

    ??? success "Output"

        ```
         5614 outputs/input_seqs.conf_0.1.out
         5614 outputs/input_seqs.out
        11228 total
        ```

The number of results in each output file is the same, and also the same as the number of sequences in the input file. This is different to `BLAST`, where an unclassified or unannotated sequence is not reported in the output.

A quick way to look for differences in the classification rate of the unfiltered and filtered `kraken2` output is to apply the following commands:

!!! terminal "code"

    ```bash
    cut -f ${n} outputs/input_seqs.out | sort | uniq -c
    ```

If we point this towards the correct column, we can count how many sequences were successfully classified or not.

!!! question "Exercise"

    Look through these tutorial notes, or the `kraken2` documentation, to find the appropriate column index for the command above. Once you have determined the correct value for `${n}`, complete the command above and run it for both `kraken2` outputs.

    How many sequences were classified for each run?

    ??? circle-check "Solution"

        The first column (C or U values) is the easiest one to work with, although you could also get the same result from column 3.

        ```bash
        n=1
        cut -f ${n} outputs/input_seqs.out | sort | uniq -c
        cut -f ${n} outputs/input_seqs.conf_0.1.out | sort | uniq -c
        ```

        ??? success "Output"

            ```
            # outputs/input_seqs.out
            5205 C
             409 U

            # outputs/input_seqs.conf_0.1.out
              10 C
            5604 U
            ```

As you can see, even a small degree of filtering has a drastic impact on the classification rate. If you are applying a filter value, using a low threshold is recommended to avoid discarding too much usable data.

---

## How reliable are these results?

The input sequences in this exercise are simulated Oxford Nanopore reads from the following organisms;

1. Tomato (*Solanum lycopersicum*) (5,447 sequences)
1. *Xyllela fastidiosa* (166 sequences)
1. Potato virus X (1 sequence)

If you skim through the results of both output files, you should note the following:

??? file-code "outputs/input_seqs.out"

    |Lineage|Observations|
    |:---|:---:|
    |*Solanum lycopersicum*|4,976|
    |unclassified|409|
    |*Xylella fastidiosa*|162|
    |*Homo sapiens*|37|
    |*Solanum* subgen. *Lycopersicon*|6|
    |*Solanum pennellii*|2|
    |*Triticum aestivum*|2|
    |*Asparagus officinalis*|1|
    |*Flavobacterium* sp.|1|
    |*Hordeum vulgare*|1|
    |*Lolium rigidum*|1|
    |*Oryza brachyantha*|1|
    |*Pandoraea thiooxydans*|1|
    |Potato virus X|1|
    |*Setaria viridis*|1|
    |*Streptomyces cinnabarinus*|1|

??? file-code "outputs/input_seqs.conf_0.1.out"

    |Lineage|Observations|
    |:---|:---:|
    |unclassified|5,604|
    |*Homo sapiens*|6|
    |*Solanum lycopersicum*|2|
    |Potato virus X|1|

What should be apparent here is that while filtering did remove a lot of the noise from our sample (we lost species we know are not truly present) we also lost a signal that **_was_** in the sample.

!!! warning "Recommendation for `kraken2`"

    Don't let this put you off using `kraken2` - it is an excellent tool for rapidly sifting through HTS data for quick summaries. However, results should only be considered indicative in nature, not definitive.

    Use `kraken2` to make a quick check for names of interest (regulated organisms, host, contamination) but follow up with more thorough investigation.

---
