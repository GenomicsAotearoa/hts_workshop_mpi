# Filtering and trimming Illumina data

!!! clock "time"

    * Teaching: 15 minutes
    * Exercises: 20 minutes

!!! circle-info "Learning objectives"

    **Objectives**
    
    * Be able to perform quality filtering to remove adapter and barcode regions, as well as low quality sequence regions.
    
    **Key points**
    
    * Sequences can be filtered to remove low quality sequences and sequencing artefacts using `fastp`.

---

## Basic execution of `fastp`

The last session was quite a lot of talking, not a lot of typing, so lets start with a quick exercise using `fastp` just to see it in action. Navigate to `/nesi/project/nesi03181/phel/<username>/level1/quality_illumina/` and try run `fastp` in a few different ways:

!!! note "fastp reporting"

    `fastp` produces a lot of output text while running, and creates some report files. You can ignore these for now.

!!! terminal "code"

    ```bash
    module load fastp/0.23.2-GCC-11.3.0
    ```

!!! question "Exercise"

    Load the `fastp` tool and have a quick test run of the different behaviours. Run the tool a few different times to achieve the following outputs with the `miseq_R*.fq.gz` files.

    1. Filter with default settings, retaining paired outputs only.
    1. Filter with default settings, writing unpaired sequences to individual files.
    1. Filter with default settings, writing unpaired sequences to the same file.

    ??? circle-check "Solution"

        1)
        !!! terminal "code"

            ```bash
            fastp -i reads/miseq_R1.fq.gz -I reads/miseq_R2.fq.gz -o results/ex1_1.fq.gz -O results/ex1_2.fq.gz
            ```

        2)
        !!! terminal "code"

            ```bash
            fastp -i reads/miseq_R1.fq.gz -I reads/miseq_R2.fq.gz \
                  -o results/ex2_1.fq.gz -O results/ex2_2.fq.gz \
                  --unpaired1 results/ex2_1_single.fq.gz \
                  --unpaired2 results/ex2_2_single.fq.gz
            ```

        3)
        !!! terminal "code"

            ```bash
            fastp -i reads/miseq_R1.fq.gz -I reads/miseq_R2.fq.gz \
                  -o results/ex3_1.fq.gz -O results/ex3_2.fq.gz \
                  --unpaired1 results/ex3_single.fq.gz
            ```

Now that we understand the importance of keeping our fastq outputs ordered, and how to achieve this with `fastp`, let's look at the actual filtering process. From the command line, you can view the manual for `fastp` by running it either without any parameters, or with the `--help` parameter.

!!! terminal "code"

    ```bash
    fastp
    # OR
    fastp --help
    ```

There are quite a few options we can modify here, and we are going to ignore most of them in the interest of time. The two main aspects of trimming that we want to work through are the adapter removal options, and the quality filtering parameters.

---

## Adapter removal

Traditionally, sequencing trimming tools require you to provide the adapters in your data set via a fasta file. The sequences in this file are then matched against the start and end of each fastq sequence read by the tool and trimming occurs when there is a match. `fastp` does provide this option, via the `--adapter_fasta` parameter, i.e.:

!!! terminal "code"

    ```bash
    fastp --adapter_fasta my_adapter_list.fasta -i ... -o ...
    ```

But it introduces two quality of life improvements that are often easier to work with. Since there should only be a single adapter in your data, rather than create a file of adapter sequences, the sequences can be passed directly to the program, saving the effort of creating the file and also making your work more transparent when viewing `history` logs:

!!! terminal "code"

    ```bash
    fastp --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
          --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
          -i ... -o ...
    ```

If we do not know what adapters were used, `fastp` has the ability to auto-detect adapters in the sequences. When running in this way, `fastp` will read the first 1,000,000 sequences in the file and search for common nucleotide patterns at the start of the read. The sequences that `fastp` identifies as the adapter are reported as part of the console output, and can then be compared against known lists of sequencing adapters to confirm they are correct.

!!! terminal "code"

    ```bash
    fastp --detect_adapter_for_pe -i ... -o ...
    ```

This is a very helpful feature, but be warned that if your library has already been trimmed of adapters and you force `fastp` to find them, results can be quite weird.

!!! question "Exercise"

    Run `fastp` against the `ERR4179828` files with and without adapter trimming. Run the command again with automatic adapter detection. Look at the outputs in `FastQC` and see how they differ. Don't worry about capturing the singleton sequences for this exercise.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            fastp --disable_adapter_trimming \
                  -i reads/ERR4179828_1.fq.gz -I reads/ERR4179828_2.fq.gz \
                  -o results/ERR4179828_1.no_trim.fq.gz -O results/ERR4179828_2.no_trim.fq.gz

            fastp --detect_adapter_for_pe \
                  -i reads/ERR4179828_1.fq.gz -I reads/ERR4179828_2.fq.gz \
                  -o results/ERR4179828_1.auto_trim.fq.gz -O results/ERR4179828_2.auto_trim.fq.gz

            fastqc -O results/ results/ERR4179828_1.*_trim.fq.gz
            ```
---

## Quality filtering

Although it's the last feature we address in this exercise, quality filtering is really the main reason we use tools like `fastp`. There are several ways to screen a sequence for quality, once we have determined what we want our lower limit of acceptable quality to be.

1. Remove a sequence with an average quality below this value.
1. Cut the regions of a sequence with quality below this value.

The later of these is definitely the better option, as a sequence with quite high average quality can still have large spans of low-quality sequence which will interfere with out analysis.

This option is usually refered to as 'sliding window' evaluation and consists of the sequence being read from one end to the other, reading a small number of bases (the window) and assessing the average quality over these nucleotides. The window progresses (slides) along the sequence, and at the first instance of average quality dropping below the threshold, the sequence is cut and the remainder discarded.

In this kind of analysis, the default threshold is a Q-score of 20, and the window size is 4. These are fairly standard values, but can be changed if desired. This method is enabled by the `--cut_right` parameter in `fastp`, as the sliding window moves from the left-hand side (start) of the sequence to the right-hand side (end), and everything discarded after the cut point is therefore the right half of the sequence.

After any form of trimming, there is also the option to remove sequences which fall below a required minimum sequence length. Length filtering is appplied *after* trimming, so that the length of a sequence is evaluated once low quality regions are removed.

`fastp` does have other ways of trimming low quality regions, in which it begins from the start or end of the sequence and removes nucleotides below the quality threshold until a high quality base is encountered. When working with Illumina data this method is probably not very valuable, as sequence quality usually starts off quite high with Illumina data and tends to degrade ove time. However, there are some sequencing platforms that can yield poor quality at the start of a sequence, so it is worth knowing that this behaviour is possible.

!!! question "Exercise"

    Inspect the help options for `fastp` and formulate a command to perform sliding window filtering with a window length of 8 nucleotides, and a quality threshold of 15.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            fastp --cut_right --cut_window_size 8 --cut_mean_quality 15 -i ... -I ... -o ... -O ...
            ```

---

## Other features

There are just a few more features of `fastp` which can be useful to know.

1. `--thread`: Increase the number of computer processes (threads) used to process the sequences. Using more threads typically makes the job faster to complete. The default value is 2.
1. `--html`/`--json`: Change the names of the report files generated during a run. This is very useful when processing many samples, as the default output names (`fastp.html` and `fastp.json`) are recycled on each use and successive runs of `fastp` overwrite the output of previous runs.
1. `--length_required`: Discard sequences shorter than the specified length (after trimming) regardless of their quality.
1. `--n_base_limit`: Discard sequences with more than the specified limit of `N` characters, regardless of the sequence quality.

---

````
