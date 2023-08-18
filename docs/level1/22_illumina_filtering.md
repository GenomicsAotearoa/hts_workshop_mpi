# 2.2 - Filtering Illumina data

!!! clock "time"

    * Teaching: 10 minutes
    * Exercises: 20 minutes

!!! circle-info "Learning objectives"

    **Objectives**
    
    * Be familiar with the consequences of paired-end sequencing data, compared with single-end sequencing.
    * Be able to perform quality filtering to remove adapter and barcode regions, as well as low quality sequence regions.
    
    **Key points**
    
    * When working with paired-end sequence data, maintaining the order of sequences between forward and reverse sequence files is critical.
    * Be aware of how to account for singleton (orphan) sequences when working with paired-end sequencing data.
    * Sequences can be filtered to remove low quality sequences and sequencing artefacts using `fastp`.

---

## Considerations when working with paired-end data

Once we have visualised our sequence quality, we need to make some decisions regarding whether we perform quality filering or not, how severely the data must be trimmed, and whether or not we need to remove adapter sequences.

There are many tools available for trimming sequence data, and today we are only going to work with one - [fastp](https://github.com/OpenGene/fastp). This tool was selected for its speed and performance, as well as having a good suite of quality-of-life features, but it there are many alternate tools out there which can perform the same job.

Before we actually start to work on our data, there are some important factors of the nature of paired-end sequencing data that we need to consider.

At it's most basic level, a tool like `fastp` simply takes one input sequence file, trims the sequences according to some specifications, then saves the output to a new file. Sequences from the input file are read and trimmed for regions that do or do not match filtering criteria.

For example, adapter regions may be removed along with regions where the average sequence quality falls below a required threshold. This means that the entirity of a sequence falls below the selection criteria, the resulting sequence will not be written to the output file. Often, short sequences will also be rejected from the output.

When working with Illumina data, each fragment sequenced is read from the forward and reverse direction so we have two corresponding **_paired sequences_** representing each piece of DNA in our library.

When working with this paired data, it is critical to preserve the ordering between the forwrad and reverse sequences so that their relationship to each other is maintained.

??? note "The consequences of treating forwad and reverse sequence files as independent data"

    If our sequences consist of paired reads you cannot just run the trimming tool over each sequence file independently because the forward/reverse reads are ordered such that the first sequence in the reverse file is the partner to the first sequence in the forward file. If a sequence in one of these files is rejected by the filtering system but it's partner passes filtering, then there will be an uneven number of sequences written into the two output files.

    This is a massive problem, because when we perform tasks like assembly and read mapping, the tools will assume that the sequence order is identical between forward and reverse files and that there is a spatial relationship between these paired sequences. For example, if we have a pair of files where the sequences align like so:

    !!! file-code "Observed input ordering"

        ```
        @Seq1-----------------@Seq1
        @Seq2-----------------@Seq2
        @Seq3-----------------@Seq3
        @Seq4-----------------@Seq4
        @Seq5-----------------@Seq5
        @Seq6-----------------@Seq6
        @Seq7-----------------@Seq7
        @Seq8-----------------@Seq8
        @Seq9-----------------@Seq9
        ```

    And after filtering we might observe an output like the following, where `Seq4_R1` and `Seq7_R2` are removed.

    !!! file-code "Expected output ordering"

        ```diff
        + @Seq1-----------------@Seq1
        + @Seq2-----------------@Seq2
        + @Seq3-----------------@Seq3
        - ----------------------@Seq4
        + @Seq5-----------------@Seq5
        + @Seq6-----------------@Seq6
        - @Seq7----------------------
        + @Seq8-----------------@Seq8
        + @Seq9-----------------@Seq9
        ```

    Any downstream analysis tools we use will not be aware of these gaps, and read the pairing information as:

    !!! file-code "Observed output ordering"

        ```diff
        + @Seq1-----------------@Seq1
        + @Seq2-----------------@Seq2
        + @Seq3-----------------@Seq3
        - @Seq5-----------------@Seq4
        - @Seq6-----------------@Seq5
        - @Seq7-----------------@Seq6
        + @Seq8-----------------@Seq8
        + @Seq9-----------------@Seq9
        ```

    The consequence of these gaps should be obvious - `Seq5_R1` pairs with `Seq4_R2`, `Seq6_R1` with `Seq5_R2`, and `Seq7_R1` with `Seq6_R2`. The impact this has varies depending on what is being done with the next stage of analysis, but whether these errors cause chimeric assembly or simply a failure to successfully map a sequence pair to a reference genome, it is entirely an aretfact of our analysis procedure and not the biology of the sample.

Any worthwhile trimming tool is aware of this problem, and has a built-in capacity for filtering these issues.

---

## Performing paired-end sequence trimming with `fastp`

We are going to use the tool `fastp` to practice filtering some Illumina sequence data, making sure to preserve read order between the forward and reverse sequences. Begin by navitating to the `quality_illumina/` folder and loading the module for `fastp`.

!!! terminal "code"

    ```bash
    module load fastp/0.23.2-GCC-11.3.0

    cd /nesi/project/nesi03181/phel/<username>/level1/quality_illumina/
    ```

We are going to build a command for filtering the paired files `reads/miseq_R1.fq.gz` and `reads/miseq_R2.fq.gz`. The most basic implementation of `fastp` is to take the input sequences, apply and trimming, and then save the paired sequences where **_both forward and reverse sequence_** passed the filtering criteria.

!!! terminal "code"

    ```bash
    fastp -i reads/miseq_R1.fq.gz -I reads/miseq_R2.fq.gz -o results/miseq_R1.fq.gz -O results/miseq_R2.fq.gz
    ```

!!! note "Interpreting the command parameters"

    We have just provided four parameters to `fastp`. There are two conventions you need to be aware of when interpreting these:

    * The `-i` and `-o` flags are used to specify _I_nput and _O_utput files. The filtered sequences read from the file provided with `-i` will be written to the destination specified with `-o`.

    * The lowercase and uppercase convention is used to signify to `fastp` that these are paired-end files. The file specified with `-I` is understood to be the reverse partner of the file specified with `-i`, and its output is written to `-O`.

!!! question "Exercise"

    Why can we not just run `fastp` in the following way?
    
    !!! terminal "code"

        ```bash
        fastp -i reads/miseq_R1.fq.gz -o results/miseq_R1.fq.gz
        fastp -i reads/miseq_R2.fq.gz -o results/miseq_R2.fq.gz
        ```

    ??? circle-check "Solution"

        The paired-end relationship between sequences in each file will not be maintained. If a sequence passes quality filtering in the forward file, but its partner in the reverse file fails the result ill still be written to the output file (and *vice versa*).

As mentioned above, this will only write out sequences where both the forward and reverse partner made it throug the trimming criteria. In any situation where a read passed filtering but its partner failed, that read is lost even though it contains 'accceptable' data.

This is a fairly common occurance, so filtering tools have an additional output type for these unpaired (also called '*orphan*' or '*singleton*') sequences which pass filtering when their partner do not. We can modify our `fastp` command to retain data when these cases happen, directing the sequences to new files:

!!! terminal "code"

    ```bash
    fastp -i reads/miseq_R1.fq.gz -I reads/miseq_R2.fq.gz -o results/miseq_R1.fq.gz -O results/miseq_R2.fq.gz --unpaired1 results/miseq_R1s.fq.gz --unpaired2 results/miseq_R2s.fq.gz
    ```

??? note "Schematic of how the output files will be populated"

    Calling back to the description above of the hypothetical data set of 10 sequences, we would now expect our output files to look something like:

    !!! file-code "miseq_R1.fq.gz / miseq_R2.fq.gz"

        ```diff
        + @Seq1---------------------@Seq1
        + @Seq2---------------------@Seq2
        + @Seq3---------------------@Seq3
        + @Seq5---------------------@Seq5
        + @Seq6---------------------@Seq6
        + @Seq8---------------------@Seq8
        + @Seq9---------------------@Seq9
        ```

    !!! file-code " miseq_R1s.fq.gz"

        ```diff
        + @Seq7_R1
        ```

    !!! file-code "miseq_R2s.fq.gz"

        ```diff
        + @Seq4_R2
        ```

    We have retained those individual reads, without compromising the pairing order.

`fastp` has one nice feature that we can exploit - if we don't set a value for the `--unpaired2` parameter, the sequences which *would* have gone to this file are instead sent to the same place as `--unpaired1`. This means that we can capture all of our unpaired reads in a single file with a more concise command:

!!! terminal "code"

    ```bash
    fastp -i reads/miseq_R1.fq.gz -I reads/miseq_R2.fq.gz -o results/miseq_R1.fq.gz -O results/miseq_R2.fq.gz --unpaired1 results/miseq_s.fq.gz
    ```

??? note "Schematic of how the output files will be populated"

    !!! file-code "miseq_R1.fq.gz / miseq_R2.fq.gz"

        ```diff
        + @Seq1---------------------@Seq1
        + @Seq2---------------------@Seq2
        + @Seq3---------------------@Seq3
        + @Seq5---------------------@Seq5
        + @Seq6---------------------@Seq6
        + @Seq8---------------------@Seq8
        + @Seq9---------------------@Seq9
        ```

    !!! file-code "miseq_s.fq.gz"

        ```diff
        + @Seq7_R1
        + @Seq4_R2
        ```

---

## Customising the quailty filtering parameters

We talked a lot above about how our reads were being filtered but not what this actually means. In general, there are two ways to screen a sequence for quality, once we have determined what we want our lower limit of acceptable quality to be.

1. Remove a sequence with an average quality below this value.
1. Cut the regions of a sequence with quality below this value.

The later of these is definitely the better option, as a sequence with quite high average quality can still have large spans of low-quality sequence which will interfere with out analysis.

A very popular method of achieving this is using a 'sliding window' approach. This consists of the sequence being read from one end to the other, reading a small number of bases (the window) and assessing the average quality over these nucleotides. The window progresses (slides) along the sequence, and at the first instance of average quality dropping below the threshold, the sequence is cut and the remainder discarded.

In this kind of analysis, the default threshold is a Q-score of 20, and the window size is 4. These are fairly standard values, but can be changed if desired.

??? note "Changing the default trimming values"

    This method is enabled by the `--cut_right` parameter in `fastp`, as the sliding window moves from the left-hand side (start) of the sequence to the right-hand side (end), and everything discarded after the cut point is therefore the right half of the sequence.

    An example of how you could customise a sliding window approach in `fastp` is:

    !!! terminal "code"

        ```bash
        fastp --cut_right --cut_window_size 8 --cut_mean_quality 15 -i ... -I ... -o ... -O ...
        ```

    Which would assess an 8 nucleotide window size and trim at the first instance of the average quality falling below Q15.

After any form of trimming, there is also the option to remove sequences which fall below a required minimum sequence length. Length filtering is appplied *after* trimming, so that the length of a sequence is evaluated once low quality regions are removed.

---

## Other useful features of `fastp`

There are just a few more features of `fastp` which can be useful to know.

* `--thread`
   * Increase the number of computer processes (threads) used to process the sequences. Using more threads typically makes the job faster to complete. The default value is 2.
* `--html`/`--json`
   * Change the names of the report files generated during a run.
   * This is very useful when processing many samples, as the default output names (`fastp.html` and `fastp.json`) are recycled on each use and successive runs of `fastp` overwrite the output of previous runs.
* `--length_required`
  * Discard sequences shorter than the specified length (after trimming) regardless of their quality.
* `--n_base_limit`
  * Discard sequences with more than the specified limit of `N` characters, regardless of the sequence quality.
* `--adapter_sequence`
  * Specify a sequencing adapter sequence to be identified and removed from the start and end of reads.
  * *Normally* this is performed by the sequencing facility, but not always (as we saw in the `FastQC` exercises).
* `--detect_adapter_for_pe`
  * As above, but with this parameter `fastp` will attempt to identify the adapter sequence by screening the input sequences and looking for conserved regions at the beginning of each read.

### Example of manual adapter removal

!!! terminal "code"

    ```bash
    fastp --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  -i ... -o ...
    ```

### Example of automatic adapter removal

!!! terminal "code"

    ```bash
    fastp --detect_adapter_for_pe  -i ... -o ...
    ```

---
