# 4.3 - Interpreting BLAST results

!!! clock "time"

    * Teaching: 30 minutes
    * Exercises: 30 minutes

!!! circle-info "Objectives and Key points"

    #### Objectives

    * Learn how to interpret `BLAST` results 
    * Understand the meaning of *identity*, *coverage*, *e-value*, and *bitscore* metrics
    
    #### Keypoints
    
    * Interpreting the results of BLAST alignments requires thought and attention 

---

## Interpretting the results of BLAST queries

It is important to remember, like most bioinformatics tools, `BLAST` has a specific job. In this case sequence alignment to a set of references. `BLAST` is really good at this job, but it **_does not offer interpretation of its alignments_**.

Interpretation is completely up to the user on a case by case basis. It is therefore important to know your data and to understand the output metrics given by `BLAST` to help you make a biologically usefull interpretation of the results. 

There are typically four main metrics that we need ot check when reviewing a `BLAST` result:

!!! info "Key `BLAST` metrics"

    **Coverage**

    * This value tells us as a percentage how much of our query aligns with the database match. 
    * A small coverage value means only a small part of our sequence has matched. A perfect match would have a coverage of 100. 

    **Identity**

    * This represents the percent of bases which are identical between the query and the database hit, *over the aligned region*.

    **E-value**

    * This is the number of hits equivalent to this hit that we would expect to see by chance alone.
    * Smaller E-values represent better hits, but an exact E-value cut off needs to be decided on a case by case basis.
    * E-values take into account the coverage and identity scores for each hit, and also the size of the database queried.

    **Bit score** 

    * Similarly to E-values, bit scores summarise the sequence similarity between the query and database hit.
    * Bit scores are calculated independently from the query sequence length and the database size, as databases are constantly evolving this makes bit scores a constant statistical indicator.
    * A higher bit score indicates a better hit. 

---

## Examining our file

With this in mind lets look at the results from our `BLAST` job.

Return to your `blast_annotation/` folder, if you left it, and examine the new output file. Take a look at your results using the `less` or `head` command.

!!! terminal "code"

    ```bash
    head output.txt
    ```

    ??? success "Output"

        ```bash
        seq1    gi|1607238104|dbj|AP019558.1|   91.750  1794    148     0       1       1794    366818  3650250.0      2494    Mycoplasma bovis KG4397 DNA, complete genome    28903
        seq1    gi|1441442372|gb|CP022588.1|    91.695  1794    149     0       1       1794    690905  6926980.0      2488    Mycoplasmopsis bovis strain MJ4 chromosome, complete genome     28903
        seq1    gi|1315670167|emb|LT578453.1|   91.695  1794    149     0       1       1794    655989  6577820.0      2488    Mycoplasma bovis isolate JF4278 genome assembly, chromosome: I  28903
        seq1    gi|2507795645|gb|CP058524.2|    91.695  1794    149     0       1       1794    401349  4031420.0      2488    Mycoplasmopsis bovis strain Mb49 chromosome     28903
        seq1    gi|2507793515|gb|CP058496.2|    91.695  1794    149     0       1       1794    320995  3192020.0      2488    Mycoplasmopsis bovis strain Mb222 chromosome    28903
        seq1    gi|2507792460|gb|CP058473.2|    91.695  1794    149     0       1       1794    1055619 10574120.0     2488    Mycoplasmopsis bovis strain VK22 chromosome     28903
        seq1    gi|2507791648|gb|CP058453.2|    91.695  1794    149     0       1       1794    9966    11759 0.0      2488    Mycoplasmopsis bovis strain Mb287 chromosome    28903
        seq1    gi|2507791648|gb|CP058453.2|    91.695  1794    149     0       1       1794    1120377 11221700.0     2488    Mycoplasmopsis bovis strain Mb287 chromosome    28903
        seq1    gi|2506302979|gb|CP058448.2|    91.695  1794    149     0       1       1794    763992  7657850.0      2488    Mycoplasmopsis bovis strain Mb1 chromosome      28903
        seq1    gi|2506302187|gb|CP058432.2|    91.695  1794    149     0       1       1794    354705  3564980.0      2488    Mycoplasmopsis bovis strain Mb216 chromosome    28903
        ```

It looks like a table, but awkwardly there are no column names. However, the names of the columns correspond to the values we provided in our `slurm` script at the start of this session.

|Column header|Meaning|
|:---|:---|
|qseqid|Sequence ID of the query sequence (input file)|
|sseqid|Sequence ID of the target sequence (reference database)|
|pident|Percentage of identical positions between query and target|
|length|Alignment length (sequence overlap) of the common region between query and target|
|mismatch|Number of mismatches between query and target|
|gapopen|Number of gap openings in the alignment|
|qstart|Position in the query sequence where alignment begins|
|qend|Position in the query sequence where alignment ends|
|sstart|Position in the target sequence where alignment begins|
|send|Position in the target sequence where alignment ends|
|evalue|The E-value for the query/target match, as described above|
|bitscore|The bit score for the query/target match, as described above|
|salltitles|Display All Subject Title(s) for the target sequence|
|staxids|Display the NCBI taxid value(s) for the target sequence|

Applying these to the layout, we get something more sensible:

|qseqid|sseqid|pident|length|mismatch|gapopen|qstart|qend|sstart|send|evalue|bitscore|salltitles|staxids|
|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|:---|
|seq1|gi\|1607238104\|dbj\|AP019558.1\||91.750|1794|148|0|1|1794|366818|3650250.0|2494|Mycoplasma bovis KG4397 DNA, complete genome|28903|
|seq1|gi\|1441442372\|gb\|CP022588.1\||91.695|1794|149|0|1|1794|690905|6926980.0|2488|Mycoplasmopsis bovis strain MJ4 chromosome, complete genome|28903|
|seq1|gi\|1315670167\|emb\|LT578453.1\||91.695|1794|149|0|1|1794|655989|6577820.0|2488|Mycoplasma bovis isolate JF4278 genome assembly, chromosome: I|28903|
|seq1|gi\|2507795645\|gb\|CP058524.2\||91.695|1794|149|0|1|1794|401349|4031420.0|2488|Mycoplasmopsis bovis strain Mb49 chromosome|28903|
|seq1|gi\|2507793515\|gb\|CP058496.2\||91.695|1794|149|0|1|1794|320995|3192020.0|2488|Mycoplasmopsis bovis strain Mb222 chromosome|28903|
|seq1|gi\|2507792460\|gb\|CP058473.2\||91.695|1794|149|0|1|1794|1055619 10574120.0|2488|Mycoplasmopsis bovis strain VK22 chromosome|28903|
|seq1|gi\|2507791648\|gb\|CP058453.2\||91.695|1794|149|0|1|1794|9966|11759 0.0|2488|Mycoplasmopsis bovis strain Mb287 chromosome|28903|
|seq1|gi\|2507791648\|gb\|CP058453.2\||91.695|1794|149|0|1|1794|1120377 11221700.0|2488|Mycoplasmopsis bovis strain Mb287 chromosome|28903|
|seq1|gi\|2506302979\|gb\|CP058448.2\||91.695|1794|149|0|1|1794|763992|7657850.0|2488|Mycoplasmopsis bovis strain Mb1 chromosome|28903|
|seq1|gi\|2506302187\|gb\|CP058432.2\||91.695|1794|149|0|1|1794|354705|3564980.0|2488|Mycoplasmopsis bovis strain Mb216 chromosome|28903|

!!! question "Exercise"

    Download the `output.txt` file do your computer and open it in `Excel`. Look through the sequence annotations, and determine the most likely organism that each sequence was obtained from.

    Are there are values in the results table that you are skeptical about? If so, raise them as a discussion with the group.

    ??? circle-check "Solution"

        * seq1 = *Mycoplasma bovis*
        * seq2 = *Bactrocera ritsemai*
        * seq3 = Pepino mosaic virus
        * seq4 = No hit
        * seq5 = *Fusarium oxysporum* (reversed)

---
