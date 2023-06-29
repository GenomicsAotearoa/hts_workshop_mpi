# Key points before filtering and trimming Illumina data

* Teaching: 45 minutes
* Exercises: 15 minutes

#### Objectives

* Be familiar with the consequences of paired-end sequencing data, compared with single-end sequencing.

#### Keypoints

* When working with paired-end sequence data, maintaining the order of sequences between forward and reverse sequence files is critical.
* Be aware of how to account for singleton (orphan) sequences when working with paired-end sequencing data.

---

## Trimming paired data and preserving read order

Once we have visualised our sequence quality, we need to make some decisions regarding whether we perform quality filering or not, how severely the data must be trimmed, and whether or not we need to remove adapter sequences. There are many tools available for trimming sequence data, and today we are only going to work with one - [fastp](https://github.com/OpenGene/fastp). This tool was selected for its speed and performance, as well as having a good suite of quality-of-life features, but it there are many alternate tools out there which can perform the same job.

Before we actually start to work on our data, there are some important factors of the nature of paired-end sequencing data that we need to consider. At it's most basic level, a tool like `fastp` simply takes one input sequence file, trims the sequences according to some specifications, then saves the output to a new file:

```bash
$ fastp -i in.R1.fq.gz -o out.R1.fq.gz
```

Sequences from the input file are read and trimmed for regions that do or do not match filtering criteria. For example, reject adapter regions will be removed, as will regions of the sequence where the average quality falls below a required threshold. If this means that the entirity of a sequence falls below the selection criteria, the resulting sequence will not be written to the output file. Often, short sequences will also be rejected from the output.

If our sequences consist of paired reads, then filtering is a slightly more complicated process than what we saw above. You cannot just run the trimming tool over each sequence file:

```bash
$ fastp -i in.R1.fq.gz -o out.R1.fq.gz
$ fastp -i in.R2.fq.gz -o out.R2.fq.gz
```

This is because that between the paired sequence files, the forward/reverse reads are ordered such that the first sequence in the `R2` (reverse) file is the partner to the first sequence in the `R1` (forward file). If a sequence in one of these files is rejected by the filtering system but it's partner passes filtering, then there will be an uneven number of sequences written into the two output files.

This is a massive problem, because when we perform tasks like assembly and read mapping, the tools will assume that the sequence order is identical between forward and reverse files and that there is a spatial relationship between these paired sequences. For example, if we have a pair of files where the sequences align like so:

```
R1-----------------------R2
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

And after filtering we might get an output like:

```diff
  R1-----------------------R2
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

Where `Seq4_R1` and `Seq7_R2` are removed. Any downstream analysis tools we use will not be aware of these gaps, and read the pairing information as:

```diff
  R1-----------------------R2
+ @Seq1-----------------@Seq1
+ @Seq2-----------------@Seq2
+ @Seq3-----------------@Seq3
- @Seq5-----------------@Seq4
- @Seq6-----------------@Seq5
- @Seq7-----------------@Seq6
+ @Seq8-----------------@Seq8
+ @Seq9-----------------@Seq9
```

Where there are now some mismatched pairs - `Seq5_R1` pairs with `Seq4_R2`, `Seq6_R1` with `Seq5_R2`, and `Seq7_R1` with `Seq6_R2`. The impact this has varies depending on what is being done with the next stage of analysis, but whether these errors cause chimeric assembly or simply a failure to successfully map a sequence pair to a reference genome, it is entirely an aretfact of our analysis procedure and not the biology of the sample.

>**NOTE:** This is just a simple example where a single sequence was removed from each direction. In practice the problem would be much more severe, as there would almost certainly be an uneven number of reads filtered from each direction so the order would almost certainly never recover, as it did in this example.

Any worthwhile trimming tool is aware of this problem, and has a built-in capacity for filtering these issues. In the case of `fastp`, this is the correct way to filter the sequences:

```bash
$ fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
```

With this usage, `fastp` takes a second input with the `-I` parameter which is understood to be the pair of the `-i` input. Similarly, the `-O` output file is the output for the `-I` file, and pairs with the `-o` output file. If the example sequences above were run in this manner, we would get the following output:

```diff
  out.R1.fq.gz---out.R2.fq.gz
+ @Seq1-----------------@Seq1
+ @Seq2-----------------@Seq2
+ @Seq3-----------------@Seq3
+ @Seq5-----------------@Seq5
+ @Seq6-----------------@Seq6
+ @Seq8-----------------@Seq8
+ @Seq9-----------------@Seq9
```

Which is great. We have lost a few sequences, but those that we retained are still ordered.

However, although the `Seq4_R1` and `Seq7_R2` sequence were of poor quality and did not pass quality filtering, their partner sequence did pass. In the usage above, the `Seq4_R2` and `Seq7_R1` sequence are still filtered out of the data even though they are of sufficient quality. This is a fairly common occurance, so filtering tools have an additional output type for these unpaired (also called '*orphan*' or '*singleton*') sequences which pass filtering when their partner do not.

The command is getting a bit long, so we will use the `\` character to break the command over several lines:

```bash
$ fastp -i in.R1.fq.gz -I in.R2.fq.gz \
        -o out.R1.fq.gz -O out.R2.fq.gz \
        --unpaired1 out.R1.unpaired.fq.gz \
        --unpaired2 out.R2.unpaired.fq.gz
```

This time the output would look like:

```
out.R1.fq.gz---out.R2.fq.gz
@Seq1-----------------@Seq1
@Seq2-----------------@Seq2
@Seq3-----------------@Seq3
@Seq5-----------------@Seq5
@Seq6-----------------@Seq6
@Seq8-----------------@Seq8
@Seq9-----------------@Seq9

out.R1.unpaired.fq.gz
@Seq7_R1

out.R2.unpaired.fq.gz
@Seq4_R2
```

This is a much better solution, as we are now retaining more high quality data. In practice, singleton sequences are rare but they are still valuable data and worth retaining where possible.

`fastp` has one more nice feature that we can exploit - if we don't set a value for the `--unpaired2` parameter, the sequences which *would* have gone to this file are instead sent to the same place as `--unpaired1`. This means that we can capture all of our unpaired reads in a single file.

```bash
$ fastp -i in.R1.fq.gz -I in.R2.fq.gz \
        -o out.R1.fq.gz -O out.R2.fq.gz \
        --unpaired1 out.unpaired.fq.gz \
```

This time the output would look like:

```
out.R1.fq.gz---out.R2.fq.gz
@Seq1-----------------@Seq1
@Seq2-----------------@Seq2
@Seq3-----------------@Seq3
@Seq5-----------------@Seq5
@Seq6-----------------@Seq6
@Seq8-----------------@Seq8
@Seq9-----------------@Seq9

out.unpaired.fq.gz
@Seq4_R2
@Seq7_R1
```

---
