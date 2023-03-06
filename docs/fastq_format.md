# The FASTQ file format

In several fo the exercises in this training course we are dealing with fastq files, although we never actually discuss what they contain. This is a file format inteded to convery nucleotide sequence information, similar to a fasta file, as well as per-nucleotide quality information which is expressed as the confidence that a position in the sequence was assigned correctly by the sequencing platform. A comprehensive description of the fastq format and its variatiosn can be found in [this article](https://en.wikipedia.org/wiki/FASTQ_format) but for the purposes of these workshops there are only a few key points you need to know.

Although it looks complicated, the format is actually pretty easy to understand once you know the pattern. Each fastq file is comprised of repeating blocks of 4 lines, where each line describes a single sequence read. Within these four lines:

|Line|Description|
|:---:|:---|
|1|Is the name of the sequence, and optionally some metadata about the read.<br />This line always begins with a `@` symbol.|
|2|The nucleic acid sequence.|
|3|A placeholder line which divides the nucleic acid sequence and quality information.<br />In some old files this line will contain a `+` symbol followed by the sequence name.<br />Newer tools tend to just use the `+` alone to keep the file smaller.|
|4|A string of characters which represent the quality scores.<br />This string must be the same length as Line 2.|

Using one of the example files from the training exercises (for example, [Level 1, module 3](../level1/13_shell_manipulation.md#viewing-the-contents-of-files)), if we run the `head` command with a `-n` parameter of 4 we can view the first complete read in one of the files.

```bash
$ head -n4 SRR097977.fastq
@SRR097977.1 209DTAAXX_Lenski2_1_7:8:3:710:178 length=36
TATTCTGCCATAATGAAATTCGCCACTTGTTAGTGT
+SRR097977.1 209DTAAXX_Lenski2_1_7:8:3:710:178 length=36
CCCCCCCCCCCCCCC>CCCCC7CCCCCCACA?5A5<
```

Line 4 shows the quality for each nucleotide in the read. Quality is interpreted as the probability of an incorrect base call (e.g. 1 in 10) or, equivalently, the base call accuracy (e.g. 90%). These probability values are the results from the base calling algorithm and dependent on how much signal was captured for the base incorporation.  To make it possible to line up each individual nucleotide with its quality score, the numerical score is converted into a code where each individual character represents the numerical quality score for an individual nucleotide. In the line above, the quality score line is `CCCCCCCCCCCCCCC>CCCCC7CCCCCCACA?5A5<`.

The numerical value assigned to each of these characters depends on the sequencing platform that generated the reads. The sequencing machine used to generate our data uses `PHRED+33` encoding, which is what is used in all modern Illumina and Nanopore sequencing, as well as Sanger sequencing. Each character is assigned a quality score between 0 and 42 as shown in the chart below.

```
Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJK
                  |         |         |         |         |
Quality score:    0........10........20........30........40..
```

If you line it up, you can see that the `C` character has a qualtiy score of 34, which is the most common value in this sequence. The probability of a position with a Q of 34 can be calculated as:

$$ P_{incorrect} = 10^{\frac{Q}{-10}} = 0.000398107 $$

Which can be converted to the probabilty that the read is correct by

$$ P_{correct} = 1 - P_{incorrect} = 0.9996 $$

Therefore, for any of these characters the base call accuracy is 99.96%. You will almost never need to calculate the per-position score, but it is handy to know the following numbers for when we are quality filtering:

|Q|Accuracy|
|:---:|:---:|
|10|90.00%|
|20|99.00%|
|30|99.90%|
|40|99.99%|

---
