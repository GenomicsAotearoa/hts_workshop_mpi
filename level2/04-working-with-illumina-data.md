
# Working with Illumina data

* Teaching: ??? minutes
* Exercises: ??? minutes

#### Objectives

* Know how to assess the quality of Illumina sequence data using visualisation tools such as `FastQC`
* Be able to perform quality filtering to remove adapter and barcode regions, as well as low quality sequence spans

#### Keypoints

* ...

---

## Contents

1. [Assessing sequence quality](#assessing-sequence-quality)

---

## Assessing sequence quality

When we obtain data from a sequencing facility, it is always important to check the overall quality of the sample, and to confirm whether or not sequencing constructs have been removed from the seqence data. Im particular, we want to know:

1. Are there adapter and/or barcode sequences attached to the reads
1. Are there any obvious low-quality regions of sequence
1. Is there a quality drop-off towards the end of read-pair sequence which might necessitate trimming

A very useful tool for answering these questions is `FastQC`. This tool takes a set of fastq files as input and produces reports for each one to allow us to answer the questions above, as well as examine over features of the sequences such as compositional bias, *k*-mer frequency profiles, and sequence duplication levels.

To activate `FastQC` on NeSI, we need to load the `slurm` module using the command

```bash
$ module purge
$ module load FastQC/0.11.7
```

>**NOTE:** Remember that it is always useful to begin a new session with the `module purge` command, particularly when working through the JupyterHub portal. If you have already run this command without loading new software, you can skip this command.

We can run `FastQC` from the command line, either one file at a time, or by passing a list of files using the wildcard (`*`) operator:

```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/shell_data/untrimmed_fastq/
$ mkdir -p fastq_output/

$ fastqc -o fastq_output/ SRR097977.fastq SRR098026.fastq
```

`FastQC` generates output reports in `.html` files and a `.zip` file that contains the main display resources. These can be viewed in a standard web browser. Since we are connection to NeSI using the JupyterHub system, we can view these directly:

1. Click on the folder icon in the top left to open the folder navigator pane (if not already open).
1. Use the file browsing system to navigate through to `/nesi/project/nesi03181/phel/USERNAME/shell_data/untrimmed_fastq/fastq_output/`
1. Double click on the output `.fastqc.html` files to open them in the a new tab

Let's now look at some of the outputs, starting with the summary for sample `SRR097977`.

![](../img/02_fastqc_overview.png)

This is the basic view for `FastQC` output. At the left-hand side of the tab is a navigation menu, which can move you quickly through the pages of summary information. Alternatively, you can simply scroll down the page to find the section you are most interested in.

From the summary view, the main point of interest is the 'Basic Statistics' table. This gives you some brief summary information for your input file, such as the name of the file read (this can be important if you are working with many files), the fastq encoding, the numbers of sequences and the average length. You should already have a rough expectation for these numbers based from correspondence with your sequencing provider.

Scrolling down (or clicking on the 'Per base sequence quality' link) will take us to the main piece of information we wish to know about the samples - the overall quality of the sequences.

![](../02_fastqc_quality.png)

This view provides us with a nice graphical summary of the average sequence quality along the length of our reads. Fastq Q-scores are ranked on the y-axis and the nucleotide position in the read (or range of positions, for reads which are several hundred nucleotides in length) are plotted sequentially along the x-axis. The sample `SRR097977` shows a severe, but expected pattern of quality, whereby the sequence quality degrades as the reading window moves towards the right-hand side of the sequence.

This view provides us with two pieces of information - how strictly we need to trim our sequences, and what effect we could expect to see on our number of sequences and sequence length after quality filtering.

![](../02_fastqc_content.png)


![](../02_fastqc_overrepresented.png)


![](../02_fastqc_adapters.png)



>**NOTE:** `FastQC` does not load the forward and reverse pairs of a library in the same window, you need to be mindful of how your samples relate to each other and in which order you have opened them.

# TODO - ADD FIGURE





