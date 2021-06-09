
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
$ mkdir -p miseq_data/

$ fastqc -o miseq_data/ sample1_R1.fastq.gz
$ fastqc -o miseq_data/ sample1_R2.fastq.gz
$ fastqc -o miseq_data/ sample2_R1.fastq.gz
$ fastqc -o miseq_data/ sample2_R2.fastq.gz

$ fastqc -o miseq_data/ sample*.fastq.gz
```

`FastQC` generates output reports in `.html` files and a `.zip` file that contains the main display resources. These can be viewed in a standard web browser. Since we are connection to NeSI using the JupyterHub system, we can view these directly:

1. Click on the folder icon in the top left to open the folder navigator pane (if not already open).
1. Use the file browsing system to navigate through to `/nesi/project/nesi03181/phel/USERNAME/miseq_data/`
1. Double click on the output `.fastqc.html` files to open them in the a new tab

>**NOTE:** `FastQC` does not load the forward and reverse pairs of a library in the same window, you need to be mindful of how your samples relate to each other and in which order you have opened them.

# TODO - ADD FIGURE





