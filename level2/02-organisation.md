# Project organisation

* Teaching: 15 minutes
* Exercises: 15 minutes

#### Objectives

* Create a file system for a bioinformatics project.
* Understand the value in splitting different types of analysis across different directories.
* Use the `history` command and a text editor like `nano` to document your work on your project.

#### Keypoints

* Spend the time to organise your file system when you start a new project. Your future self will thank you!
* Always save a write-protected copy of your raw data.

---

## Contents

1. [Getting your project started](#getting-your-project-started)
1. [Organising your files](#organising-your-files)
1. [Preparing the data folder](#preparing-the-data-folder)
3. [References](#references)

---

## Getting your project started

Project organisation is one of the most important parts of a sequencing project, and yet is often overlooked amidst the excitement of getting a first look at new data. Of course, while it's best to get yourself organised before you even begin your analyses, it's never too late to start, either.  

You should approach any bioinformatics analysis similarly to how you do a diagnostic sample, and this ideally begins with experimental design. For all of the laboratory steps (collecting specimens, extracting DNA, preparing your samples) you have kept records of what methods, solutions, and reagent batches were used at each step. Documentation doesn't stop just because we're out of the lab.

Genomics projects can quickly accumulate dozens or even hundreds of files spread across multiple folders, as different tools and parameters are tested. Good documentation is key tracking what has been done for future reporting, and which files were generated as part of each stage of analysis.

There are several published guides to the best practices for documenting a genomics project. Due to the nature of our file storage, documentation, and reporting requirements we will not be adhering to these perfectly, but the underlying principle is valuable.

In this exercise we will setup a file system for the work we will be perofrming for the rest of this workshop. First navigate to your `nesi03181` directory, then confirm that you are in the correct directory using the `pwd` command.

```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/
$ pwd
```

> ### Exercise
>
> Use the `mkdir` command to make the following directories:   
> 1. `fastq_processing/`
> 1. `fastq_processing/docs/`
> 1. `fastq_processing/data/`
> 1. `fastq_processing/results/`
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> $ mkdir fastq_processing/
> $ mkdir fastq_processing/docs/
> $ mkdir fastq_processing/data/
> $ mkdir fastq_processing/results/
> ```
> </details>

Use `ls -R` to verify that you have created these directories. The `-R` option for `ls` stands for recursive. This option causes `ls` to return the contents of each subdirectory within the directory iteratively. 

```bash
$ ls -R fastq_processing/
```

You should see the following output:

```
fastq_processing/:
data  docs  results

fastq_processing/data:

fastq_processing/docs:

fastq_processing/results: 
```

---

## Organising your files

Before beginning any analysis, it's important to save a copy of your raw data. The raw data should never be changed. Regardless of how sure you are that you want to carry out a particular data cleaning step, there's always the chance that you'll change your mind later or that there will be an error in carrying out the data cleaning and you'll need to go back a step in the process. Having a raw copy of your data that you never modify guarantees that you will always be able to start over if something goes wrong with your analysis. How your team stores the raw data might vary, but for this exercise we will be creating a copy on NeSI using the `fastq_processing/data/` directory.

You can store any results that are generated from your analysis in the `results/` folder. This guarantees that you won't confuse results file and data files.

The `docs` folder is the place to store any written analysis of your results, notes about how your analyses were carried out, and documents related to your eventual publication. These are extremely important for reproducing your work, or reminding yourself what analysis was performed if you have to revisit the work weeks or months later.

When carrying out wet-lab protocols, we work from a standardised worksheet and keep a hard copy recording variable values such as extraction kit, lot numbers, and reagent expirey dates. This detailed record-keeping process is just as important when doing computational analyses, where we document software names, versions, and the exact parameteres used when running the program. We have already worked with the `history` command to access a list of our recent commands, and this is generally the most reliable way of documenting our work although it is sometimes helpful to curate the information that it returns.

> ### Exercise
>
> Using your knowledge of the shell, use the append redirect `>>` to create a file in the `fastq_processing/docs/` directory with the name `file_layout.XXXX_XX_XX.txt`. Replace the X's with the four-digit year, two-digit month, and two digit day values for the day of the workshop, e.g. `file_layout.2021_06_21.txt`).
> 
> As there will probably be a lot of information in your `history` output which is not relevant to this workshop, use the `tail` command to filter the results and capture only the last 10 or so lines into the output file. 
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> $ history | tail -n10 >> fastq_processing/docs/file_layout.XXXX_XX_XX.txt
> ```
>
> **NOTE:** Here we used the last 10 lines as an example, the number of lines may vary depending on the content of your `history` log.
> </details>

You may have noticed that your history contains the `history` command itself. To remove this redundancy from our log, let's use the `nano` text editor to fix the file:  

```bash
$ nano fastq_processing/docs/file_layout.XXXX_XX_XX.txt
```

(Remember to replace the `XXXX_XX_XX` with your workshop date)

From the `nano` screen, you can use your cursor to navigate, type, and delete any redundant lines.   

> ### Navigating in Nano
>
> Although `nano` is useful, it can be frustrating to edit documents, as you can't use your mouse to navigate to the part of the document you would like to edit. Here are some useful keyboard shortcuts for moving around within a text document in `nano`. You can find more information by typing <kbd>Ctrl</kbd>-<kbd>G</kbd> within `nano`.
> 
> | key     | action |
> | ------- | ---------- |
> | <kbd>Ctrl</kbd>-<kbd>Space</kbd> OR <kbd>Ctrl</kbd>-<kbd>→</kbd> | Move forward one word |
> | <kbd>Alt</kbd>-<kbd>Space</kbd> OR <kbd>Esc</kbd>-<kbd>Space</kbd> OR <kbd>Ctrl</kbd>-<kbd>←</kbd> | Move back one word |
> | <kbd>Ctrl</kbd>-<kbd>A</kbd> | Move to the beginning of the current line |
> | <kbd>Ctrl</kbd>-<kbd>E</kbd> | Move to the end of the current line |
> | <kbd>Ctrl</kbd>-<kbd>K</kbd> | Cut (remove) the current line |
> | <kbd>Ctrl</kbd>-<kbd>U</kbd> | Paste (un-cut) a line removed with cut |
> | <kbd>Ctrl</kbd>-<kbd>W</kbd> | Search |

Add a date line and comment to the line where you have created the directory. Recall that any text on a line after a `#` is ignored by bash when evaluating the text as code. For example:   

```
# 2021/06/21
# Created sample directories for the workshop  
```

Next, remove any lines of the history that are not relevant by navigating to those lines and using your delete key. Save your file and close `nano`.

Your file should look something like this: 

```
# 2021/06/21
# Created sample directories for the workshop

mkdir fastq_processing/
mkdir fastq_processing/docs/
mkdir fastq_processing/data/
mkdir fastq_processing/results/
```

---

## Preparing the data folder

Later in this workshop we will begin to work with some exemplar fastq files generated from Illumina MiSeq and Oxford Nanopore MinION sequencing. We will create a copy the data to treat as our 'raw' data in these exercises.

> ### Exercise
>
> The data we will use later in this workshop is a mixture of the data we used in the last workshop, as well as some new sequence data obtained from the NCBI Sequence Read Archive, and consists of sequence data obtained from a number of *Mycoplasma bovis* isolates. All data can be found in the `/nesi/project/nesi03181/phel/module_2/` directory.
>
> Copy the following files from this directory into your `fastq_processing/data/` folder:
>
> 1. SRR098026.fastq
> 1. SRR097977.fastq
> 1. Mb1_1.fastq.gz
> 1. Mb1_2.fastq.gz
> 1. Mb1.fastq
>
> You can do this in whatever way makes the most sense to you, either one file at a time, with a `for loop`, or by using wild cards. Once you have created a personal copy of the data, use `chmod` to remove the write permissions from **your** verison of the data.
>
> <details>
> <summary>Solution</summary>
>
> **NOTE:** There are many ways to achieve this. The answer below is one method, but it is not the only way to perofrm this task.
>
> ```bash
> $ cd /nesi/project/nesi03181/phel/dwaite/fastq_processing/data/
> $ cp ../../../module_2/*.fastq ./
> $ cp ../../../module_2/*.gz ./
> $ chmod -w *.fastq *.gz
> ```
> </details>

---

## References

The workflow used here doesn't have to be the exact same as what will be used by your team, particularly when we start working with the `slurm` management system in later workshops. You can read the references below for some ideas about what is worth documenting.

1. [Noble WS (2009) A Quick Guide to Organizing Computational Biology Projects. PLoS Comput Biol 5(7): e1000424.](https://doi.org/10.1371/journal.pcbi.1000424)
1. [Sandve GK, Nekrutenko A, Taylor J, Hovig E (2013) Ten Simple Rules for Reproducible Computational Research. PLoS Comput Biol 9(10): e1003285.](https://doi.org/10.1371/journal.pcbi.1003285)

---

[Next lesson](03-intro-to-nesi.md)
