# Project organisation

* Teaching: 15 minutes
* Exercises: 15 minutes

#### Objectives

* Create a file system for a bioinformatics project.
* Explain what types of files should go in your `docs`, `data`, and `results` directories.
* Use the `history` command and a text editor like `nano` to document your work on your project.

#### Keypoints

* Spend the time to organize your file system when you start a new project. Your future self will thank you!
* Always save a write-protected copy of your raw data.

---

## Contents

1. [Getting your project started](#getting-your-project-started)
1. [Organizing your files](#organising-your-files)
1. [Documenting your activity on the project](#documenting-your-activity-on-the-project)
1. [References](#references)

---

## Getting your project started

Project organization is one of the most important parts of a sequencing project, and yet is often overlooked amidst the excitement of getting a first look at new data. Of course, while it's best to get yourself organized before you even begin your analyses, it's never too late to start, either.  

You should approach any bioinformatics analysis similarly to how you do a diagnostic sample, and this ideally begins with experimental design. For all of the laboratory steps (collecting specimens, extracting DNA, preparing your samples) you have kept records of what methods, solutions, and reagent batches were used at each step. Documentation doesn't stop just because we're out of the lab.

Genomics projects can quickly accumulate dozens or even hundreds of files spread across multiple folders, as different tools and parameters are tested. Good documentation is key tracking what has been done for future reporting, and which files were generated as part of each stage of analysis.

With this in mind, let's have a look at the best practices for documenting your genomics project. Your future self will thank you.  

In this exercise we will setup a file system for the project we will be working on during this workshop.  

We will start by creating a directory that we can use for the rest of the workshop. First navigate to your `nesi03181` directory, then confirm that you are in the correct directory using the `pwd` command.

```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/
$ pwd
```

> ### Exercise
>
> Use the `mkdir` command to make the following directories:   
> 1. `dc_workshop/`
> 1. `dc_workshop/docs/`
> 1. `dc_workshop/data/`
> 1. `dc_workshop/results/`
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> $ mkdir dc_workshop/
> $ mkdir dc_workshop/docs/
> $ mkdir dc_workshop/data/
> $ mkdir dc_workshop/results/
> ```
> </details>

Use `ls -R` to verify that you have created these directories. The `-R` option for `ls` stands for recursive. This option causes `ls` to return the contents of each subdirectory within the directory iteratively. 

```bash
$ ls -R dc_workshop/
```

You should see the following output:

```
dc_workshop/:
data  docs  results

dc_workshop/data:

dc_workshop/docs:

dc_workshop/results: 
```

---

## Organising your files

Before beginning any analysis, it's important to save a copy of your raw data. The raw data should never be changed. Regardless of how sure you are that you want to carry out a particular data cleaning step, there's always the chance that you'll change your mind later or that there will be an error in carrying out the data cleaning and you'll need to go back a step in the process. Having a raw copy of your data that you never modify guarantees that you will always be able to start over if something goes wrong with your analysis.

When starting any analysis, you can make a copy of your raw data file and do your manipulations on that file, rather than the raw version. We learned in a previous exercise how to prevent overwriting our raw data files by setting restrictive file permissions. 

You can store any results that are generated from your analysis in the `results/` folder. This guarantees that you won't confuse results file and data files in six months or two years when you are looking back through your files.

The `docs` folder is the place to store any written analysis of your results, notes about how your analyses were carried out, and documents related to your eventual publication.

---

## Documenting your activity on the project

When carrying out wet-lab protocols, we work from a standardised worksheet and keep a hard copy recording variable values such as extraction kit, lot numbers, and reagent expirey dates. This detailed record-keeping process is just as important when doing computational analyses. Fortunately, it's easy to record the bioinformatic steps you've carried out during an analysis.

The `history` command is a convenient way to document all the commands you have used while analysing and manipulating your project files. Let's document the work we have done on our project so far. 

View the commands that you have used so far during this session using `history`:

```bash
$ history
```

The history likely contains many more commands than you have used for the current exercise. Let's view the last several commands that focus on just what we need for this project.   

View the last 10 lines of your `history` by redirecting the output of `history` into the `tail` command.

```bash   
$ history | tail -n 10
```

> ### Exercise
>
> Using your knowledge of the shell, use the append redirect `>>` to create a file called `dc_workshop_log_XXXX_XX_XX.sh` (Use the four-digit year, two-digit month, and two digit day, e.g. `dc_workshop_log_2021_05_11.sh`)  
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> $ history | tail -n10 >> dc_workshop_log_XXXX_XX_XX.sh
> ```
>
> **NOTE:** Here we used the last 10 lines as an example, the number of lines may vary depending on the content of your `history` log.
> </details>

You may have noticed that your history contains the `history` command itself. To remove this redundancy from our log, let's use the `nano` text editor to fix the file:  

```bash
$ nano dc_workshop_log_XXXX_XX_XX.sh
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
# 2021/05/11
# Created sample directories for the workshop  
```

Next, remove any lines of the history that are not relevant by navigating to those lines and using your delete key. Save your file and close `nano`.

Your file should look something like this: 

```
# 2021/05/11
# Created sample directories for the workshop

mkdir dc_workshop/
mkdir dc_workshop/docs/
mkdir dc_workshop/data/
mkdir dc_workshop/results/
```

If you keep this file up to date, you can use it to re-do your work on your project if something happens to your results files. To demonstrate how this works, first delete your `dc_workshop/` directory and all of its subdirectories. Look at your directory contents to verify the directory is gone. 

```bash
$ rm -r dc_workshop
$ ls
shell_data	dc_workshop_log_XXXX_XX_XX.sh
```

Then run your workshop log file as a bash script. You should see the `dc_workshop/` directory and all of its subdirectories reappear. 

```bash
$ bash dc_workshop_log_2017_10_27.sh
$ ls
shell_data	dc_workshop dc_workshop_log_XXXX_XX_XX.sh
```

It's important that we keep our workshop log file outside of our `dc_workshop/` directory if we want to use it to recreate our work. It's also important for us to keep it up to date by regularly updating with the commands that we used to generate our results files.

---

## References

The workflow used here doesn't have to be the exact same as what is used in your team, particularly when we start working with the `slurm` management system in later workshops. You can read the references below for some ideas about what is worth documenting.

1. [Noble WS (2009) A Quick Guide to Organizing Computational Biology Projects. PLoS Comput Biol 5(7): e1000424.](https://doi.org/10.1371/journal.pcbi.1000424)
1. [Sandve GK, Nekrutenko A, Taylor J, Hovig E (2013) Ten Simple Rules for Reproducible Computational Research. PLoS Comput Biol 9(10): e1003285.](https://doi.org/10.1371/journal.pcbi.1003285)

---

[Next lesson](03-intro-to-nesi.md)
