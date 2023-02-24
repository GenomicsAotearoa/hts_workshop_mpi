# Manipulating files in the shell

* Teaching: 30 minutes
* Exercises: 15 minutes

#### Objectives

* View and search within files for basic pieces of text.
* Copy, move, and rename filesm and create/remove directories.
* Make a file read only.
* Use the `history` command to view and repeat recently used commands.

#### Keypoints

* You can view file contents using `less`, `cat`, `head` or `tail`.
* The commands `cp`, `mv`, and `mkdir` are useful for manipulating existing files and creating new directories.
* You can view file permissions using `ls -l` and change permissions using `chmod`.
* The `history` command and the up arrow on your keyboard can be used to repeat recently used commands.

---

## Contents

1. [Viewing the contents of files](#viewing-the-contents-of-files)
1. [Basic file manipulation](#basic-file-manipulation)
1. [Setting file permissions](#setting-file-permissions)
1. [Command history](#command-history)

---

## Viewing the contents of files

From the previous exercises we know how to move around the file system of NeSI, but how do we look at the contents of files? One way to examine a file is to print out all of the contents using the program `cat`.

To look at a text file, navigate to the `shell_data/` directory in your training folder and try to run the following command:

```bash
$ cat SRR097977.fastq
```

That wasn't very helpful. It printed the full file content to the screen. It's a very small file, as bioinformatic files go, but still far to much to scan by eye. For smaller files, `cat` is a terrific tool but when the file is really big, it can be annoying to use.

Fortunately there is another handy tool to read large files in a more manageable way. The command `less` can be used to open a file and navigate through it line by line. Enter the following command:

```bash
$ less SRR097977.fastq
```

This will load the content of the file into your terminal, but rather than print every line instantly only those that can fit on one page are shown. Since this is a tool designed to run from the command line only we generally need to navigate using the keyboard. Some of the commonly used navigation commands are:

|Key|Action|
|:---:|:---|
|<kbd>↓</kbd>|Go forward one line|
|<kbd>↑</kbd>|Go back one line|
|<kbd>Space</kbd>|Go forward one page|
|<kbd>b</kbd>|Go back one page|
|<kbd>g</kbd>|Return to the beginning of the file|
|<kbd>G</kbd>|Jump to the end of the file|
|<kbd>q</kbd>|Quit|

That said, you can also scroll forward and back theough the file using the mouse scroll wheel.

`less` can also be used to search through files. Use the <kbd>/</kbd> key to begin a search. Enter the word you would like to search for and press <kbd>Enter</kbd>. The screen will jump to the next location where that word is found. You can seach for the next word by pressing <kbd>/</kbd> repeatedly. Each time, `less` searches from the current location forward. If you need to go back one entry, use <kbd>?</kbd>.

> ### Exercise
>
> As an example, let's search forward for the sequence `TTTTT` in our file What are the next three nucleotides (characters) after the first instance of this sequence?
> 
> <details>
> <summary>Solution</summary
>
> `CAC`
> </details>

Sometimes we want to strike a balance between `cat` and `less`. We need to see a bit of the file, but we don't want to look at it line by line. This is uaually when we need to process a file in some way, and we just need to remind outselves how it's formatted.

The commands `head` and `tail` are for this task. They let you look at the beginning and end of a file, respectively.

```bash
$ head SRR098026.fastq
@SRR098026.1 HWUSI-EAS1599_1:2:1:0:968 length=35
NNNNNNNNNNNNNNNNCNNNNNNNNNNNNNNNNNN
+SRR098026.1 HWUSI-EAS1599_1:2:1:0:968 length=35
!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!
@SRR098026.2 HWUSI-EAS1599_1:2:1:0:312 length=35
NNNNNNNNNNNNNNNNANNNNNNNNNNNNNNNNNN
+SRR098026.2 HWUSI-EAS1599_1:2:1:0:312 length=35
!!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!
@SRR098026.3 HWUSI-EAS1599_1:2:1:0:570 length=35
NNNNNNNNNNNNNNNNANNNNNNNNNNNNNNNNNN
```

```bash
$ tail SRR098026.fastq
+SRR098026.247 HWUSI-EAS1599_1:2:1:2:1311 length=35
#!##!#################!!!!!!!######
@SRR098026.248 HWUSI-EAS1599_1:2:1:2:118 length=35
GNTGNGGTCATCATACGCGCCCNNNNNNNGGCATG
+SRR098026.248 HWUSI-EAS1599_1:2:1:2:118 length=35
B!;?!A=5922:##########!!!!!!!######
@SRR098026.249 HWUSI-EAS1599_1:2:1:2:1057 length=35
CNCTNTATGCGTACGGCAGTGANNNNNNNGGAGAT
+SRR098026.249 HWUSI-EAS1599_1:2:1:2:1057 length=35
A!@B!BBB@ABAB#########!!!!!!!######
```

By default the first/last 10 lines are printed. This can be changed by adding the `-n` option to the command to change the number. 

```bash
$ head -n1 SRR098026.fastq
@SRR098026.1 HWUSI-EAS1599_1:2:1:0:968 length=35
```

```bash
$ tail -n1 SRR098026.fastq
A!@B!BBB@ABAB#########!!!!!!!######
```

> <details>
> <summary>Details of the FASTQ format</summary
>
> In this session we have been peaking inside fastq files, although we haven't actually discussed what these files contain. Although it looks complicated, the format is actually pretty easy to understand.
> The [fastq format](https://en.wikipedia.org/wiki/FASTQ_format) is used for writing nucleic acid sequences with a score for the confidence in each nucleotide in the sequence. This is typically the result you get back from a sequencing platform. The file is broken into a repeating pattern of 4 lines, where each line describes a single sequence read. Within these four lines:
>
> |Line|Description|
> |:---:|:---|
> |1|Is the name of the sequence, and optionally some metadata about the read. This line always begins with a `@` symbol, similar to how a fasta file begins a sequence name with `>`.|
> |2|The nucleic acid sequence.|
> |3|A placeholder line which divides the nucleic acid sequence and quality information. In some older applications this line will contain a `+` symbol followed by the sequence name from Line 1, although newer tools tend to just use the `+` alone to save file space.|
> |4|A string of characters which represent the quality scores. This string  must be the same length as Line 2.|
>
> To make it easy to understand, running the `head` command with a `-n` parameter of 4 we can view the first complete read in one of the files.
>
> ```bash
> $ head -n4 SRR097977.fastq
> @SRR097977.1 209DTAAXX_Lenski2_1_7:8:3:710:178 length=36
> TATTCTGCCATAATGAAATTCGCCACTTGTTAGTGT
> +SRR097977.1 209DTAAXX_Lenski2_1_7:8:3:710:178 length=36
> CCCCCCCCCCCCCCC>CCCCC7CCCCCCACA?5A5<
> ```
>
> Line 4 shows the quality for each nucleotide in the read. Quality is interpreted as the probability of an incorrect base call (e.g. 1 in 10) or, equivalently, the base call accuracy (e.g. 90%). These probability values are the results from the base calling algorithm and dependent on how much signal was captured for the base incorporation.  To make it possible to line up each individual nucleotide with its quality score, the numerical score is converted into a code where each individual character represents the numerical quality score for an individual nucleotide. In the line above, the quality score line is `CCCCCCCCCCCCCCC>CCCCC7CCCCCCACA?5A5<`.
>
> The numerical value assigned to each of these characters depends on the sequencing platform that generated the reads. The sequencing machine used to generate our data uses `PHRED+33` encoding, which is what is used in all modern Illumina and Nanopore sequencing, as well as Sanger sequencing. Each character is assigned a quality score between 0 and 42 as shown in the chart below.
>
> ```
> Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJK
>                   |         |         |         |         |
> Quality score:    0........10........20........30........40..
> ```
>
> If you line it up, you can see that the `C` character has a qualtiy score of 34, which is the most common value in this sequence. The probability of a position with a Q of 34 can be calculated as:
>
> $$ P_{incorrect} = 10^{\frac{Q}{-10}} = 0.000398107 $$
>
> Which can be converted to the probabilty that the read is correct by
>
> $$ P_{correct} = 1 - P_{incorrect} = 0.9996 $$
>
> Therefore, for any of these characters the base call accuracy is 99.96%. You will almost never need to calculate the per-position score, but it is handy to know the following numbers for when we are quality filtering:
>
> |Q|Accuracy|
> |:---:|:---:|
> |10|90.00%|
> |20|99.00%|
> |30|99.90%|
> |40|99.99%|

---

## Basic file manipulation

We now know how to read files through the command line, but before we worry about any complex bioinformatic work the most basic thing we need to be able to do on the command line, other than navigating directories, is moving, copying, and deleting files and directories. These are operations which you probably perform daily on your desktop computer and these all exist on the command line as well.

There are commands available on the command line which do all of these things. These typically have short names which are derived from the word that represent:

|Command|Action|
|:---:|:---|
|`mkdir`|Make a new directory.|
|`rmdir`|Remove (delete) an empty directory.|
|`cp`|Copy a file, either to a new location or into a new file.|
|`mv`|Move a file from one location to another. Can also be used to rename files by 'moving' them to a file with a different name.|
|`rm`|Remove a file. This is how you delete.|

### Creating and removing directories

The `mkdir` ("make directory") command is used to make a new directory. Enter `mkdir` followed by a space, then the directory name you want to create:

```bash
$ mkdir backup/
```

If you want to create multiple directories at once, you can specify multiple names:

```bash
$ mkdir other_backup/ another_backup/
```

As long as these directories are empty, they can be removed with the `rmdir` ("remove directory") command:

```bash
$ rmdir another_backup/
```

If you try to remove a directory with files in it, you will recieve an error and the directory will remain intact.

### Copying and moving files

The `cp` ("copy") and `mv` ("move") commands are mostly identical in how they work. Each command can be used in one of two ways. We can either copy/move a file from one name to another, or we can copy/move them into a new directory without changing the file names.

```bash
$ cp SRR097977.fastq SRR097977.fq_backup
$ mv SRR097977.fq_backup SRR097977.fq_bkup
```

In this pair on commands, we first create a copy of the file `SRR097977.fastq` with the name `SRR097977.fq_backup`. These are identical in their content. We then move/rename the `SRR097977.fq_backup` file (effectively, moving the contents of the file into a new file) to a different file named `SRR097977.fq_bkup`. If you run the `ls` command you will see that `SRR097977.fastq` is still present, as it was only copied, but `SRR097977.fq_backup` no longer exists.

>**Note:** we can also use the `mv` command to rename directories. However, the `cp` command by default does not work for directories, it must be invoked with a specific parameter to copy directory and it's contents.

The other way we can use these commands is to copy/move one or more files into a different directory. This can be handy when creating backups of data we do not wish to risk losing (copy), or when we want to organise data into different folders to make navigation easier (move).

```bash
$ cp SRR097977.fastq SRR098026.fastq backup/
$ mv SRR097977.fq_bkup other_backup/
```

In these cases, we either create duplicates of the target files in a new location, or move an existing file into a new location.

### Removing files

To remove a file, it's super easy. Just use the `rm` ("remove") command:

```bash
$ rm other_backup/SRR097977.fq_bkup
```

Boom. It's gone. There is no Recycle Bin on the command line, there is no way to get that file back...

We need to be really careful with the `rm` command - be very careful in which files you are removing. Because of this, the `rmdir` command only removes empty directories and by default the `rm` command will not delete directories either. At this level, we will not expand upon this further.

---

## Setting file permissions

If it's that easy to permanently delete a file, how can be put some checks in place to prevent it from happening accidentally. The answer to this is through file permissions. In your `shell_data/` folder, run the following command:

```bash
$ ls -l
```

This should return a view similar to:

```
drwxrws---+ 2 dwaite comm00008  4096 Feb 24 15:50 backup
drwxrws---+ 2 dwaite comm00008  4096 Feb 24 15:50 other_backup
-rw-rw-r--+ 1 dwaite comm00008 47552 Jun 25  2021 SRR097977.fastq
-rw-rw-r--+ 1 dwaite comm00008 43332 Jun 25  2021 SRR098026.fastq
```

What we interested in is the first part of the output, the strings of characters which look like `-rw-rw-r--+`. This represents the current permission state of the file. If we ignore the first and last characters, what we are left with are 9 characters, which represent 3 file permissions for 3 different user groups, like so:

![Permissions breakdown](../img/level1_02_rwx_figure.svg)

We're going to concentrate on the three positions that deal with your permissions (as the file owner), which will have the values `rw-`. The `r` means that you have permission to read the file, the `w` indicates that you have permission to write to (i.e. modify) the file, and the third position is set to `-`, indicating that you don't have permission to carry out the ability encoded by that space. This is the space where `x` or executable ability is stored, we'll talk more about this in a later lesson.

Our goal for now is to change permissions on this file so that you no longer have `w` or write permissions. We can do this using the `chmod` ("change mode") command and subtracting (`-`) the write permission `-w`. 

```bash
$ chmod -w SRR097977.fastq SRR098026.fastq
$ ls -l 
drwxrws---+ 2 dwaite comm00008  4096 Feb 24 15:50 backup
drwxrws---+ 2 dwaite comm00008  4096 Feb 24 15:50 other_backup
-r--r--r--+ 1 dwaite comm00008 47552 Jun 25  2021 SRR097977.fastq
-r--r--r--+ 1 dwaite comm00008 43332 Jun 25  2021 SRR098026.fastq
```

You can see that the write permissions for these files have been removed. For other users, they will not be completely unable to remove or modify these files. For you, *as the file owner* it is still possible to delete the file, but now we get a confirmation asking if we wish to proceed:

```bash
$ rm SRR098026.fastq
rm: remove write-protected regular file ‘SRR098026.fastq’?
```

---

## Command history

We've done a lot on the command line, but how do we keep track of everything so far? What if we forget something we've done recently, and want to check exactly what we did?

You can view previous commands using the up arrow on your keyboard to go back through your recent commands. Likewise, the down arrow takes you forward in the command history. If what you're looking for is only a few commands ago this is the best way to see the information. IF you are looking for something from a long time aho, you can view a list of your last ~1,000 commands with the `history` command:

```bash
$ history
```

This will print a numbered list of recent commands. You can reuse one of these commands directly by referring to the number of that command.

For example, if your history looked like this:

```
1053  2023-02-24 15:58:35 ll
1054  2023-02-24 16:01:35 ls
1055  2023-02-24 16:01:42 ls -sh
1056  2023-02-24 16:01:45 history
```

You could repeat command #1055 by entering:

```bash
$ !1055
```

This can be really useful when you are taking notes on a complicated set of commands you have run, or if you are trying to remember what you did last time you were logged into NeSI.

---
