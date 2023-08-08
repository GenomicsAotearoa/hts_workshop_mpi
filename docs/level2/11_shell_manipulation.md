# Manipulating files in the shell

* Teaching: 30 minutes
* Exercises: 30 minutes

#### Objectives

* Search the contents of basic text files for specific strings.
* Copy and remove directories, create tarballs and compress files.

#### Keypoints

* You can view search file contents within a `less` instance, or using the `grep` command.
* The commands `cp` and `rm` can be applied to directories, with the correct parameters.
* You can stitch multiple files together using the `tar` command. This is useful when you have to upload or download a lot of files in a single command.
* You can compress files using the `gzip` and `zip` commands, and decompress them using the `gunzip` and `unzip` commands.
* You can view file permissions using `ls -l` and change permissions using `chmod`.

---

## Contents

1. [Searching files using `less`](#searching-files-using-less-refresher)
1. [Searching files using `grep`](#searching-files-using-grep)
1. [Copying and removing folders of files](#copying-and-removing-folders-of-files)
1. [Grouping files into tarballs](#grouping-files-into-tarballs)
1. [Extracting files from a tarball archive](#extracting-files-from-a-tarball-archive)
1. [Compressing and uncompressing large files](#compressing-and-uncompressing-large-files)
1. [Setting file permissions](#setting-file-permissions)

---

## Searching files using `less` (refresher)

In previous training we have used the `less` command to peak inside text files and read their contents. This is great for small files where we can see the full content in a single terminal but sometimes we need to search larger documents for specific keywords or phrases. If you need a reminder of how this is performed, expand the refresher section below.

> <details>
> <summary>Refresher for searching with `less`</summary>
> 
> Some navigation commands in `less`:
> 
> |Key|Action|
> |:---:|:---|
> |<kbd>↓</kbd>|Go forward one line|
> |<kbd>↑</kbd>|Go back one line|
> |<kbd>Space</kbd>|Go forward one page|
> |<kbd>b</kbd>|Go back one page|
> |<kbd>g</kbd>|Return to the beginning of the file|
> |<kbd>G</kbd>|Jump to the end of the file|
> |<kbd>q</kbd>|Quit|
> 
> We can also use `less` to  search through files. Use the <kbd>/<kbd> key to begin a search. Enter the word you would like to search for and press <kbd>Enter<kbd>. The screen will jump to the next location where that word is found. 
> 
> If you hit <kbd>/<kbd> then <kbd>Enter<kbd> again, `less` will  repeat the previous search starting from the current location. This means that if you are near the end of a file and search returns no matches the term may be in your document, just prior to the location you searched from. You can use <kbd>?<kbd> instead of <kbd>/<kbd> if you wish to search backwards through the document.
> 
> </details>

---

## Searching files using `grep`

As powerful as `less` can be, at some point it becomes impractical to use it to screen documents as even with the search tools we are retrieving too many search hits, or covering too much content to reasonably interpret it by eye. For such situations, we can use another command line tool to search through documents without opening them and report the matches directory to the terminal. This tool is called `grep`.

`grep` is a command-line utility for searching plain-text files for lines matching a specific set of characters (sometimes called a string) or a particular pattern. We are going to work with one of the fastq files to practice using `grep` and to demonstrate some of the tool features. For these exercises we will searching some fastq files based on theri sequence content. As we are probably all aware, the four nucleotides that appear in DNA are abbreviated `A`, `C`, `T` and `G`. According to the IUPAC (International Union of Pure and Applied Chemistry) code, unknown nucleotides are represented with the letter `N` (a**N**y base). An `N` appearing in a sequencing file represents a position where the sequencing machine was not able to confidently determine the nucleotide in that position.

We'll search for strings inside of our fastq files. Let's first make sure we are in the correct directory:

```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/level2/shell_data/
```

Suppose we want to see how many reads in our file have really bad segments containing 10 consecutive unknown nucleotides (`N`s).

Let's search for the string `NNNNNNNNNN` in the `SRR098026.fastq` file:

```bash
$ grep NNNNNNNNNN SRR098026.fastq
```

This command returns a lot of output to the terminal. Every single line in the `SRR098026.fastq` file that contains at least 10 consecutive `N`s is printed to the terminal, regardless of how long or short the file is. We may be interested not only in the actual sequence which contains this string, but in the name (or identifier) of that sequence. We discussed in a previous lesson that the identifier line immediately precedes the nucleotide sequence for each read in a fastq file. We may also want to inspect the quality scores associated with each of these reads. To get all of this information, we will return the line immediately before each match and the two lines immediately after each match (see [this description of the fastq format](../docs/fastq_format.md) if you are unsure why we are using these numbers).

We can use the `-B` argument for grep to return a specific number of lines before each match. The `-A` argument returns a specific number of lines after each matching line. Here we want the line *before* and the two lines *after* each matching line, so we add `-B1 -A2` to our grep command:

```bash
$ grep -B1 -A2 NNNNNNNNNN SRR098026.fastq
```

This will also return a lot of text, but the bottom four lines should look like:

```
@SRR098026.177 HWUSI-EAS1599_1:2:1:1:2025 length=35
CNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+SRR098026.177 HWUSI-EAS1599_1:2:1:1:2025 length=35
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
```

> **Exercise**
>
> 1. Search for the sequence `GNATNACCACTTCC` in the `SRR098026.fastq` file. Have your search return all matching lines and the name (or identifier) for each sequence that contains a match.
> 
> <details>
> <summary>Solution</summary>
>  
> ```bash
> grep -B1 GNATNACCACTTCC SRR098026.fastq
> ```
> ```
> @SRR098026.245 HWUSI-EAS1599_1:2:1:2:801 length=35
> GNATNACCACTTCCAGTGCTGANNNNNNNGGGATG
> ```
> </details>
>
> 2. Search for the sequence `AAGTT` in both FASTQ files. Have your search return all matching lines and the name (or identifier) for each sequence that contains a match.
>
> <details>
> <summary>Solution</summary>
>  
> ```bash
> grep -B1 AAGTT *.fastq
> ```
> ```
> SRR097977.fastq-@SRR097977.11 209DTAAXX_Lenski2_1_7:8:3:247:351 length=36
> SRR097977.fastq:GATTGCTTTAATGAAAAAGTCATATAAGTTGCCATG
> --
> SRR097977.fastq-@SRR097977.67 209DTAAXX_Lenski2_1_7:8:3:544:566 length=36
> SRR097977.fastq:TTGTCCACGCTTTTCTATGTAAAGTTTATTTGCTTT
> --
> SRR097977.fastq-@SRR097977.68 209DTAAXX_Lenski2_1_7:8:3:724:110 length=36
> SRR097977.fastq:TGAAGCCTGCTTTTTTATACTAAGTTTGCATTATAA
> --
> SRR097977.fastq-@SRR097977.80 209DTAAXX_Lenski2_1_7:8:3:258:281 length=36
> SRR097977.fastq:GTGGCGCTGCTGCATAAGTTGGGTTATCAGGTCGTT
> --
> SRR097977.fastq-@SRR097977.92 209DTAAXX_Lenski2_1_7:8:3:353:318 length=36
> SRR097977.fastq:GGCAAAATGGTCCTCCAGCCAGGCCAGAAGCAAGTT
> --
> SRR097977.fastq-@SRR097977.139 209DTAAXX_Lenski2_1_7:8:3:703:655 length=36
> SRR097977.fastq:TTTATTTGTAAAGTTTTGTTGAAATAAGGGTTGTAA
> --
> SRR097977.fastq-@SRR097977.238 209DTAAXX_Lenski2_1_7:8:3:592:919 length=36
> SRR097977.fastq:TTCTTACCATCCTGAAGTTTTTTCATCTTCCCTGAT
> --
> SRR098026.fastq-@SRR098026.158 HWUSI-EAS1599_1:2:1:1:1505 length=35
> SRR098026.fastq:GNNNNNNNNCAAAGTTGATCNNNNNNNNNTGTGCG
> ```
> </details>

---

## Copying and removing folders of files

In the level 1 training we covered how to copy files and remove empty folders (directories), but it was noted that the copy operation would not work on folders and the `rmdir` command did not work on folders containing files.

If you are in a position where you need to copy a directory of files to a new location, such as when creating backups this can be achieved by adding the `-r` (recursive) flag to the copy (`cp`) command:

```bash
$ cp -r shell_data/ shell_data_backup/
```

In the level 1 training the way to remove the `shell_data_backup/` folder would be to navigate into the directory, remove all files, then navigate out of the directory and remove it with the `rmdir` command. If you are careful this can be compressed into a single command, by adding the `-r` flag to the remove (`rm`) command.

```bash
$ rm -r shell_data_backup/
```

This is subject to the usual `bash` warning that contents removed in this way is not recoverable. If you try to perform this operation on a folder with contents which have had their write permissions removed you will receive a confirmation prompt for each file in the directory. This may or may no be a problem depending on the number of files that this applies to.

> **Exercise**
>
> Use the `rm` help content to find the way to delete a folder of files bypassign the need to approach each file for detetion.
> 
> <details>
> <summary>Solution</summary>
>  
> ```bash
> $ rm -rf shell_data_backup/
> ```
> </details>

---

## Grouping files into tarballs

When creating backup directories for files, it may the be the case that we do not expect to need the contents of the backup frequently. Given that NeSI imposes both disk usage and file count limitations on our directories it is a good idea to reduce the data footprint of our backup files as much as possible. 

For this tutorial, navigate to the `tarball_extractions/` folder.

```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/level2/tarball_extractions/
```

The first part of this that we can adddress is reducing the number of files existing on the hard drive. This can be quite easily acheived using the `tar` command which glues a number of individual input files into a single 'tarball' file. This file is not compressed in any way, but represents the contents of multiple input files in a single entry. Manipulating existing tarball files can be quite complicated but the basic premise for creating one is simple.

```bash
$ tar -cf shell_data.tar shell_data/
```

This command is comprised of four elements, described blow:

|Term|Meaning|
|:---:|:---|
|`tar`|The `tar` command, invoked from the terminal|
|`-cf`|The run flags for the `tar` command<br>In this instance the flags are to create (`-c`) and new file (`-f`)|
|`shell_data.tar`|The name of the tarball file to be created|
|`shell_data/`|A folder or list of files to be added to the tarball file|

If you do not want to add all of the contents of a directory you can specify individual files to be added by changing the command to individual files rather than a directory:

```bash
$ tar -cf other_shell_data.tar shell_data/SRR097977.fastq shell_data/SRR098026.fastq
```

We can then check the contents of the tarball files using the `tar` command with the list (`-t`) flag to view the contents of the tarball.

```bash
$ tar -tf shell_data.tar
$ tar -tf other_shell_data.tar
```

As you can see here, there are differences in the contents of each tarball as the `shell_data/` folder actually contains a hidden folder which gets picked up with the directory input method but is ignored when we explicitly state target files to be added to the tarball.

--

## Extracting files from a tarball archive

Once you have obtained a tarball file, you will probably need to unpack the contents at some point. In cases where we simply want to extract the complete contents of the tarball file this is easily achieved with the extract (`-x`) flag.

Navigate to the `tarball_extractions/` folder and then perform the following operation:

```bash
$ tar -xf shell_data.tar
```

Using the `ls` command with the 'all files' flag (`-a`) you will be able to see that both of the fastq files and the hidden directory have been extracted to a new folder named `shell_data/`. We were given no choice in which directory the output files were sent to as the file paths are encoded in the tarball file from when it was created. When we perform the extract operation the file paths of the tarball built in the directory from which we invoke the command. If folders with the appropriate names exist in the directory the tarball files will be extracted into these directories, and if the correct directory does not exist it will be created by `tar`.

Sometimes this is not the behaviour we want to see - it may be the case that we only want certain tarball files to be extracted, or we wish to send them to a different output path. Before we demonstrate these alternate behaviours, deleted the `shell_data/` directory and its contents.

> <details>
> <summary>How do I delete a folder again?</summary>
>  
> ```bash
> $ rm -r shell_data/
> ```
> </details>

Firstly, lets say that we want to extract just the fastq files from the tarball, but still into a new directory named `shell_data/`. This can be achieved in exactly the same way as creating a tarball from a specific list of files, just swapping the `-c` flag with the `-x` flag:

```bash
$ tar -xf shell_data.tar shell_data/SRR097977.fastq shell_data/SRR098026.fastq
```

You can run a quick `ls` to confirm that only the desired files have been extracted.

---

## Compressing and uncompressing large files

We have now seen how we can use `tar` to glue multiple files together into a single output file, and how we can be either broad or selective in which files are added to or extracted from the tarball. This does not address the issue of file size, however. To reduce the file space required by the output file we need to use a different tool. There are two commonly used command line compression tools, `gzip` and `zip`. We will only be covering `gzip` in this tutorial.

Each of these tools applies a specific compression algorithm to your file contents to create a smaller, less readable, verison of the file. There is no particular preference to use either of these tools on NeSI other than the fact that the `zip` format is used on Windows devices and may be your default compression tool if you are uploading or downloading compressed files to and from NeSI. Alternatively, if you are sourcing files from online the `gzip` compression standard is very popular and you may find that you need to use this tool for reference files you download from resources such as the NCBI Sequence Read Archive.

>Note: The file encoding schema used by `gzip` is not compatible with the `zip` standard used by default on your Windows device. If you are uploading compressed file contents from your work computer you may need to use the `zip` and `unzip` tools to compress and extract these files instead of `gzip`.

Navigate to the `file_compression/` folder for this set of exercises.

```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/level2/file_compression/
```

The file we have there is the results of some publicly available MinION sequencing of a *Helicobacter pylori* isolate. The sequences were obtained from the NCBI Sequence Read Archive ([SRR23286833](https://www.ncbi.nlm.nih.gov/sra/SRX19230161[accn])).

We can quickly check the size of this downloaded file using the `ls` command with the `-s` (file size) and `-h` (human readable) like so:

```bash
$ ls -sh
```

```
128M SRR23286833.fastq
```

To compress this file, we simply need to run the `gzip` command and provide it with the target file to be compressed:

```bash
$ gzip SRR23286833.fastq
$ ls -sh
```

```
59M SRR23286833.fastq.gz
```

So we can see that this has reduced the file size sequences by about half. That's a great saving in disk space but what has happened to the original file? If we check the manual for `gzip`, you can see the following text:

> DESCRIPTION
>
> Gzip reduces the size of the named files using Lempel-Ziv coding (LZ77). Whenever possible, **_each file is replaced_** by one with the extension .gz, while keeping the same ownership modes, access and modification times.

Is this a problem? It can be. When we run `gzip` it is operating under the assumption that we do not need the original file anymore. This is probably true if we were to compress a tarball file, as producing a tarball leaves the original files in place. However, what if we just wanted to compress the sequence files but then continue to work with them?

When we are working with [fastq files](../docs/fastq_format.md) many of the standard bioinformatic tools can read compressed sequence information natively, so you can pass them the raw fastq file or a compressed representation of the file without issue. When we are working with other file formats though, like fasta or the outputs of a `BLASTn` annotation our other tools might not be able to read the compressed file format. In these cases, we might need to uncompress the file to recover the original, plain text output.

To start with, we can decompress the `gz` file using `gzip` with a particular flag or with the tool `gunzip`, which is the partner to `gzip`. Both options give the same result:

```bash
$ gzip --decompress SRR23286833.fastq.gz
# OR
$ gunzip SRR23286833.fastq.gz
$ ls -sh
```

```
128M SRR23286833.fastq
```

In either case, we get the same output file.

> <details>
> <summary>How do I compress a file and leave the original in place?</summary>
> 
> In some cases you may wish to compress a file but leave the original untouched. The most common way to achieve this is to run `gzip` with the `-c` flag, which sends the results of the compression algorim to the *stdout* channel rather than to a file. We will explore the *stdout* channel in a [later session](./13_shell_redirection.md), but for now this is the command to create a compressed file alongside the original:
> 
> ```bash
> $ gzip -c SRR23286833.fastq > SRR23286833.fastq.gz
> ```
> </details>

---

## Setting file permissions

If it's that easy to permanently delete a file, how can be put some checks in place to prevent it from happening accidentally? The answer to this is through file permissions. In your `shell_data/` folder, run the following command:

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

![Permissions breakdown](../img/level2_11_rwx_figure.svg)

We're going to concentrate on the three positions that deal with your permissions (as the file owner), which will have the values `rw-`. This shows that you have permission to read (`r`) the file as well as edit the contents (write, `w`). The third position is set to `-`, indicating that you don't have permission to carry out the ability encoded by that space. This is the space where the ability to execute the file (`x`) is set, we'll talk more about this in a later lesson.

Our goal for now is to change permissions on this file so that you no longer have `w` or write permissions. We can do this using the `chmod` ("change mode") command and subtracting (`-`) the write permission `-w`. 

```bash
$ chmod -w SRR097977.fastq SRR098026.fastq
$ ls -l
```

```
drwxrws---+ 2 dwaite comm00008  4096 Feb 24 15:50 backup
drwxrws---+ 2 dwaite comm00008  4096 Feb 24 15:50 other_backup
-r--r--r--+ 1 dwaite comm00008 47552 Jun 25  2021 SRR097977.fastq
-r--r--r--+ 1 dwaite comm00008 43332 Jun 25  2021 SRR098026.fastq
```

You can see that the write permissions for these files have been removed. For other users, they will be completely unable to remove or modify these files. For you, *as the file owner* it is still possible to delete the file, but you will first get a confirmation asking if we wish to proceed:

```bash
$ rm SRR098026.fastq
```

```
rm: remove write-protected regular file ‘SRR098026.fastq’?
```

---
