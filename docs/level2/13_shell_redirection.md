# Working with redirection

* Teaching: 30 minutes
* Exercises: 30 minutes

#### Objectives

* Understand the basic idea or shell redirection and what can be done with it.
* Understand the differences between redirecting to a file and another command.
* Know that different command line tools can be linked together to produce complex operations.

#### Keypoints

* When a command on the shell is printing information to the terminal, we can capture and redirect that information to a new location or into a new command.
* When redirecting to a file, we can create a new file, overwrite an existing file, or append to an existing file.
* Each channel can be independently captured and directed to a new location. SUch locations include:

---

## Contents

1. [Redirecting output to a new file](#redirecting-output-to-a-new-file)
   1. [Applying redirection to loops](#applying-redirection-to-loops)
      1. [Using appending redirection](#using-appending-redirection)
      1. [Using variables in the output file name](#using-variables-in-the-output-file-name)
1. [Redirecting output to a different command](#redirecting-output-to-a-different-command)

---

## Redirecting output to a new file

In the previous session we were using the `grep` command allowed us to identify sequences in fastq files that matched particular patterns. All of these sequences were printed to our terminal screen which wasn't particularly helpful so what we're now going to do is learn how to capture that output in a more useful manner.

We can do this with something called **redirection**. The idea is that we are taking what would ordinarily be printed to the terminal screen and redirecting it to another location. In our case, we want to move the information printed on our terminal into a file so that we can look at it later.

Navigate to the `/nesi/project/nesi03181/phel/USERNAME/level2/redirection/` folder and we'll get started.

```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/level2/redirection/
```

The command for redirecting output to a file is <kbd>></kbd>. Let's revisit our example using the `${MOTIF}` variable to represent the sequence motif we want to capture and we'll redirect the output of the search into a new file.

```bash
$ MOTIF="NNNNNNNNNN"
$ grep ${MOTIF} SRR098026.fastq > my_file.txt
```

This might take a second or two to complete, but when it is done a quick `ls` should show you a new file called `my_file.txt`. This file did not exist before, but the <kbd>></kbd> symbol when performing redirection is understood by the computer to mean that a new file must be created with the name corresponding to the text which follows the <kbd>></kbd> character.

>Note: If you redirect into an existing file using this method, the old file contents will be replaced with the new values. You will lose the data permanently so make sure you are redirecting data to the correct place.

You might notice that we have used the file extension `.txt` in the output file rather than `.fastq`. This is largely semantics, in the command line environment the extensions to files generally have no meaning to the computer, they are for the user to help identify the expected contents of the file in question. In this instance, we are using a different extension as the `grep` command will only extract the lines (sequences) with the `${MOTIF}` value and not the sequence name or quality information so it will no longer be recognised as a `fastq` file.

> <details>
> <summary>Example of file content differences</summary>
>
> First 10 lines of `SRR098026.fastq`:
>
> ```
> @SRR098026.1 HWUSI-EAS1599_1:2:1:0:968 length=35
> NNNNNNNNNNNNNNNNCNNNNNNNNNNNNNNNNNN
> +SRR098026.1 HWUSI-EAS1599_1:2:1:0:968 length=35
> !!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!
> @SRR098026.2 HWUSI-EAS1599_1:2:1:0:312 length=35
> NNNNNNNNNNNNNNNNANNNNNNNNNNNNNNNNNN
> +SRR098026.2 HWUSI-EAS1599_1:2:1:0:312 length=35
> !!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!
> @SRR098026.3 HWUSI-EAS1599_1:2:1:0:570 length=35
> NNNNNNNNNNNNNNNNANNNNNNNNNNNNNNNNNN
> ```
>
> First 10 lines of `my_file.txt`:
>
> ```
> NNNNNNNNNNNNNNNNCNNNNNNNNNNNNNNNNNN
> NNNNNNNNNNNNNNNNANNNNNNNNNNNNNNNNNN
> NNNNNNNNNNNNNNNNANNNNNNNNNNNNNNNNNN
> NNNNNNNNNNNNNNNNANNNNNNNNNNNNNNNNNN
> NNNNNNNNNNNNNNNNGNNNNNNNNNNNNNNNNNN
> NNNNNNNNNNNNNNNNGNNNNNNNNNNNNNNNNNN
> NNNNNNNNNNNNNNNNGNNNNNNNNNNNNNNNNNN
> NNNNNNNNNNNNNNNNANNNNNNNNNNNNNNNNNN
> NNNNNNNNNNNNNNNNGNNNNNNNNNNNNNNNNNN
> NNNNNNNNNNNNNNNNGNNNNNNNNNNNNNNNNNN
> ```
> </details>

However, this is not a very useful way of reporting the bad sequences, so you're going to modify the `grep` command to extract enough information from the `SRR098026.fastq` file that the new file created using the redirection is a valid `fastq` file.

> **Exercise**
>
> Search through the `grep` manual to find the commands you will need to extract a complete sequence record from entries in the `SRR098026.fastq` file which match the `${MOTIF}`. Apply these to your `grep` command and verify if your output file is correct or not.
>
> *__Hint 1:__ There is a primer on the structure of a `fastq` file available in the course notes [here](../docs/fastq_format.md).*
>
> *__Hint 2:__ Remember that you can use the command `man grep` to get the full user manual for the `grep` command.*
>
> <details>
> <summary>Solution</summary>
>  
> ```bash
> grep -B1 -A2 ${MOTIF} SRR098026.fastq > my_file.fastq
> ```
> </details>

---

### Applying redirection to loops

We have now created a dynamic loop that can search a single, hardcoded file and create new `fastq` files with the contents which match the `${MOTIF}` value. what happens if we start to use a loop to specify changing values for `${MOTIF}`?

> **Exercise**
>
> Revise the for loop from the previous session ([notes here](./12_shell_variables.md#extending-our-loop-with-multiple-variables)) to loop through a number of values for `${MOTIF}` with the redirection parameters from your current work.
>
> After the loop completes and the `grep` search has run three times, what are the final contents of the output file and why?
> 
> 1. The results from all `grep` searches.
> 1. The results from the first `grep` search.
> 1. The results from the second `grep` search.
> 1. The results from the last `grep` search.
>
> <details>
> <summary>Solution</summary>
> 
> ```bash
> $ FILENAME=SRR098026.fastq
> $ for MOTIF in "NNNNNNNNNN" "GCTGGCGNNN" "TTTTTTTTTT";
> > do
> >     grep -B1 -A2 ${MOTIF} ${FILENAME} > my_file.fastq
> > done
> ```
> 
> Ony the results from the last search will be shown in `my_file.fastq`, as each iteration of the loop overwrites the contents of the file using the `>` operator.
>
> </details>

The is a problem because we may need to record the results of multiple searches, or even search multiple `fastq` files and report each set of results independently. There are two ways to deal with this and the best option depends on your requirements.

#### Using appending redirection

The first modification we can make is to replace the redirection operator `>` with the appending operator `>>`. This is a small change but now the results from each iteration of the loop will be added to the end of the output `fastq` file instead of replacing the current contents.

Write the following loop and see how the outputs differ.

```bash
$ FILENAME=SRR098026.fastq
$ for MOTIF in "NNNNNNNNNN" "GCTGGCGNNN" "TTTTTTTTTT";
> do
>     grep -B1 -A2 ${MOTIF} ${FILENAME} >> my_appended_file.fastq
> done
```

#### Using variables in the output file name

If you wanted an output where all of your search motifs were placed into a single file, that approach is fine. However, it is equally likely that we will not want this outcome and will instead want to create different output files for each motif searched. Try to incorporate this into your existing loop.

> **Exercise**
>
> Modify your loop above to create individual files according to the value of `${MOTIF}`. You can use either the overwriting (`>`) of appending (`>>`) operator.
>
> <details>
> <summary>Solution</summary>
> 
> ```bash
> $ FILENAME=SRR098026.fastq
> $ for MOTIF in "NNNNNNNNNN" "GCTGGCGNNN" "TTTTTTTTTT";
> > do
> >     grep -B1 -A2 ${MOTIF} ${FILENAME} > my_file.${MOTIF}.fastq
> > done
> ```
>
> </details>

---

## Redirecting output to a different command

Capturing outputs in a file is a very useful technique when working on the command line. However, when we are using a command line tool to process files we often run into the problem that there is no single tool powerful enough to do all the things we need in a single operation. This is largely part of the design philosophy behind command line tools - they should be very good at one thing, but only really do that one thing. There is then a form of redirection to pass the output of one tool directly into another, which allows us to chain together multiple different commands into a more powerful pipeline.

We will not cover all of the command line tools that exist on NeSI, but here are a few of the commonly used ones and when you might need to use them.

|Tool|Purpose|
|:---:|:---|
|`grep`|Search a file for a keyword, or set of keywords|
|`sed`|Search a file for a keyword, or set of keywords, and replace them with a new value|
|`cat`|Print the content of a file to the terminal|
|`head`|Print the first N lines of a file to the terminal|
|`tail`|Print the last N lines of a file to the terminal|
|`sort`|Sort the contents of a file alphabetically|
|`uniq`|Print the unique occurences in a stream of text entries*|
|`cut`|Extract a specific set of columns from a tab (default) or overwise delimited text file|

> <details>
> <summary>More about `uniq`</summary>
>
> The `uniq` command is basically as described, from a set of entries it will print out the unique occurences but it works by comparing the current entry in a sequence and testing whether or not it is the same as the last entry. This is not necessarily how you might expect it to work - a common assumption is that it will compare the current entry against *all previously observed entries* but this is not the case. To see the difference in action, consider a file with the following values:
>
> ```
> a
> a
> a
> b
> b
> c
> ```
>
> If we were to run this through `uniq` we would see the following:
>
> ```bash
> $ uniq my_file.txt
> ```
>
> ```
> a
> b
> c
> ```
>
> Which is pretty much as expected. However, if the content was instead:
>
> ```
> a
> a
> a
> b
> b
> c
> a
> ```
>
> Then we would instead see:
>
> ```bash
> $ uniq my_file.txt
> ```
>
> ```
> a
> b
> c
> a
> ```
>
> This is because `uniq` keeps no record that it has previously seen an `a` value. When it hits the last line of the file and encounters the `a` it only tests if that `a` matches the value of the previous line, which was a `c`. Since they do not match, the `a` is reported.
>
> </details>

In this exercise we are going to write a command that chains together three of these commands to produce a basic summary of a table of information. In the `redirection/` folder there is a table called `ncbi_viruses.txt` which is a non-exhaustive list of the virus names recorded in the [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy/) database. It's a big file (22,387 lines) so we're going to build a set of commands to perform a quick summary of some of the information it contains.

If you run a quick `less` (remembering <kdb>Q</kdb> to quit) or `head` command you will see that the contents look something like

```
taxid   Kingdom Phylum  Class   Order   Family  Genus   Species
2716741 Viruses Pisuviricota    Pisoniviricetes Picornavirales  Iflaviridae     Iflavirus       ACT flea iflavirus
1244521 Viruses Negarnaviricota Ellioviricetes  Bunyavirales    Hantaviridae    Orthohantavirus ANAJ Hantavirus
1482734 Viruses Pisuviricota    Pisoniviricetes Picornavirales  Picornaviridae  Aalivirus       Aalivirus A
2320189 Viruses Uroviricota     Caudoviricetes  unclassified    Autographiviridae       Aarhusvirus     Aarhusvirus dagda
```

The rendering of the contents is a bit uneven, but each text value is separated by a <kdb>Tab</kdb> character so if you copied this data into Excel you would see the values as a table:

|taxid|Kingdom|Phylum|Class|Order|Family|Genus|Species|
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|2716741|Viruses|Pisuviricota|Pisoniviricetes|Picornavirales|Iflaviridae|Iflavirus|ACT flea iflavirus|
|1244521|Viruses|Negarnaviricota|Ellioviricetes|Bunyavirales|Hantaviridae|Orthohantavirus|ANAJ Hantavirus|
|1482734|Viruses|Pisuviricota|Pisoniviricetes|Picornavirales|Picornaviridae|Aalivirus|Aalivirus A|
|2320189|Viruses|Uroviricota|Caudoviricetes|unclassified|Autographiviridae|Aarhusvirus|Aarhusvirus dagda|

Our aim here is to use the commands mentioned above to perform a quick tally of the genera within the family *Potyviridae*. This will require us to perform the following steps:

1. Select only the rows which contain the value 'Potyviridae'.
1. Extract the values from the column corresponding to 'Genus'.
1. Identify the unique genera within the column.

> **Exercise**
>
> Before we proceed have a quick look through the tools listed above and try to identify which tools can be used for each step.
>
> <details>
> <summary>Solution</summary>
>
> 1. `grep`
> 1. `cut`
> 1. `uniq`, combined with `sort` (see notes on `uniq` above)
>
> </details>

When we redirect from one tool to another, we use the pipe character (`|`). This is probably not a key on your keyboard you use very much, so let's take a minute to find that key. For the standard QWERTY keyboard layout, the `|` character can be found using the key combination <kbd>Shift</kbd>+<kbd>\</kbd>. Rather than direct to a file location, this now captures the printed output and directs it to an input channel of the command that follows the pipe operator.

For example, if we are passing the output of the `grep` command to find lines in the table with the term 'Potyviridae' and want to pass it into the `cut` command, it will look like:

```bash
> grep Potyviridae ncbi_viruses.txt | cut
```

From there we can provide `cut` with the parameters we need and then either direct its output to a file (`>`) or to another command (`|`). In this instance we only have a single parameter to provide `cut` and that is which column(s) we want it to extract from the input data. This is done with the `-f` parameter which takes the column position (not label) that we want to retain. You can provide multiple values to `-f`, separated with a comma, but in this case we will only need the one.

```bash
$ grep Potyviridae ncbi_viruses.txt | cut -f7 | sort | uniq
```

And that's most of the command we need. The initial file is filtered for lines that contain the text 'Potyviridae', those lines are passed into the `cut` command which extracts just the information in the genus column. The results are then alphabetically sorted and `uniq` reports the unique names found.

The only thing we need to do to complete the exercise is to report how many of each genus were found.

> **Exercise**
>
> Search the manual for `uniq` and find the parameter required to print out the number of observations for each unique entry.
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> $ grep Potyviridae ncbi_viruses.txt | cut -f7 | sort | uniq -c
> ```
>
> ```
>   2 Arepavirus
>   1 Bevemovirus
>   1 Brambyvirus
>   7 Bymovirus
>   1 Celavirus
>   9 Ipomovirus
>  14 Macluravirus
>   3 Poacevirus
> 316 Potyvirus
>   2 Roymovirus
>   3 Rymovirus
>   6 Tritimovirus
> ```
>
> </details>

---
