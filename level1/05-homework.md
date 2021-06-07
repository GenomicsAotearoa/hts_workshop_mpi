# Level 1 revision

This revision material consists of five short exercises covering the main points of the Level 1 training. When you have completed all exercises, please contact one of the training representatives to have them check your work.

If you are struggling with any exercises, remember the following:

1. Every exercise here was covered in the training material. Refer to this for hints.
2. It is not cheating to use Google to find help on how various programmes work.
3. You are welcome to contact the trainers for advice on a particular exercise, but please attempt the first two options before resorting to this.

---

## Exercise 1 - Navigating the file system

In the workshop we were mostly working in a diretory called `untrimmed_fastq/`. However, there was a second folder in the `shell_data/` directory. Using your knowledge of the `cd` and `ls` commands, navigate to this other folder and find out what it contains. Give a brief desription of the folder contents to the trainer, using whatever commands you feel best suited to viewing the file contents.

>**NOTE:** If you do not remember the location of the data, all work for this workshop is located in the `/nesi/project/nesi03181/phel/` directory.

---

## Exercise 2 - Relative pathing

Once you have found the contents of the other folder, use your knowledge of relative pathing to write the shortest `cd` command to move from the folder found in **Exercise 1** into the `untrimmed_fastq/`.

---

## Exercise 3 - Searching a file with `grep`

In the folder you discovered in **Exercise 1** the file contains information regarding the two fastq files we worked with in the workshop as well as a number of other samples. Use the `grep` tool to find the entry for sample *SAMN00205533*, reporting not just the line contents but also the line number for this entry.

>**NOTE:** The `man` command can be used to find additional options and filters for the grep command. Alternatively, you can run the command as
>```bash
>$ grep --help
>```
>To see the options for `grep`.

---

## Exercise 4 - Writing a `for` loop

Navigate back into the `untrimmed_fastq/` folder. Imagine a situation where we wanted to capture the first and last sequence entry for a set of fastq files. Use your knowledge of the `head` and `tail` commands to write a `for` loop which will take the first and last sequence from a set of fastq files and redirect them into a new file.

Some details to consider:

1. Each sequence in a fastq file spans 4 lines.
2. You will need to run the `head` and `tail` commands as separate commands in the `for` loop. Make sure you use the correct redirection command to prevent the `tail` command overwriting the results of the `head` command.
3. Each iteration of the loop should write output to a different file. It is your choice what to call the output file, but if you choose to use the `basename` command to generate an output name, the command that we used in the workshop was `o=$(basename ${i} .fastq)`, where `i` was the variable containing the input file name.

---

## Exercise 5 - Creating backups and changing file permissions

We now want to create a backup of the outputs from **Exercise 4**. Create a new *hidden* directory with a name of your choice, then copy the files created in **Exercise 4** to this location. Once the files have been copied, restrict access to the versions inside the hidden directory by removing **write** and **execute** permissions from the files.

---

## Reporting your work

When you have completed all 5 exercises create a log of the remaining exercises, perform the following commands:

```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/shell_data/

$ history > module_1.2.homework.txt
```

Once this is complete, email the trainers telling them your answers to **Exercise 1** and to let them know that the homework log file is ready to be checked.
