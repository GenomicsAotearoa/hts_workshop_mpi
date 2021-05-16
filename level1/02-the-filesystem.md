# Navigating Files and Directories

* Teaching: 30 minutes
* Exercises: 20

#### Objectives

* Use a single command to navigate multiple steps in your directory structure, including moving backwards (one level up).
* Perform operations on files in directories outside your working directory.
* Work with hidden directories and hidden files.
* Interconvert between absolute and relative paths.
* Employ navigational shortcuts to move around your file system.

#### Keypoints

* The `/`, `~`, and `..` characters represent important navigational shortcuts.
* Hidden files and directories start with `.` and can be viewed using `ls -a`.
* Relative paths specify a location starting from the current location, while absolute paths specify a location from the root of the file system."

---

## Contents

1. [Moving around the file system](#moving-around-the-file-system)
1. [Finding hidden directories](#finding-hidden-directories)
1. [Examining the contents of other directories](#examining-the-contents-of-other-directories)
1. [Full vs. relative paths](#full-vs-relative-paths)
1. [Navigational shortcuts](#navigational-shortcuts)

---

## Moving around the file system

We've learned how to use `pwd` to find our current location within our file system. We've also learned how to use `cd` to change locations and `ls` to list the contents of a directory. Now we're going to learn some additional commands for moving around  within our file system.

Use the commands we've learned so far to navigate to the `shell_data/untrimmed_fastq/` directory, if you're not already there.

> <details>
> <summary>Solution</summary>
> 
> ```bash
> $ cd /nesi/nobackup/nesi03181/username/shell_data/untrimmed_fastq/
> ```

What if we want to move back up and out of this directory and to our top level  directory? Can we type `cd shell_data/`? Try it and see what happens.

```bash
$ cd shell_data
-bash: cd: shell_data: No such file or directory
```

Your computer looked for a directory or file called `shell_data/` within the directory you were already in. It didn't know you wanted to look at a directory level above the one you were located in. 

We have a special command to tell the computer to move us back or up one directory level. 

```bash
$ cd ..
```

Now we can use `pwd` to make sure that we are in the directory we intended to navigate to, and `ls` to check that the contents of the directory are correct.

```bash
$ pwd
/nesi/nobackup/nesi03181/username/shell_data
```

```bash
$ ls
sra_metadata  untrimmed_fastq
```

From this output, we can see that `..` did indeed take us back one level in our file system. You can chain these together. For example:

```bash
$ ls ../../
```

Prints the contents of `/nesi/nobackup/nesi03181/`.

---

## Finding hidden directories

Return to the `shell_data` directory, if you left it in the previous exercise. There is a hidden directory within this directory.

> ### Exercise
>
> Explore the options for `ls` to find out how to see hidden directories. List the contents of the directory and identify the name of the text file in that directory.
> 
> Hint: in Unix file systems, hidden files and directories are prefixed with the `.` character.
> 
> <details>
> <summary>Solution</summary>
> 
> First use the `man` command to look at the options for `ls`. 
> 
> ```bash
> $ man ls
> ```
> 
> The `-a` option is short for `all` and says that it causes `ls` to "not ignore entries starting with ." This is the option we want. 
> 
> ```bash
> $ ls -a
> .  ..  .hidden	sra_metadata  untrimmed_fastq
> ```
> 
> The name of the hidden directory is `.hidden`. We can navigate to that directory using `cd`.
> </details>

Once you have found the hidden directory, list the contents of the directory using `ls`. 

```bash
 $ ls
youfoundit.txt
```

---

## Examining the contents of other directories

By default, the `ls` commands lists the contents of the working directory (i.e. the directory you are in). You can always find the directory you are in using the `pwd` command. However, you can also give `ls` the names of other directories to view. Navigate to your home directory if you are not already there.

```bash
$ cd ~
```

Then enter the command:

```bash
$ ls /nesi/nobackup/nesi03181/username/shell_data/
sra_metadata  untrimmed_fastq
```

This will list the contents of the `shell_data/` directory without you needing to navigate there. You might have already spotted that the `cd` command works in a similar way.

```bash
$ cd /nesi/nobackup/nesi03181/username/shell_data/untrimmed_fastq/
```

---

## Full vs. relative paths

The `cd` command takes an argument which is a directory name. Directories can be specified using either a *relative* path or a full *absolute* path. The directories on the computer are arranged into a hierarchy. The full path tells you where a directory is in that hierarchy. Dependending on where you currently are in the file system when you enter the `pwd` command you will see something like

```bash
$ pwd
/nesi/nobackup/nesi03181/username/
```

or

```bash
$ pwd
/home/username/
```

This is the full name of your current directory. Assuming you saw the first output, this tells you that you are in a directory called `username/`, which sits inside a directory called `nesi03181/` which in turn sits inside a directory `nobackup/`. At the very top of the hierarchy is a directory called `/` which is usually referred to as the **root directory**.

Return to your home directory, then navigate to the `.hidden/` folder using the following command:

```bash
$ cd /nesi/nobackup/nesi03181/username/shell_data/.hidden/
```

Now return to your home directory again, and navigate back to the `.hidden/` folder using the following commands:

```bash
$ cd /
$ cd nesi/
$ cd nobackup/
$ cd nesi03181/
$ cd username/
$ cd shell_data/
$ cd .hidden/
```

These two commands have the same effect, they both take us to the `.hidden` directory. However, the first uses the absolute path, giving the full address from the home directory. The second uses a series of relative paths, with the directory specified in each command contingent on the current working directory.

A relative path is like getting directions from someone on the street. They tell you to "go right at the stop sign, and then turn left on Main Street". That works great if you're standing there together, but not so well if you're trying to tell someone how to get there from another country. A full path is like GPS coordinates. It tells you exactly where something is no matter where you are right now.

You can usually use either a full path or a relative path depending on what is most convenient. If we are in the home directory, it is more convenient to enter the full path. If we are in the working directory, it is more convenient to enter the relative path since it involves less typing.

Over time, it will become easier for you to keep a mental note of the structure of the directories that you are using and how to quickly navigate amongst them.

A full path always starts with a `/` (the root directory). A relative path does not. This is a helpful difference to remember so that you always know which type of path you are working with.

> ### Exercise
>
> Using the filesystem diagram below, if `pwd` displays `/Users/thing`, what will `ls ../backup` display?
> 1.  `../backup: No such file or directory`
> 2.  `2012-12-01 2013-01-08 2013-01-27`
> 3.  `2012-12-01/ 2013-01-08/ 2013-01-27/`
> 4.  `original pnas_final pnas_sub`
> 
> ![](../img/01_filesystem_challenge.svg)
> 
> <details>
> <summary>Solution</summary>
> 
> 1. No: there *is* a directory `backup/` in `/Users`.
> 2. No: this is the content of `Users/thing/backup`, but with `..` we asked for one level further up.
> 3. No: see previous explanation. Also, we did not specify `-F` to display `/` at the end of the directory names.
> 4. Yes: `../backup` refers to `/Users/backup`.
> 
> </details>

---

## Navigational shortcuts

The root directory is the highest level directory in your file system and contains files that are important for your computer to perform its daily work. While you will be using the root (`/`) at the beginning of your absolute paths, it is important that you avoid working with data in these higher-level directories, as your commands can permanently alter files that the operating system needs to function.

In many cases, including when working on NeSI, trying to run commands in root directories will require special permissions which are not available to you as a regualar user.

Dealing with the home directory is very common. The tilde character, `~`, is a shortcut for your home directory. On a Linux operating system the root directory is **two** levels above our home directory, so `cd` or `cd ~` will take you to `/home/username/` and `cd /` will take you to `/`.

~~The commands `cd`, and `cd ~` are very useful for quickly navigating back to your home directory. We will be using the `~` character in later lessons to specify our home directory.~~
