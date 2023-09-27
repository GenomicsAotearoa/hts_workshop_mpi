# Manipulating files in the shell

!!! clock "time"

    * Teaching: 15 minutes
    * Exercises: 15 minutes
    
!!! circle-info "Objectives and Key points"

    #### Objectives
    
    * Search the contents of basic text files for specific strings.
    * Copy and remove directories, create tarballs and compress files.
    
    #### Keypoints
    
    * You can view search file contents using the `grep` command.
    * The commands `cp` and `rm` can be applied to directories, with the correct parameters.
    * You can perform find and replace operations in files using the `sed` command.

---

## Searching files using `grep`

As powerful as `less` can be, at some point it becomes impractical to use it to screen documents as even with the search tools we are retrieving too many search hits, or covering too much content to reasonably interpret it by eye. For such situations, we can use another command line tool to search through documents without opening them and report the matches directory to the terminal. This tool is called `grep`.

`grep` is a command-line utility for searching plain-text files for lines matching a specific set of characters (sometimes called a string) or a particular pattern. We are going to work with one of the fastq files to practice using `grep` and to demonstrate some of the tool features. For these exercises we will searching some fastq files based on theri sequence content. As we are probably all aware, the four nucleotides that appear in DNA are abbreviated `A`, `C`, `T` and `G`. According to the IUPAC (International Union of Pure and Applied Chemistry) code, unknown nucleotides are represented with the letter `N` (a**N**y base). An `N` appearing in a sequencing file represents a position where the sequencing machine was not able to confidently determine the nucleotide in that position.

We'll search for strings inside of our fastq files. Let's first make sure we are in the correct directory:

!!! terminal "code"

    ```bash
    cd /nesi/project/nesi03181/phel/USERNAME/level2/shell_data/
    ```

Suppose we want to see how many reads in our file have really bad segments containing 10 consecutive unknown nucleotides (`N`s).

Let's search for the string `NNNNNNNNNN` in the `SRR098026.fastq` file:

!!! terminal "code"

    ```bash
    grep NNNNNNNNNN SRR098026.fastq
    ```

This command returns a lot of output to the terminal. Every single line in the `SRR098026.fastq` file that contains at least 10 consecutive `N`s is printed to the terminal, regardless of how long or short the file is. We may be interested not only in the actual sequence which contains this string, but in the name (or identifier) of that sequence. We discussed in a previous lesson that the identifier line immediately precedes the nucleotide sequence for each read in a fastq file. We may also want to inspect the quality scores associated with each of these reads. To get all of this information, we will return the line immediately before each match and the two lines immediately after each match (see [this description of the fastq format](../supplementary/fastq_format.md) if you are unsure why we are using these numbers).

We can use the `-B` argument for grep to return a specific number of lines before each match. The `-A` argument returns a specific number of lines after each matching line. Here we want the line *before* and the two lines *after* each matching line, so we add `-B1 -A2` to our grep command:

!!! terminal "code"

    ```bash
    grep -B1 -A2 NNNNNNNNNN SRR098026.fastq
    ```

??? success "Output"

    ```
    @SRR098026.177 HWUSI-EAS1599_1:2:1:1:2025 length=35
    CNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    +SRR098026.177 HWUSI-EAS1599_1:2:1:1:2025 length=35
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ```

!!! question "Exercise"

    1. Search for the sequence `GNATNACCACTTCC` in the `SRR098026.fastq` file. Have your search return all matching lines and the name (or identifier) for each sequence that contains a match.

    ??? circle-check "Solution"
 
        !!! terminal "code"
        
            ```bash
            grep -B1 GNATNACCACTTCC SRR098026.fastq
            ```

        ??? success "Output"

            ```
            @SRR098026.245 HWUSI-EAS1599_1:2:1:2:801 length=35
            GNATNACCACTTCCAGTGCTGANNNNNNNGGGATG
            ```

!!! question "Exercice"

    2. Search for the sequence `AAGTT` in both FASTQ files. Have your search return all matching lines and the name (or identifier) for each sequence that contains a match.

    ??? circle-check "Solution"
 
        !!! terminal "code"

            ```bash
            grep -B1 AAGTT *.fastq
            ```

        ??? success "Output"
        
            ```
            SRR097977.fastq-@SRR097977.11 209DTAAXX_Lenski2_1_7:8:3:247:351 length=36
            SRR097977.fastq:GATTGCTTTAATGAAAAAGTCATATAAGTTGCCATG
            --
            SRR097977.fastq-@SRR097977.67 209DTAAXX_Lenski2_1_7:8:3:544:566 length=36
            SRR097977.fastq:TTGTCCACGCTTTTCTATGTAAAGTTTATTTGCTTT
            --
            SRR097977.fastq-@SRR097977.68 209DTAAXX_Lenski2_1_7:8:3:724:110 length=36
            SRR097977.fastq:TGAAGCCTGCTTTTTTATACTAAGTTTGCATTATAA
            --
            SRR097977.fastq-@SRR097977.80 209DTAAXX_Lenski2_1_7:8:3:258:281 length=36
            SRR097977.fastq:GTGGCGCTGCTGCATAAGTTGGGTTATCAGGTCGTT
            --
            SRR097977.fastq-@SRR097977.92 209DTAAXX_Lenski2_1_7:8:3:353:318 length=36
            SRR097977.fastq:GGCAAAATGGTCCTCCAGCCAGGCCAGAAGCAAGTT
            --
            SRR097977.fastq-@SRR097977.139 209DTAAXX_Lenski2_1_7:8:3:703:655 length=36
            SRR097977.fastq:TTTATTTGTAAAGTTTTGTTGAAATAAGGGTTGTAA
            --
            SRR097977.fastq-@SRR097977.238 209DTAAXX_Lenski2_1_7:8:3:592:919 length=36
            SRR097977.fastq:TTCTTACCATCCTGAAGTTTTTTCATCTTCCCTGAT
            --
            SRR098026.fastq-@SRR098026.158 HWUSI-EAS1599_1:2:1:1:1505 length=35
            SRR098026.fastq:GNNNNNNNNCAAAGTTGATCNNNNNNNNNTGTGCG
            ```

---

## Copying and removing folders of files

In the level 1 training we covered how to copy files and remove empty folders (directories), but it was noted that the copy operation would not work on folders and the `rmdir` command did not work on folders containing files.

If you are in a position where you need to copy a directory of files to a new location, such as when creating backups this can be achieved by adding the `-r` (recursive) flag to the copy (`cp`) command:

!!! terminal "code"

    ```bash
    cp -r shell_data/ shell_data_backup/
    ```

In the level 1 training the way to remove the `shell_data_backup/` folder would be to navigate into the directory, remove all files, then navigate out of the directory and remove it with the `rmdir` command. If you are careful this can be compressed into a single command, by adding the `-r` flag to the remove (`rm`) command.

!!! terminal "code"

    ```bash
    rm -r shell_data_backup/
    ```

This is subject to the usual `bash` warning that contents removed in this way is not recoverable. If you try to perform this operation on a folder with contents which have had their write permissions removed you will receive a confirmation prompt for each file in the directory. This may or may no be a problem depending on the number of files that this applies to.

!!! question "Exercise"
    Use the `rm` help content to find the way to delete a folder of files bypassign the need to approach each file for detetion.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            rm -rf shell_data_backup/
            ```

??? circle-info "Setting file permissions"

    If it's that easy to permanently delete a file, how can be put some checks in place to prevent it from happening accidentally? The answer to this is through file permissions. In computer systems such as NeSI, the ability to read and write files is stored as metadata in each file. This can be viewed and modified using the appropriate commands.
    
    For example, you can run a modified version of the `ls` command to see the permissions of files:

    !!! terminal "code"

        ```bash
        ls -l
        ```

    ??? success "Output"

        ```
        drwxrws---+ 2 dwaite comm00008  4096 Feb 24 15:50 backup
        drwxrws---+ 2 dwaite comm00008  4096 Feb 24 15:50 other_backup
        -rw-rw-r--+ 1 dwaite comm00008 47552 Jun 25  2021 SRR097977.fastq
        -rw-rw-r--+ 1 dwaite comm00008 43332 Jun 25  2021 SRR098026.fastq
        ```

    What we interested in is the first part of the output, the strings of characters which look like `-rw-rw-r--+`. This represents the current permission state of the file. If we ignore the first and last characters, what we are left with are 9 characters, which represent 3 file permissions for 3 different user groups, like so:

    ![Permissions breakdown](../img/level2_11_rwx_figure.svg)

    We're going to concentrate on the three positions that deal with your permissions (as the file owner), which will have the values `rw-`. This shows that you have permission to read (`r`) the file as well as edit the contents (write, `w`). The third position is set to `-`, indicating that you don't have permission to carry out the ability encoded by that space. This is the space where the ability to execute the file (`x`) is set.

    It is possible to modify these permissions, to expand the read/write access to other groups, or restrict the ability to read or write a file using the `chmod` command. However, we generally discourage the use of this within the PHEL environment as our data are a shared resource and restricting the ability of others to work with out data can lead to problems when handing over data related to diagnostic work.

---

## Modifying file contents using `sed`

Sometimes making manual changes to text files is tedious, doubly so when working through command line text editors. Fortunately, there is a command line tool available which can do find-and-replace operations on files or data. The tool `sed` (**s**tream **ed**itor) is used to make changes to the contents of files, either in place or by printing the modified contents to a new file.

Navigate to the `shell_data/` folder, and inspect the contents of the file `SRR097977.small.fq`.

!!! terminal "code"

    ```bash
    cd /nesi/project/nesi03181/phel/USERNAME/level2/shell_data/
    cat SRR097977.small.fq
    ```

This is just the first 20 lines (5 sequences) from the `SRR097977.fastq` file - we are going to work with this today as we have not yet learned the appropriate commands for saving the output of command line operations in to files, so working with the full verison will be an information overload on the terminal.

When running `sed`, there are two pieces of information we need - the file stream to be modified, and an *expression* telling the tool how to make the replacement. This expression is the most complicated part of working with `sed`, but it is a very powerful tool when mastered. For this exercise we are going to just change a few pieces of information in the sequence headers just to get a feel for how the tool operates.

As a worked example, take a look at the following command:

!!! terminal "code"

    ```bash
    sed "s/A/B/" SRR097977.small.fq
    ```

    ??? success "Output"

        ```
        @SRR097977.1 209DTBAXX_Lenski2_1_7:8:3:710:178 length=36
        TBTTCTGCCATAATGAAATTCGCCACTTGTTAGTGT
        +SRR097977.1 209DTBAXX_Lenski2_1_7:8:3:710:178 length=36
        CCCCCCCCCCCCCCC>CCCCC7CCCCCCBCA?5A5<
        @SRR097977.2 209DTBAXX_Lenski2_1_7:8:3:365:371 length=36
        GGTTBCTCTTTTAACCTTGATGTTTCGACGCTGTAT
        +SRR097977.2 209DTBAXX_Lenski2_1_7:8:3:365:371 length=36
        CC:?:CC:?CCCCC??C?:?C-&:C:,?<&*?+7?<
        @SRR097977.3 209DTBAXX_Lenski2_1_7:8:3:663:569 length=36
        TTGTTCGCTTTTGGTBATTAATCCCGGAAATAATAA
        +SRR097977.3 209DTBAXX_Lenski2_1_7:8:3:663:569 length=36
        CCCCCCCCCCCC&9BACCC,C>CCAA&0?4A9&A<6
        @SRR097977.4 209DTBAXX_Lenski2_1_7:8:3:715:205 length=36
        TBTCACTAAAGATCAAATCATTGAAGCAGTTGCAGC
        +SRR097977.4 209DTBAXX_Lenski2_1_7:8:3:715:205 length=36
        CCCCCCC:CCC:CCC:CCC9CC??CCCC?0?*?1--
        @SRR097977.5 209DTBAXX_Lenski2_1_7:8:3:639:209 length=36
        TBTCTATCAAAGCCAGGCAATGGAAGACCTACTCCC
        +SRR097977.5 209DTBAXX_Lenski2_1_7:8:3:639:209 length=36
        CCCCCCCCC?C?CC3C?CC5C?C1C<?CC8BA+AA%
        ```

In that command, there should be two obvious pieces of information - the name of the `sed` tool, and the name of the file we are working with. The remaining part, enclosed by quotation marks, is the expression for what to change. This takes the following form;

!!! terminal "code"

    ```
    s/TEXT_TO_BE_REPLACED/TEXT_TO_REPLACE_WITH/
    ```

The `/` characters are used to separate the `s` flag, and the find/replace values. you can use any character here as long as they are used consistently. For example, the following three commands all create the same output:

!!! terminal "code"

    ```bash
    sed "s/A/B/" SRR097977.small.fq
    sed "s|A|B|" SRR097977.small.fq
    sed "s,A,B," SRR097977.small.fq
    ```

The `/` character is typically used in online examples, but it can be useful to know it can be replaced if you find the notation confusing or are ever in a situation when you need to find and replace a `/` character.

If you run any of these, you will probably spot that even though we are specifying that the `A` character be replaced with `B`, there are plenty of `A`s in each line of the output. By default, `sed` only replaces the first instance of the text it is searching for. We can modify the expression to include all instances by making it the following:

!!! terminal "code"

    ```bash
    sed "s/A/B/g" SRR097977.small.fq
    ```

    ??? success "Output"

        ```
        @SRR097977.1 209DTBBXX_Lenski2_1_7:8:3:710:178 length=36
        TBTTCTGCCBTBBTGBBBTTCGCCBCTTGTTBGTGT
        +SRR097977.1 209DTBBXX_Lenski2_1_7:8:3:710:178 length=36
        CCCCCCCCCCCCCCC>CCCCC7CCCCCCBCB?5B5<
        @SRR097977.2 209DTBBXX_Lenski2_1_7:8:3:365:371 length=36
        GGTTBCTCTTTTBBCCTTGBTGTTTCGBCGCTGTBT
        +SRR097977.2 209DTBBXX_Lenski2_1_7:8:3:365:371 length=36
        CC:?:CC:?CCCCC??C?:?C-&:C:,?<&*?+7?<
        @SRR097977.3 209DTBBXX_Lenski2_1_7:8:3:663:569 length=36
        TTGTTCGCTTTTGGTBBTTBBTCCCGGBBBTBBTBB
        +SRR097977.3 209DTBBXX_Lenski2_1_7:8:3:663:569 length=36
        CCCCCCCCCCCC&9BBCCC,C>CCBB&0?4B9&B<6
        @SRR097977.4 209DTBBXX_Lenski2_1_7:8:3:715:205 length=36
        TBTCBCTBBBGBTCBBBTCBTTGBBGCBGTTGCBGC
        +SRR097977.4 209DTBBXX_Lenski2_1_7:8:3:715:205 length=36
        CCCCCCC:CCC:CCC:CCC9CC??CCCC?0?*?1--
        @SRR097977.5 209DTBBXX_Lenski2_1_7:8:3:639:209 length=36
        TBTCTBTCBBBGCCBGGCBBTGGBBGBCCTBCTCCC
        +SRR097977.5 209DTBBXX_Lenski2_1_7:8:3:639:209 length=36
        CCCCCCCCC?C?CC3C?CC5C?C1C<?CC8BB+BB%
        ```

!!! question "Exercise"

    1. Create a `sed` expression to replace the `SRR097977` text in each line with your name in the `SRR097977.small.fq` file.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            sed "s/SRR097977/Bob/" SRR097977.small.fq
            ```

    2. Create a `sed` expression to replace the `length=36` text with nothing.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            sed "s/ length=36//" SRR097977.small.fq
            ```

---
