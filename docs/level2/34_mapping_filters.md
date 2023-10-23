# 3.4 - Filtering and compressing sam files

!!! clock "time"

    * Teaching: 20 minutes
    * Exercises: 20 minutes

!!! circle-info "Objectives and Key points"

    #### Objectives
    
    * Use `samtools` to sort and compress a raw `sam` file into the `bam` format.
    * Use `samtools` to filter a `bam` file into either the successfully mapped, or unmapped reads.
    * Use `samtools` to recover reads in `fastq` format from a `bam` file.
    
    #### Keypoints
    
    * Understand the reasons for sorting and compressing files in the `sam` and `bam` formats.
    * Understand the situations in which you may wish to filtering a `sam`/`bam` file and what the downstream applications of the output would be.

---

## Why we need to compress and filter `sam` files

Now that we have created some basic mapping data in the `sam` files , let's take a look at the size of these files. Navigate to the `/nesi/project/nesi03181/phel/USERNAME/level2/mapping/` directory and run the following command:

!!! terminal "code"

    ```bash
    ls -sh *.sam
    ```

    ??? success "Output"

        ```
        1.3G Mbovis_87900.16S_rRNA.bowtie2.sam
        1.5G Mbovis_87900.genome.bowtie2.sam
        52M Mbovis_87900.genome.nanopore.sam
        ```

These are quite large files - especially in the case of the 16S rRNA mapping file where the mapping success rate was only around 0.4% of the Illumina sequences. Why is this file so far? There are two reasons.

!!! info "When mapping is performed, *all* of the reads in the input file are recorded"

    This is a deliberate design feature of the mapping tools as depending on our context we may prefer to have the mapped or unmapped reads. For example, if we are trying to perform a reference-based genome assembly or determine the sequencing coverage of the genome (or a particular region) having the mapped reads is critical.

    However, if we are trying to subtract the host reads from a set of reads, such as when we are trying to find a pathogen, the reads which do not map to the host genome are actually the ones we are interested in examining in further detail.

Although there are only about 8,700 reads actually mapped to the reference sequence in the 16S rRNA file (we'll see how to calculate that number in the next session), the file still contains all 2 million reads from the Illumina library.

---

## Sorting and compressing `sam` files

The first step of reducing the file size is to efficienctly compress the contents of the sam file. Fortunately there is a built-in solution for this - the sam file specification also has a binary-encoded equivalent which records the exact same information, in a much more efficieny format. This is the Binary Alignment/Map (`bam`) format.

When performing this compression from sam to bam we also use this opportunity to sort the mapped reads in terms of their starting position in the reference sequence. This sorting is important as it increases the speed of many of the operations we need to perform using a `sam` file, particular when producing coverage statistics. Because the sorting and compression can be performed in a single command line operation we tend to do these things together once and never worry about it again.

!!! info "We should always sort our sequences!"

    In the situations where you need your mapped reads sorted, operations will fail if the reads are not sorted. In situations where you do not need them sorted, operations will succeed.

    It is generally just easier to sort as soon as mapping completes then never worry about it again.

The main tool used for handling sam and bam files is called `samtools`. Load the module and execute the file below. While it is running, refresh yourself on what the `|` operator in the command is doing.

!!! terminal "code"

    ```bash
    module purge
    module load SAMtools/1.16.1-GCC-11.3.0

    samtools view -bS Mbovis_87900.16S_rRNA.bowtie2.sam | samtools sort -o Mbovis_87900.16S_rRNA.bowtie2.bam
    ```

??? question "What is the pipe for..."

    This is redirection, taking the output of the first `samtools` command and passing it as input into the second command.

    The reason we need to redirect the data from one `samtools` command to another is due to the behaviour of the `samtools sort` subcommand. If you examine the manual for this command, you will see that while the *output* of the command can be written in sam or bam format the *input* must be bam.

    We therefore need to use the `samtools view` subcommand to first convert the sam file into bam format. Rather than write the results to NeSI's hard drive then perform `samtools sort` as a second command we can redirect between the commands to save on hard drive space and speed up the operation.

Compare the file sizes between the `sam` and `bam` files:

!!! terminal "code"

    ```bash
    ls -sh Mbovis_87900.16S_rRNA.bowtie2.sam Mbovis_87900.16S_rRNA.bowtie2.bam
    ```

    ??? success "Output"

        ```
        1.3G Mbovis_87900.16S_rRNA.bowtie2.sam
        470M Mbovis_87900.16S_rRNA.bowtie2.bam
        ```

!!! warning "Do we need both files..."

    There is no point in retaining the original sam file at this point, as the information it contains is more efficiently encoded within the bam file.

    If working with real data, this is the point you should delete your sam file. In this training exercise, it does not matter whether you delete the sam file or not.

By compressing the data of the sam file into the bam file we have already reduced its size to about one third of the original. Not only will this save us space on NeSI, it also makes downloading the data much quicker, and also makes subsequent analyses of the file faster as it requires less time to read the smaller bam file.

!!! question "Exercise"

    Sort and compress any other `sam` files you have produced during the previous mapping tutorials.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            for i in bowtie2 nanopore;
            do
                samtools view -bS Mbovis_87900.genome.${i}.sam | samtools sort -o Mbovis_87900.genome.${i}.bam
            done
            ```

---

## Splitting `bam` files to separate mapped and unmapped reads

Most of the time when we are mapping against a reference sequence, we are interested in the sequences which successfully mapped to the target. In these cases, having a bam file which contains all of the unmapped reads is not particularly useful so we will now practice filtering bam files according to the mapping state of the reads in the file.

This is simple to do in terms of the command which needs to be run, but it the meaning of the command can be a bit confusing. Within the sam and bam file format specification is a numeric flag which may be assign to an unmapped read, denoting its status as such. Because it is a *positive* marker of unmapped state, if we want to filter for mapped reads we need to filter for sequences which do not have the marker.

When running the `samtools view` command there is a pair of optional paramters which allow us to filter reads by their mapping flags. We use the `-f` to *keep* reads with the flag(s) we specify, or `-F` to *reject* reads with the flag(s) we specify. As the unmapped flag is only applied to reads which fail to map to the reference we use the following logic:

* `-f 4` = Include reads with the unmapped flag = Keep unmapped reads
* `-F 4` = Reject reads with the unmapped flag = Keep mapped reads

!!! info "The value `4` is the numeric flag for ummapped reads"

We will run `samtools` to filter the `Mbovis_87900.16S_rRNA.bowtie2.bam` file to only keep mapped reads:

!!! terminal "code"

    ```bash
    samtools view -h -F 4 -b Mbovis_87900.16S_rRNA.bowtie2.bam > Mbovis_87900.16S_rRNA.bowtie2.mapped.bam
    ```

!!! warning ""

    In this command, the difference between keeping mapped or unmapped reads is the case of the `-f` flag. It is really easy to get confused about these values so be very careful when applying this command.

In this case we are expecting only a handful of reads to have actually mapped to the reference sequence so a good acid test for whether we used the right commnad will be to check the size of the output file. If it is roughly the same size as the input then we have probably kept the unmapped reads. If it is a lot smaller, we have most likely kept only the mapped reads.

!!! terminal "code"

    ```bash
    ls -sh Mbovis_87900.16S_rRNA.bowtie2.bam Mbovis_87900.16S_rRNA.bowtie2.mapped.bam
    ```

    ??? success "Output"

        ```
        470M Mbovis_87900.16S_rRNA.bowtie2.bam
        1.8M Mbovis_87900.16S_rRNA.bowtie2.mapped.bam
        ```

!!! question "Exercise"

    Write a loop to filter the genome mapping files, creating an output file for both the mapped and unmapped reads for each case.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            for i in bowtie2 nanopore;
            do
                samtools view -h -F 4 -b Mbovis_87900.genome.${i}.bam > Mbovis_87900.genome.${i}.mapped.bam
                samtools view -h -f 4 -b Mbovis_87900.genome.${i}.bam > Mbovis_87900.genome.${i}.unmapped.bam
            done
            ````

---
