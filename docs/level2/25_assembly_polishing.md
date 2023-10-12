# 2.5 - Polishing of Oxford Nanopore assemblies

!!! clock "time"

    * Teaching: 15 minutes
    * Exercises: 15 minutes
    
!!! circle-info "Objectives and Key points"

    #### Objectives

    * Use `racon` to polish a Nanopore assembly

    #### Keypoints

    * Rounds of refinement are sometimes required to reduce error and can improve the initial assembly when working with long rad data.
    * Sometimes the above sentence does not apply and polishing does not really improve the assembly. The only way to know is trying!

---

## Preparing to polish our assembly

Once we have our draft assembly, we want to revisit it and attempt to improve regions by aligning the original reads against the assembled contigs and trying to improve regions which might have been modelled incorrectly by ther assembler. This is an important process as long read assemblers have to make assumptions about the data they are working with, and do not always make the 'correct' call.

The process of using our raw data to re-call areas of the assembly is called polishing, and there are many tools available for performing this task. Similar to the case with assembly, there are specific tools for using short read or long read data for polishing an assembly.

Today we are going to use a long read polishing tool called `racon` ([Vaser *et al.*, 2017](https://doi.org/10.1101/gr.214270.116)).

??? Question "What are our other options?"

    1. medaka ([GitHub](https://github.com/nanoporetech/medaka))
    1. FMLRC ([Wang *et al.*, 2018](https://doi.org/10.1186/s12859-018-2051-3))
    1. LoRDEC ([Salmela *et al.*, 2014](https://doi.org/10.1093/bioinformatics/btu538))

    A comprehensive comparison of polishing tools was published in 2019 ([Fu *et al.*, 2019](https://doi.org/10.1186/s13059-018-1605-z)), which is still a useful reference for getting started if you are using short reads to polish a long read assembly. When working with exclusively long read data, `racon` is probably the best (most robust and universally applied) tool for a first attempt at error correction.

!!! question "Do we always need to polish an assembly?"

    Whether or not our assembly will benefit from polishing is hard to predict. When working through this process it is a good idea to make copies of your data as you perform different correction procedures (or combinations of procedures) and to evaluate each outcome. You often need to compare each iteration of the assembly to a reference to see how different regions are affects.

---

## Performing the first round of polishing

Strictly speaking, `racon` is designed for polishing assemblies which have been obtained from tools that do not perform extensive error correction themselves. In practice, it rarely has a negative impact on assembly quality so while it doesn't hurt to apply it to as assembly from a tool like `Flye`, it is not always worthwhile.

We will start working with one of several draft genome assemblies, of varying quality (in terms of mismatches to the reference). We will start with the `draft_moderate.fna` genome, which has a 2% rate of mismatch with the reference genome.

Before running `racon` we must produce a mapping file of the quality filtered sequences against the assembly. We can do this with `minimap2`. We will work with `minimap2` more in the mapping exercises in the next session, so will not explain its parameters and workflow today.

!!! terminal "code"

    ```bash
    module load minimap2/2.24-GCC-11.3.0

    minimap2 -t 4 -ax map-ont draft_genomes/draft_moderate.fna reads/Mbovis_87900.nanopore.fq.gz > draft_moderate.sam
    ```

??? success "Output"

    ```
    [M::mm_idx_gen::0.023*0.62] collected minimizers
    [M::mm_idx_gen::0.031*0.79] sorted minimizers
    [M::main::0.031*0.79] loaded/built the index for 1 target sequence(s)
    [M::mm_mapopt_update::0.032*0.80] mid_occ = 10
    [M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 1
    [M::mm_idx_stat::0.033*0.80] distinct minimizers: 47527 (99.61% are singletons); average occurrences: 1.004; average spacing: 5.316; total length: 253631
    [M::worker_pipeline::6.513*1.92] mapped 1766 sequences
    [M::main] Version: 2.24-r1122
    [M::main] CMD: minimap2 -t 4 -ax map-ont draft_genomes/draft_moderate.fna reads/Mbovis_87900.nanopore.fq.gz
    [M::main] Real time: 6.517 sec; CPU: 12.478 sec; Peak RSS: 0.187 GB
    ```

We can then use this mapping file as the input for `racon`:

!!! terminal "code"

    ```bash
    module load Racon/1.5.0-GCC-11.3.0

    racon -t 4 reads/Mbovis_87900.nanopore.fq.gz draft_moderate.sam draft_genomes/draft_moderate.fna > draft_moderate.racon_1.fna
    ```

??? success "Output"

    ```
    [racon::Polisher::initialize] loaded target sequences 0.019926 s
    [racon::Polisher::initialize] loaded sequences 0.415303 s
    [racon::Polisher::initialize] loaded overlaps 0.360286 s
    [racon::Polisher::initialize] aligning overlaps [====================] 0.159391 s
    [racon::Polisher::initialize] transformed data into windows 0.027524 s
    [racon::Polisher::polish] generating consensus [====================] 21.940342 s
    [racon::Polisher::] total = 22.926985 s
    ```

Before we assess the results of this, we will run `racon` over a few different genomes, so that when we assess the final qualities, we can generate a single report for all genomes.

!!! question "Exercise"

    Run `racon` polishing on either the `draft_mild.fna` or `draft_severe.fna` genome (or both!).

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            for assembly in draft_mild draft_severe;
            do
                minimap2 -t 4 -ax map-ont draft_genomes/${assembly}.fna reads/Mbovis_87900.nanopore.fq.gz > ${assembly}.sam
                racon -t 4 reads/Mbovis_87900.nanopore.fq.gz ${assembly}.sam draft_genomes/${assembly}.fna > ${assembly}.racon_1.fna
            done
            ```

---

## Performing additional rounds of polishing

It is possible to perform the `racon` process iteratively, remapping reads to the output and then running the polishing cycle again. There is some data ([link here](https://nanoporetech.github.io/medaka/draft_origin.html#discussion)) which suggests that up to four rounds of `racon` polishing, in conjunction with `medaka`, produces better quality output than running a single polishing step.

However, there are costs associated with this approach both in terms of time invested and over-zealous correction to repeat regions. Whether or not improvement with multiple rounds will be seen in your data is unclear, and ultimately it is your decision whether or not to perform this approach so although this can work, it is a judgement call as to whether or not it is necessary.

Running `racon` the second time is pretty much the same as the first time, except that instead of mapping to the original draft genome, we now map our reads to the output of the first `racon` run and use that alignment for correction.

!!! terminal "code"

    ```bash
    # Map the reads, overwriting our previous sam file
    minimap2 -t 4 -ax map-ont draft_moderate.racon_1.fna reads/Mbovis_87900.nanopore.fq.gz > draft_moderate.sam

    # Perform correction, creating a new output file for the polished genome
    racon -t 4 reads/Mbovis_87900.nanopore.fq.gz draft_moderate.sam draft_moderate.racon_1.fna > draft_moderate.racon_2.fna
    ```

!!! question "Exercise"

    Complete the second rounds of `racon` polishing for the `draft_mild` and/or `draft_severe` genomes.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            for assembly in draft_mild draft_severe;
            do
                minimap2 -t 4 -ax map-ont ${assembly}.racon_1.fna reads/Mbovis_87900.nanopore.fq.gz > ${assembly}.sam
                racon -t 4 reads/Mbovis_87900.nanopore.fq.gz ${assembly}.sam ${assembly}.racon_1.fna > ${assembly}.racon_2.fna
            done
            ```

??? info "Post-polishing follow up"

    One important piece on information to note with the polishing process is that `racon` **_does not_** change the names of contigs during polishing. This is helpful, as it allows us to easily compare contigs between different polishing steps but it also means that you have to be careful when importing the data into `Geneious` as it might become hard to track which step of the analysis your contig comes from.

    As an easy solution to this is to rename your contigs using a tool like `seqtk` to add some versioning information to each sequence name:

    !!! terminal "code"

        ```bash
        seqtk/1.4-GCC-11.3.0

        seqtk rename draft_genomes/draft_moderate.fna "BASE_" > draft_moderate.rename.fna
        seqtk rename draft_moderate.racon_1.fna "RACON1_" > draft_moderate.racon_1.rename.fna
        seqtk rename draft_moderate.racon_2.fna "RACON2_" > draft_moderate.racon_2.rename.fna
        ```

---

## Assessing the results with `QUAST`

As a quick confirmation of how successful the cleaning step was, we will use `QUAST` to compare our raw and polished genomes to the reference genome.

When running this command, modify it to include the genomes your polished as part of the exercises above.

!!! terminal "code"

    ```bash
    module load QUAST/5.2.0-gimkl-2022a

    quast.py -r reference/Mbovis_87900.genome.fna --gene-finding -o quast/ \
        draft_genomes/draft_mild.fna \
        draft_mild.racon_1.fna \
        draft_mild.racon_2.fna # ...plus your genomes
    ```

??? success "Output"

    ```
    Version: 5.2.0

    System information:
      OS: Linux-3.10.0-693.2.2.el7.x86_64-x86_64-with-glibc2.17 (linux_64)
      Python version: 3.10.5
      CPUs number: 2

    Started: 2023-09-28 15:02:30

    # Text omitted...

    Finished: 2023-09-28 15:02:40
    Elapsed time: 0:00:10.422580
    NOTICEs: 4; WARNINGs: 1; non-fatal ERRORs: 0

    Thank you for using QUAST!
    ```

Compare the results of the before and after of polishing. Did this process improve the quality of your genome(s)?

---
