# 2.4 - Polishing of Oxford Nanopore assemblies

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

Today we are going to use a long read polising tool called `racon` ([Vaser *et al.*, 2017](https://doi.org/10.1101/gr.214270.116)).

??? Question "What are our other options?"

    1. medaka ([GitHub](https://github.com/nanoporetech/medaka))
    1. FMLRC ([Wang *et al.*, 2018](https://doi.org/10.1186/s12859-018-2051-3))
    1. LoRDEC ([Salmela *et al.*, 2014](https://doi.org/10.1093/bioinformatics/btu538))

    A comprehensive comparison of polishing tools was published in 2019 ([Fu *et al.*, 2019](https://doi.org/10.1186/s13059-018-1605-z)), which is still a useful reference for getting started if you are using short reads to polish a long read assembly. When working with exclusively long read data, `racon` is probably the best (most robust and universally applied) tool for a first attempt at error correction.

Whether or not our assembly will benefit from polishing is hard to predict. When working through this process it is a good idea to make copies of your data as you perform different correction procedures (or combinations of procedures) and to evaluate each outcome.

In theory, good assemblies can be obtained by applying the workflow;

1. Assemble
1. Polish with `racon`
   1. Repeat up to four times
1. Polish with `medaka`

But this is not always the case, and you often need to compare each iteration of the assembly to a reference to see how different regions are affects.

---

## Performing the first round of polishing

Strictly speaking, `racon` is designed for polishing assemblies which have been obtained from tools that do not perform extensive error correction themselves. However, we are going to apply it to the `Flye` assembly as we have limited time in this workshop and it is unlikely to negatively affect our assembly.

Before running `racon` we must produce a mapping file of the quality filtered sequences against the assembly. We can do this with `minimap2`.

!!! note ""

    We will work with `minimap2` more in the third session so will not explain it's parameters and workflow today.

!!! terminal "code"

    ```bash
    module load minimap2/2.24-GCC-11.3.0

    minimap2 -t 4 -ax map-ont reference/flye.fna reads/Mbovis_87900.nanopore.fq.gz > Mb1.sam
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
    [M::main] CMD: minimap2 -t 4 -ax map-ont reference/flye.fna reads/Mbovis_87900.nanopore.fq.gz
    [M::main] Real time: 6.517 sec; CPU: 12.478 sec; Peak RSS: 0.187 GB
    ```

We can then use this mapping file as the input for `racon`:

!!! terminal "code"

    ```bash
    module load Racon/1.5.0-GCC-11.3.0

    racon -t 4 reads/Mbovis_87900.nanopore.fq.gz Mb1.sam reference/flye.fna > flye.racon_1.fna
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

---

## Performing additional rounds of polishing

It is possible to perform the `racon` process iteratively, remapping reads to the output and then running the polishing cycle again. There is some data ([link here](https://nanoporetech.github.io/medaka/draft_origin.html#discussion)) which suggests that up to four rounds of `racon` polishing, in conjunction with `medaka`, produces better quality output than running a single polishing step.

However there are costs associated with this approach both in terms of time invested and over-zealous correction to repeat regions. Whether or not improvement with multiple rounds will be seen in your data is unclear, and ultimately it is your decision whether or not to perform this approach so although this can work, it is a judgement call as to whether or not it is necessary.

!!! question "Exercise"

    Run a second round of `racon` polishing, using the output of your first iteration as the input for the second polishing round.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            minimap2 -t 4 -ax map-ont flye.racon_1.fna reads/Mbovis_87900.nanopore.fq.gz > Mb2.sam

            racon -t 4 reads/Mbovis_87900.nanopore.fq.gz Mb2.sam flye.racon_1.fna > flye.racon_2.fna
            ```

---

## Post-polishing follow up

One important piece on information to note with the polishing process is that `racon` **_does not_** change the names of contigs during polishing. This is helpful, as it allows us to easily compare contigs between different polishing steps but it also means that you have to be careful when importing the data into `Geneious` as it might become hard to track which step of the analysis your contig comes from.

As an easy solution to this is to rename your contigs using `seqmagick` to append some versioning information to each sequence name:

!!! terminal "code"

    ```bash
    seqtk/1.4-GCC-11.3.0

    seqtk rename reference/flye.fna "BASE_" > flye.base.rename.fna
    seqtk rename flye.racon_1.fna "RACON1_" > flye.racon_1.rename.fna
    seqtk rename flye.racon_1.fna "RACON2_" > flye.racon_2.rename.fna
    ```

---
