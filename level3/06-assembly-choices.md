# Selecting assembly parameters

* Teaching: 30 minutes

#### Objectives

* Understand the different parameters that can be tuned when performing a *de novo* genome assembly.
* Know the differences between the different genome assembly tools.
* Understand the differences between assembling a genome, metagenome, and (meta)transcriptome.

#### Keypoints

* It is hard to know the best parameters for optimising a particular genomic assembly.
* There are many genome assembly tools available, but some perform better than others in terms of performance or final assembly quality.
* Your choice of assembly tool and run time parameters will vary depending on the type of assembly you are trying to perform.

---

## Contents

1. [Key considerations for optimising an assembly](#key-considerations-for-optimising-an-assembly)
   1. [Assembly modes](#assembly-modes)
   1. [Stages of assembly](#stages-of-assembly)
   1. [Input data and reference files](#input-data-and-reference-files)
   1. [Tuning the assembly](#tuning-the-assembly)
   1. [Checkpointing](#checkpointing)
   1. [Resource limitations](#resource-limitations)

---

## Key considerations for optimising an assembly

If you look through the manual of `SPAdes` you will see a **_lot_** of parameters which can be adjusted for try to optimise the final assembly. There are also a number of parameters to adjust the performance of `SPAdes` itself, to either speed up the job or to prevent it from overloading the server. Here we will discuss some of the key parameter choices which should be considered for each assembly.

```bash
$ module load SPAdes/3.15.2-gimkl-2020a

$ spades.py -h
```

### Assembly modes

In its original inception, `SPAdes` was developed as an assembly tool to address some of the problems inherent in assembly single-cell genome data, as contemporary asssemblers performed poorly on these data. As part of the testing, the authors also noted that `SPAdes` performed extremely well of regular enome sequencing data.

Since it's original release ([Bankevich *et al.*, 2012](https://dx.doi.org/10.1089%2Fcmb.2012.0021)), `SPAdes` has been constantly improved upon, adding new assembly models for a range of different data types as the needs of the community have shifted. When executing `SPAdes` the first important piece of information we must provide the assembler is which of the optimised sub-routines we wish to perform. The current general-purpose assembly options that `SPAdes` provides are as follows:

|Mode|Parameter|Input library|
|:---|:---|:---|
|Isolate sequencing|`--isolate`|Standard pure-culture sequencing library.|
|Single-cell|`--sc`|Single cell sequencing data, typically amplified through MDA (multiple displacement amplification).|
|Metagenome|`--meta`|Assemble a library produced from a mixture of organisms, which have all contributed DNA to the library.|
|Transcriptome|`--rna`|Assemble RNA transcripts/isoforms from an RNA-Seq library.|
|IonTorrent|`--iontorrent`|Assembly of IonTorrent short read data, which are noisier than Illumina libraries.|

`SPAdes` also provides several targetted-assembly modes, which are focused on performing a rapid assembly of particular genetic elements at the expense of the boarder genomic context.

|Mode|Parameter|Input library|
|:---|:---|:---|
|Biosynthetic cluster assembly|`--bio`|Perform an assembly of specific biosynthetic gene clusters, determined by an HMM profile of BGCs of interest.|
|Corona virus|`--corona`|A simplfiied version of the `--rnaviral` mode, specifically targetting a set of SARS-CoV-2 proteins collected from the Pfam database and [Phan *et al.*, 2018](https://doi.org/10.1093/ve/vey035).|
|Plasmid recovery (genome)|`--plasmid`|Targetted identification and assembly of plasmid sequences.|
|Plasmid recovery (metagenome)|`--metaplasmid`|Targetted identification and assembly of plasmid sequences from a metagenomic library.|
|Viral recovery (metagenome)|`--metaviral`|Targetted identification and assembly of viral sequences from a metagenomic library.|
|Viral recovery (transcriptome)|`--rnaviral`|Targetted identification and assembly of viral sequences from a transcriptome or metatranscriptome.|

In practice, we will only use a small number of these assembly modes in routine work. Even, it is important to know which is the correct method for your data as there can be significant differences in the assembly output if the wrong mode is selected, particular with the targetted assembly methods.

### Stages of assembly

There are a number of steps performed by `SPAdes` during assembly, but broadly speaking we are only concerned with two phases of the assembly - **error correction** and **assembly**. When a new run of `SPAdes` begins, the first stage of the workflow is to analyse the input data and identify sequences which may contain obvious sequencing error. If possible, these are corrected to the sequence which they 'should' contain. Performing *a priori* correction of the sequence data streamlines the assembly process, as there are fewer erronous paths in the assembly graph for `SPAdes` to analyse.

However, the error correction process can be quite slow and memory intensive. If desired, it is possible to skip this step using the `--only-assembler` flag when running assembly. Be warned, however, that a large degree of `SPAdes` assembly quality can be attributed to the error correction phase and while disabling it will speed up analysis, the results will usually be inferior to that of the complete pipeline.

By contrast, we can also perform the error correction by itself, producing a new set of 'improved' sequences without moving into the assembly. This can be achieved by invoking the `--only-error-correction` flag when executing `SPAdes`.

### Input data and reference files

`SPAdes` is designed as an assembler of Illumina sequence data, and as such the parameter for defining input files reflect this. Input sequences can be in either interleaved fastq format or the standard forward/reverse/single format we produced in previous exercises. The syntax for passing these files into `SPAdes` is either:

```bash
# Forward/reverse/single input files
$ spades.py --isolate -1 sample_1.fq -2 sample_2.fq -s sample_s.fq ...

# Interleaved input file
$ spades.py --isolate --12 sample_12.fq -s sample_s.fq ...
```

`SPAdes` can also accept additional sequence files to resolve spatial information in the assembly, either through [mate-pair sequencing](https://www.illumina.com/science/technology/next-generation-sequencing/mate-pair-sequencing.html), reference contigs, or long read sequencing technologies.

```
# Run with a mate-pair library to improve assembly
$ spades.py --isolate -1 sample_1.fq -2 sample_2.fq -s sample_s.fq --mp-1 MP_library_1.fq --mp-2 MP_library_2.fq ...

# Run with Nanopore sequences to resolve repeat regions
$ spades.py --isolate -1 sample_1.fq -2 sample_2.fq -s sample_s.fq --nanopore sample.minion.fq ...
```

When using long reads to aid assembly, `SPAdes` can process Nanopore, PacBio, or Sanger sequence data or trusted reference contigs from a previous assembly.

### Tuning the assembly

Modern short read assemblers most commonly use an assembly strategy based on the analysis of short *k*-mer sequence fragments and the construction of a mathematical construct known as a de Bruijn graph. Compared with the more intuitive and traditional overlap layout consensus method, whereby sequences are aligned against each other, the de Bruijn graph approach scales much more efficiently with increasing input data.

The challenge with this approach is that we are subdiving our sequences down into shorter fragments of *k* nucleotides in length and determining the optimal value for *k* is a difficult process. With values of *k* which are too short we are unable to resolve repeat regions (i.e. if the *k*-mer size is shorter than the repeat region, we cannot bridge it with our *k*-mers), yet if *k* is too large we lose sensitivity in joining sequences together.

This is resolved in all modern assemblers by assembling the data with a range of *k*-mer sizes and then aggregating the results. When assembling with `SPAdes` the *k*-mer size(s) to assemble can be provided explicitly to the assembler or we can allow the tool to automatically increase the *k*-mer size from a short starting value until the tool determines that additional effort is no longer justified.

```bash
# Automatic k-mer selection
$ spades.py --isolate -k auto ...

# Manual selection of k-mer sizes
$ spades.py --isolate -k 21,33,55,77,99 ...
```

Notice above that all manual *k*-mer values are odd-numbered. This is deliberate, as even-lengthed *k*-mers can yield palindromic *k*-mers which cannot be fitted into the assembly graph. Chosing only odd values is a simple way to prevent this issue.

In practice, we tend to find that `SPAdes` is quite quick to end the *k*-mer extension process and we can often get much better assemblies by manually specifying additonal *k*-mer sizes. However, the more *k*-mers we specify the longer assembly will take (as each *k*-mer size must go through a round of error correction and assembly). For an example of the effect of *k*-mer size on assembly, here is some data from a Genomics Aotearoa experiment assembling a marine metagenome:

|Method|*k*-mers assembled|Contigs assembled|Longest contig length|N50 (>2kbp)|L50 (>2kbp)|
|:---|:---|:---|:---|:---|
|Automatic|21, 33, 55|4,239,806|660,812|6,782|12,906|
|Manual|43, 55, 77, 99, 121|2,519,669|1,022,083|7,990|12,673|
|Manual|21, 43, 55, 77, 99, 121|3,388,682|71,022,083|7,789|13,327|

As you can see from this table, selecting additional *k*-mer sizes beyond what `-k auto` would chose resulted in a much larger longest contig. There were also modest improvement to the N50 and L50 of the assembly. However, this improvement came at a significant increase in run time.

### Checkpointing

As `SPAdes` assembly is a long and resource intensive process, the assembler creates checkpoints as it passes particular assembly milestones. This means that if a job is terminated and must be restarted, we do not need to repeat the full workflow from scratch. However, a restarted job is restricted to repeating with the same parameters as the original job so while checkpoints are a powerful quality of life improvement they do not allow us to branch our assembly process.

```bash
# Original job
$ spades.py --isolate -k 21,55,77 -1 ... -2 ... -o test_assembly/

# Restart the job
$ spades.py --continue -o test_assembly/
```

### Resource limitations

Finally, there are several parameters we can adjust with `SPAdes` to affect its resource utilisation. There are two options we will routinely need to change - the thread count, the memory ceiling.

##### Specifying the thread count

Like most of the tools we have used to date `SPAdes` is able to divide steps of its workflows onto separate compute cores allowing parts of the assembly to run in parallel with each other. By default the number of threads `SPAdes` will try to work with is 16, but this can be changed:

```bash
$ spades.py --threads 32 ...
```

##### Specifying the memory cap

`SPAdes` is notoriously memory-hungry when performing assembly, and therefore imposes an internal memory limitation on itslf to prevent itself from requesting more RMA thatn the computer can provide. This is a great feature, as it is better for `SPAdes` to end itself (crash) rather than cause the host computer to run out of memory and begin using its [swap space](https://opensource.com/article/18/9/swap-space-linux-systems). By default `SPAdes` will request up to 250 GB of RAM before terminating. This can be increased with the `--memory` parameter.

The interaction between `SPAdes` and `slurm` memory management is important to remember. If you detect that an assembly is crashing due to a lack of memory, you will need to increase the memory allocation both in your `slurm` file AND in `SPAdes` itself. For example, if a particular job needs 500 GB of RAM to complete we can encounter the following scenarios:

|slurm limit|SPAdes limit|Outcome|
|:---|:---|:---|
|250 GB|250 GB|Job terminates when memory request exceeds 250 GB|
|600 GB|250 GB|`SPAdes` self-terminates wen memory request exceeds 250 GB|
|250 GB|600 GB|`slurm` terminates `SPAdes` job when memory request execeeds 250 GB|
|600 GB|600 GB|Job completes successfully|

