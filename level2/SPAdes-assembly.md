## Introduction to the SPAdes assembly tool

The `SPAdes` assembler is a very powerful and popular assembler which utilises de Bruijn graphs to assemble short read sequence data into larger contigs.

SPAdes is particularly good at: 
* *De novo* assembly of short read sequences (i.e. Illumina or IonTorrent)
* Assembling small genomes (ideally <100 Mb, such as bacterial, viral, fungal, mitochondrial genomes)
  * Despite this, it is possible to apply `SPAdes` to larger genomes and obtain very good results.


Before beginning to work with `SPAdes` we need to ensure that our data is free from adapter sequences. The main reason for doing this is when we assemble sequences to form a genome, the assembler is looking for spans of nucleic acid sequence which are common to multiple reads, so that those reads can be joined together to create longer contigs. As the adapater sequence is an identical tag added to every read, these create regions of artificial homology between sequences which have no real connection to each other. Before we attempt assembly it is critical to remove these from our data.

In order to begin, we must first find the versions of `SPAdes` installed on NeSI and load the module of interest.

```bash
$ module load SPAdes/3.15.2-gimkl-2020a

$ spades.py -h
```

---

## Performing an assembly using SPAdes

The test data is a set of Illumina MiSeq sequencing reads from a BMSB sample, specifically isolated from mitochondria. To save time, the reads have already been quality filtered with `fastp` and the number of reads reduced to speed up analysis.

Copy the BMSB reads to your SPAdes working directory

```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/2_Quality_filtered_data/
$ cp /nesi/project/nesi03181/phel/module_3/2_Quality_filtered_data/SRR13090255_*.fq.gz ./

$ cd ../
```

Create a `slurm` script with the following contents. Be sure to replace the `YOUR_EMAIL` and `USERNAME` values with your details.

```bash
#!/bin/bash -e
#SBATCH --account       nesi03181
#SBATCH --job-name      bmsb_spades
#SBATCH --time          00:10:00
#SBATCH --cpus-per-task 10
#SBATCH --mem           4G
#SBATCH --error         bmsb_spades.%j.err
#SBATCH --output        bmsb_spades.%j.out
#SBATCH --mail-type     END
#SBATCH --mail-user     YOUR_EMAIL

module purge
module load  SPAdes/3.15.2-gimkl-2020a

# Set the path from which the script will execute SPAdes
cd /nesi/project/nesi03181/phel/USERNAME/

# Execute SPAdes
spades.py --threads 10 \
          -1 2_Quality_filtered_data/SRR13090255_1.fq.gz \
          -2 2_Quality_filtered_data/SRR13090255_2.fq.gz \
          -s 2_Quality_filtered_data/SRR13090255_s.fq.gz \
          -o 3_Assembly-mapping/bmsb_spades/ 
```
Submit this job to `slurm`:

```bash
$ sbatch bmsb_spades.sl
```
