# Level 2 revision

This revision material consists of five short exercises covering the main points of the Level 2 training. When you have completed all exercises, please contact one of the training representatives to have them check your work.

If you are struggling with any exercises, remember the following:

1. Every exercise here was covered in the training material. Refer to this for hints.
2. It is not cheating to use Google to find help on how various commands work.
3. You are welcome to contact the trainers for advice on a particular exercise, but please attempt the first two options before resorting to this.

---

## Exercise 1 - Modifying and transferring a text file

1. In your `fastq_processing/docs/` folder is a template file `module_2.1_homework.txt`. Use the `nano` text editor to modify the file and change the values (`#####`) where required.
1. Once completed, download the file through JupyterHub and email it to the trainer as part of your submission.

---

## Exercise 2 - Working with `slurm` scripts

Consider the `slurm` batch file below. It contains a number of errors. Identify at least 5 and send your observations, along with the correction to the trainer.

>**:Note:** As we have not yet covered running the `BLASTn` command from the command line, this is the same line that you used in the workshop. There are no intentional errors in this line.

```bash
#!/bin/bas -e

#SBATCH --account nesi03811
#SBATCH --job-name slurm_test
#SBATCH --time 00-30-00
#SBATCH --cpus-per-task2
#SBAYCH --output slurm_test.$j.out

# Loading the required software
module load blast
module purge

# Execute the command
blastn -db nt -query input.fasta -out output.blast.txt -outfmt 1 -evalue 1e-3
```

---

## Exercise 3 - Writing your own `slurm` script

For this task you need to create a small `slurm` script that can perform the following commands:

```bash
$ cd /nesi/project/nesi03181/phel/USERNAME/fastq_processing/data/
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/696/015/GCF_000696015.1_ASM69601v1/GCF_000696015.1_ASM69601v1_genomic.fna.gz
$ gunzip GCF_000696015.1_ASM69601v1_genomic.fna.gz
```

Create a `slurm` batch file to perform these commands. For the sake of efficiency, this job will only require 1 or 2 minutes to run and needs less than 1 GB of RAM to perform, so pick values that reflect this.

Once you are happy with your script, submit it using `sbatch` and report the *successful* job identifier to the trainers, and provide them with the location of the script and results in NeSI.

---

## Exercise 4 - Working with `fastp` and `FastQC`

In the `module_2/` resource folder there is a hidden directory `.homework_resources/` which contains the MiSeq output for two additional *M. bovis* isolates. Copy these to your `fastq_processing/data/` directory, then perform the following tasks:

1. Create a `FastQC` report of the sequences. Inspect any problems identified by `FastQC` and determine whether or not these are likely to be real issues for analysing the data.
1. Trim the sequences with `fastp`, removing any low quality regions and adapter sequences.
   * You can pick any value for the Q-score threshold.
   * Remember that these are paired-end data
1. Also process the `SRR097977` and `SRR098026` sequences with `fastp`
   * We know from our previous work with `grep` and `FastQC` that the quality of these mock libraries is much lower than a standard output. How many reads from each sample pass through your filtering? 

---

## Exercise 5 - Working with `pycoQC`

TO DO...
