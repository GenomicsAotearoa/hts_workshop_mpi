# 2023 Level 2 proficiency testing (version 1)

!!! clock "time"

    * Exercises: 150 minutes

!!! circle-info "Key information"

    * The aim is to complete each of the exercises before the end of the session.
    * This testing is open book - you may use any materials from the level 1 or level 2 training and anything online.
    * This testing is individual - please do not discuss your solutions with other trainees.
    * The tutors are here to help with significant technical issues such as loss of connection. *We are not here to help you identify tools or how to solve specific questions.*

    All test folders are located in the following location:

    !!! terminal "code"

        ```bash
        /nesi/project/nesi03181/phel/proficiency/USERNAME/
        ```

!!! warning "Output files"

    For each exercise, you are free to name your output folders/files whatever you wish.

    Please ensure that your results are in the appropriate folder for each exercise, however.

!!! warning "Reporting findings"

    Whenever an exercise asks your to in some way interpret an output file and report your results, just create a small text file to store your notes.

    Name the file according to the exercise that asked for the information.

---

## Exercise 1

Navigate to your instance of the `01_illumina_assembly/` folder. You have been provided with the following files:

!!! terminal "code"

    ```
    reads/illumina_R1.fq.gz
    reads/illumina_R2.fq.gz
    reference/reference.fna
    results/provided_assembly.fna
    ```

!!! question "Part 1"

    You are to produce a `slurm` script that uses an approriate assembly tool to assemble the Illumina MiSeq sequences in the `reads/` folder into a draft genome.

    >Your script will only need about 10 minutes to complete with 16 CPUs.

!!! question "Part 2"

    Use an appropriate tool to determine assembly statistics for your assembly file and use the file `reference/reference.fna` as your reference genome.

    >If your job does not complete in a reasonable amount of time, you can use the provided `results/provided_assembly.fna` file instead of your own output.

    Report the folllowing metrics from your assembly:

    1. Number of contigs with length greater than or equal than 50 kbp.
    1. The N50 value for the assembly.
    1. Number of indels relative to the reference genome.

---

## Exercise 2

Navigate to your instance of the `02_mapping/` folder. You have been provided with the following files:

!!! terminal "code"

    ```
    reads/nanopore.fq.gz
    reference/reference.fna
    results/provided_mapping.sam
    ```

!!! question "Part 1"

    Use whatever tool is appropriate to map the sequences in the `reads/` folder against the reference file, creating a **sam** file with the mapping results.

!!! question "Part 2"

    1. Once you have completed mapping your reads, sort and compress the mapping information from your output **sam** file.

    >If you are unable to complete your mapping job for any reason, use the file supplied at `results/provided_mapping.sam.

    2. Produce a depth report from the compressed mapping file.

    3. When you have produced your compressed file, filter the results to keep only the **mapped** reads in the mapping output.

---

## Exercise 3

Navigate to your instance of the `03_gene_calling/` folder. You have been provided with the following file:

!!! terminal "code"

    ```
    input_sequences.fna
    ```

!!! question "Question"

    You have been provided with a draft assembly of a prokaryotic genome. Create predictions of the protein coding regions of the assembly using whichever tool you believe is appropriate.

---

## Exercise 4

Navigate to your instance of the `04_classification_kraken2/` folder. You have been provided with the following files:

!!! terminal "code"

    ```
    inputs/input_sequences.fna
    results/provided_kraken2.txt
    ```

!!! question "Part 1"

    Write a `slurm` script to execute a `kraken2` classification of the sequences provided in `inputs/input_sequences.fna`.

    When selecting a database, use the `PlusPFP` database located at:

    !!! terminal "code"

        ```
        /nesi/project/nesi03181/phel/databases/k2_pluspfp_16gb_20231009/
        ```

    Check your notes carefully to make sure you understand how to direct `kraken2` towards this database.

!!! question "Part 2"

    You have been provided with an output file of the above classification job - `results/provided_kraken2.txt`. Use your knowledge of the command line to search through the output file and identify the most likely species to be in this sample.

    Create a text file containing your findings. You do not need to report every species hit in the `provided_kraken2.txt` file - just pick the candidates you believe are most likely and give a quick 1-sentence justification for each choice.

---

## Exercise 5

Navigate to your instance of the `05_classification_diamond/` folder. You have been provided with the following file:

!!! terminal "code"

    ```
    inputs/input_sequences.faa
    results/provided_diamond.txt
    ```

!!! question "Part 1"

    Write a `slurm` script to execute a `diamond` classification of the sequences provided in `inputs/input_sequences.faa`.

    When selecting a database, use the `uniprot_sprot.dmnd` database located at:

    !!! terminal "code"

        ```
        /nesi/project/nesi03181/phel/databases/swissprot_dmnd/uniprot_sprot.dmnd
        ```

!!! question "Part 2"

    Examine the outputs of the file `results/provided_diamond.txt`, remembering the [BLAST6 format](https://genomicsaotearoa.github.io/hts_workshop_mpi/level2/44_annotation_protein/#comparing-the-outputs) for column meaning.

    Based on the top hit results, determine the most likely genus or species of the organism from which these sequences were obtained. Create a text file containing your conclusion.

---


## Closing out the test

When you have completed all of the exercises, and you are happy with your results, please run the following command to log the history of your command line to a file.

!!! terminal "code"

    ```bash
    history > /nesi/project/nesi03181/phel/proficiency/USERNAME/history_log.txt
    ```

---
