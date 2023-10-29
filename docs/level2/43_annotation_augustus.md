# 4.3 - Prediction with `AUGUSTUS`

!!! clock "time"

    * Teaching: 10 minutes
    * Exercises: 30 minutes

!!! circle-info "Objectives and Key points"

    #### Objectives
    
    * Understand how to run `AUGUSTUS` for predicting protein coding regions in genomes

    #### Keypoints

    * Selecting the best prediction model is key to getting accurate predictions.
    * You can create your own training models if there are no closely related organisms for your species of interest.

---

## Getting started with the tool

Unlike prokaryotic genomes, the genes of eukaryotes carry intronic sequences which need to be spliced out of the gene sequence before undergoing translation. The detection of splicing boundaries is a difficult task, as there are many organism-specific patterns used to mark splice sites.

Protein prediction tools for this purpose typically come with  number of pre-trained models for finding protein domains within contigs, but if there is no model for your organism, or a closely related lineage, then results may not be ideal.

A recently published article ([Scalzitti *et al.*, 2020](https://doi.org/10.1186/s12864-020-6707-9)) which profiled a number of these tools found `AUGUSTUS` to be one of the best performing tools for gene prediction in eukaryotic organisms, so this is what we will use today.

>**Note:** `AUGUSUTUS` does require training against a closely related model organisms to generate accurate predictions, which we do not have for this workshop. We will instead be performing predictions with a few different models and seeing how the outputs differ.

!!! question "Exercise"

    Find the latest version of `AUGUSTUS` on NeSI and load it.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            module spider augustus

            module purge
            module load AUGUSTUS/3.5.0-gimkl-2022a

            augustus
            ```

            ??? success "Output"

                ```
                AUGUSTUS (3.5.0) is a gene prediction tool.
                Sources and documentation at https://github.com/Gaius-Augustus/Augustus

                usage:
                augustus [parameters] --species=SPECIES queryfilename

                'queryfilename' is the filename (including relative path) to the file containing the query sequence(s)
                in fasta format.

                SPECIES is an identifier for the species. Use --species=help to see a list.

                parameters:
                ...

                For a complete list of parameters, type "augustus --paramlist". A description of the important ones can be found in the file RUNNING-AUGUSTUS.md.
                ```

---

## Predicting protein coding regions

The sequence(s) we have to work with today are from the brown marmorated stink bug (*Halyomorpha halys*). Run `AUGUSTUS` with the prompt below to see what species models are available for prediction.

!!! terminal "code"

    ```bash
    augustus --species=help
    ```

    ??? success "Output"

        ```
        usage:
        augustus [parameters] --species=SPECIES queryfilename

        where SPECIES is one of the following identifiers

        identifier                               | species
        -----------------------------------------|----------------------
        pea_aphid                                | Acyrthosiphon pisum
        aedes                                    | Aedes aegypti
        amphimedon                               | Amphimedon queenslandica
        ancylostoma_ceylanicum                   | Ancylostoma ceylanicum
        adorsata                                 | Apis dorsata
        honeybee1                                | Apis mellifera
        arabidopsis                              | Arabidopsis thaliana
        ...
        (maize5)                                 | Zea mays
        ```

As you will see, there is no good model for *H. halys*, or even a closely related species from the *Pentatomidae*. 

We will use two different models for an initial round of prediction on the *Halyomorpha halys* sequences - one insect and one bacterial species. We will first run `AUGUSTUS` using the *Apis mellifera** (honey bee) model.

!!! terminal "code"

    ```bash
    augustus --genemodel=partial --protein=on --codingseq=on \
        --species=honeybee1 \
        input/NW_020110202.fna \
        > outputs/NW_020110202.aug_hb1.gff
    ```

This will take about 15 minutes to run, so while it is running set up the following exercise as well:

!!! question "Exercise"

    Select a second model organism, preferably a non-insect model, and begin a second round of prediction using `AUGUSTUS`.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            augustus --genemodel=partial --protein=on --codingseq=on \
                --species=E_coli_K12 \
                input/NW_020110202.fna \
                > outputs/NW_020110202.aug_ecoli.gff
            ```

---

## Extracting sequences from prediction files

Once your jobs have finished (about 15 minutes per attempt), we need to run a helper script that comes with `AUGUSTUS` to extract the gene and coding sequence predicitons from the output file.

!!! terminal "code"

    ```bash
    getAnnoFasta.pl outputs/NW_020110202.aug_hb1.gff
    ```

This will return no output to the console, but creates two new files, with names dervied from the file name above.

1. `*.codingseq` - Nucleotide sequences for each predicted coding region
1. `*.aa` - Amino acid residue translation from each predicted coding region

!!! question "Exercise"

    Use `grep` to count the number of protein coding sequences in your predicted files (both models) and the official annotation for this *H. halys* genome.

    ??? circle-check "Solution"

        !!! terminal "code"

            ```bash
            grep -c ">" outputs/NW_020110202.aug_hb1.aa outputs/NW_020110202.aug_ecoli.aa outputs/NW_020110202.canonical.faa 
            ```

            ??? success "Output"

                ```
                outputs/NW_020110202.aug_hb1.aa:170
                outputs/NW_020110202.aug_ecoli.aa:17
                outputs/NW_020110202.canonical.faa:110
                ```

Similar to when running predictions with `prodigal`, you can see that neither of the models identified the correct number of protein coding sequences. This is because neither of these gene models are accurate for the organism we are trying to characterise. However, the more closely related model was much closer in its prediction that the non-insect model.

---

## Creating a custom species profile for gene prediciton (optional)

Creating a new model is a slow process so we will not be running through it today, but if this is something you need to do in your own work then use the steps below to get started.

There are a few steps we need to perform in advance of the new model training. The first is to do with file permissions - the location of the prediction databases that `AUGUSTUS` uses for gene prediction are not writable to us, so we cannot add new data into them. We must create our own copy of the configuration information and point `AUGUSTUS` towards this new location in order to create new models

!!! terminal "code"

    ```bash
    cp -r ${AUGUSTUS_CONFIG_PATH} /your/location/to/store/data/

    # Update the `AUGUSTUS_CONFIG_PATH` variable to point to the new location
    AUGUSTUS_CONFIG_PATH="/your/location/to/store/data"
    ```

Now we just need to obtain a reference genome to train against. For *H. halys*, this can found on the [NCBI website](https://www.ncbi.nlm.nih.gov/assembly/GCA_000696795.3) and downloaded from the command line:

!!! terminal "code"

    ```bash
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/696/795/GCA_000696795.3_Hhal_1.1/GCA_000696795.3_Hhal_1.1_genomic.gbff.gz
    gunzip GCA_000696795.3_Hhal_1.1_genomic.gbff.gz
    ```

Once these steps are completed, training is a single command:

!!! terminal "code"

    ```bash
    autoAugTrain.pl --trainingset=GCA_000696795.3_Hhal_1.1_genomic.gbff --species=hhalys
    ```

We could then perform prediction using the new model as usual. Note that if this was performed in a new session, we would need to set the `AUGUSTUS_CONFIG_PATH` variable again.

!!! terminal "code"

    ```bash
    augustus --AUGUSTUS_CONFIG_PATH=/your/location/to/store/data --species=hhalys --genemodel=partial ...
    ```

---
