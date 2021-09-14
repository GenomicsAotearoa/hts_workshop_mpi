# Overview

This is the opening page for the HTS workshop. Content is divided into four levels of increasing complexity, and will be added throughout the year as materials are developed and the training program advances.

---

## Contents

1. Level 1 - Shell genomics
   * [Background](#level-1---shell-genomics)
   * [Common terms and references](common_terms.md)
   * [Materials](level1/01-introduction.md)
   * [Homework](level1/05-homework.md)
1. Level 2 - Working with scripts
   * [Background](#level-2---slurm-and-quality-filtering)
   * [`slurm` and `module` cheatsheet](slurm_module_guide.md)
   * [Materials (Part 1)](level2/01-writing-scripts.md)
   * [Materials (Part 2)](level2/03-intro-to-nesi.md)
   * [Homework](level2/06-homework.md)
1. Level 3 - Read mapping and assembly
   * [Background](#Level-3---read-mapping-and-assembly)
   * [Materials (part 1)](level3/01-introduction-to-mapping.md)
   * [Homework (part 1)](level3/04-homework.md)
   * [Materials (part 2)](level3/05-illumina-assembly.md)
   * [Homework (part 2)](level3/09-homework.md)
   * [Materials (part 3.1)](level3/10-gene-prediction.md)

---

## Level 1 - Shell genomics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260560.svg)](https://doi.org/10.5281/zenodo.3260560)

An introduction to the Unix shell for people working with genomics data. This material is adapted from the [Data Carpentry Genomics Workshop](http://www.datacarpentry.org/genomics-workshop/). Please see http://www.datacarpentry.org/shell-genomics/ for the original version of this material.

Command line interface (OS shell) and graphic user interface (GUI) are different ways of interacting with a computer's operating system. The shell is a program that presents a command line interface which allows you to control your computer using commands entered with a keyboard instead of controlling graphical user interfaces (GUIs) with a mouse/keyboard combination.

There are quite a few reasons to start learning about the shell:

- For most bioinformatics tools, you have to use the shell. There is no graphical interface. If you want to work in metagenomics or genomics you're going to need to use the shell.
- The shell gives you power. The command line gives you the power to do your work more efficiently and more quickly. When you need to do things tens to hundreds of times, knowing how to use the shell is transformative.
- To use remote computers or cloud computing, you need to use the shell.

This lesson is part of a workshop that is run through the [New Zealand eScience Infrastructure](https://www.nesi.org.nz/) (NeSI) platforn. Workshop participants are expected to have set up an account with the correct project access prior to attending this workshop. Please contact the workshop organisers to arrange access to the project accounts.

The starting data for these exercises can be obtained from the Data Carptentry Genomics Workshop [setup webpage](https://datacarpentry.org/genomics-workshop/setup.html) or from [this link](https://ndownloader.figshare.com/files/14417834) directly.

### Authors

Shell Genomics is authored and maintained by the [community](https://github.com/datacarpentry/shell-genomics/network/members). Specific edits have been made to make the materials consistent with the MPI software environment and working on NeSI.

### Citation

The original work can be cited as:

Erin Alison Becker, Anita Schürch, Tracy Teal, Sheldon John McKay, Jessica Elizabeth Mizzi, François Michonneau, et al. (2019, June). datacarpentry/shell-genomics: Data Carpentry: Introduction to the shell for genomics data, June 2019 (Version v2019.06.1). Zenodo. [http://doi.org/10.5281/zenodo.3260560](http://doi.org/10.5281/zenodo.3260560)

---

## Level-2 - `slurm` and quality filtering

This workshop provides a basic introduction to working with the `slurm` scheduling system, and begins working with Illumina MiSeq and Oxford Nanopore Technology sequence data. The data used in this workshop is from a real study - available at the [NCBI website](https://www.ncbi.nlm.nih.gov/) under `BioProject` number [PRJEB38523](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB38523).

In particular, we are working with the following data samples:

1. [Run ERR4179765](https://www.ncbi.nlm.nih.gov/sra/ERX4143189[accn]) - MinION sequences from *Mycoplasma bovis* strains, obtained using the `guppy` basecalling software
1. [Run ERR4179828](https://www.ncbi.nlm.nih.gov/sra/ERX4143252[accn]) - Illumina MiSeq sequences from the same strains

Additional teaching materials were sourced from the [Genomics Aoteoroa Metagenomic Summer School workshop](https://github.com/GenomicsAotearoa/metagenomics_summer_school) and the [Long-Read, long reach Bioinformatics Tutorial](https://timkahlke.github.io/LongRead_tutorials/) put together by Tim Kahlke.

---

## Level 3 - Read mapping and assembly

This workshop covers the process of mapping Illumina and Oxford Nanopore sequences against reference genomes, and covers the basics of performing *de novo* assembly.

Some material has been adapted from the [Genomics Aotearoa](https://www.genomics-aotearoa.org.nz/) metagenomics summer school ([github here](https://github.com/GenomicsAotearoa/metagenomics_summer_school)), and the training materials provided by the [Galaxy community](https://training.galaxyproject.org/training-material/).

This module also uses sequence data from `BioProject` [PRJNA678533]:

1. [Run SRR13090255](https://www.ncbi.nlm.nih.gov/sra/SRX9536177[accn]) - DNA-Seq of *Halyomorpha halys* (H1 Haplotype)

### Citation

Joachim Wolff, Bérénice Batut, Helena Rasche, (2021). Mapping (Galaxy Training Materials). [https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/mapping/tutorial.html](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/mapping/tutorial.html).
