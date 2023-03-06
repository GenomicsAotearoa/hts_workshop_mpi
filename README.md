# Overview

This is the opening page for the HTS workshop. Content is divided according to the teaching stream depending on user experience and profficiency.

---

## Contents

1. Level 1 - Beginner
   1. [Connecting to NeSI](./level1/11_nesi_connection.md)
   1. [Navigation on the command line](./level1/12_shell_navigation.md)
   1. [Working with files on the command line](./level1/13_shell_manipulation.md)
   1. [Quality filtering Nanopore data](./level1/2X_quality_filter_nanopore.md)
   1. [Quality filtering Illumina data](./level1/3X_quality_filter_illumina.md)
   1. [Annotating sequences with BLAST](./level1/4X_blastn_annotation.md)
1. Level 2 - Advanced
   1. [Shell navigation (advanced)](./level2/01_shell_manipulation.md)
   1. [Loops and variables in the command line](./level2/02_shell_variables.md)
   1. [Redirection in the command line](./level2/03_shell_redirection.md)
   1. [*De novo* assembly of sequencing data](./level2/04_assembly_de_novo.md)
   1. [Polishing of genome assemblies](./level2/05_assembly_polishing.md)
   1. [Mapping reads to a reference](./level2/06_coverage_mapping.md)
   1. [Calculating coverage statistics](./level2/07_coverage_statistics.md)
   1. [Introduction to workflow managers](./level2/08_workflows_introduction.md)
   1. [Applying a workflow rule to existing data](./level2/09_workflows_applying.md)

---

## Getting started

The work covered in this training programme are run through the [New Zealand eScience Infrastructure](https://www.nesi.org.nz/) (NeSI) platforn. Workshop participants are expected to have set up an account with the correct project access prior to attending this workshop. Please contact the workshop organisers to arrange access to the project accounts.

If you are a beginner to this work, keep in mind that the [glossary of terms](./docs/common_terms.md) and [slurm module guide](./docs/slurm_module_guide.md) which will be helpful as we progress through the materials.

---

## Background - Shell genomics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260560.svg)](https://doi.org/10.5281/zenodo.3260560)

An introduction to the Unix shell for people working with genomics data. This material is adapted from the [Data Carpentry Genomics Workshop](http://www.datacarpentry.org/genomics-workshop/). Please see http://www.datacarpentry.org/shell-genomics/ for the original version of this material.

Command line interface (OS shell) and graphic user interface (GUI) are different ways of interacting with a computer's operating system. The shell is a program that presents a command line interface which allows you to control your computer using commands entered with a keyboard instead of controlling graphical user interfaces (GUIs) with a mouse/keyboard combination.

There are quite a few reasons to start learning about the shell:

* For most bioinformatics tools, you have to use the shell as there is no graphical interface.
* The shell gives you power to do your work more efficiently and more quickly.
  * When you need to do things tens to hundreds of times, knowing how to use the shell is transformative.

Many of the exercises covered in this training programme are obtained from, or inspired by, the [Data Carptentry](https://datacarpentry.org/) initiative, particularly their [Genomics Workshop](https://datacarpentry.org/genomics-workshop/setup.html)[^1].

---

## Background - Data used in training

This workshop provides a basic introduction to working with the `slurm` scheduling system, and begins working with Illumina MiSeq and Oxford Nanopore Technology sequence data. The data used in this workshop is from a real study - available at the [NCBI website](https://www.ncbi.nlm.nih.gov/) under `BioProject` number [PRJEB38523](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB38523).

In particular, we are working with the following data samples:

1. [Run ERR4179765](https://www.ncbi.nlm.nih.gov/sra/ERX4143189[accn]) - MinION sequences from *Mycoplasma bovis* strains, obtained using the `guppy` basecalling software
1. [Run ERR4179828](https://www.ncbi.nlm.nih.gov/sra/ERX4143252[accn]) - Illumina MiSeq sequences from the same strains
1. [Run SRR13090255](https://www.ncbi.nlm.nih.gov/sra/SRX9536177[accn]) - DNA-Seq of *Halyomorpha halys* (H1 Haplotype)
1. [*H. halys* reference genome](https://www.ncbi.nlm.nih.gov/assembly/GCA_000696795.3) - Hhal_1.1 assembly

Additional teaching materials were sourced from:

1. Genomics Aoteoroa Metagenomic Summer School workshop[^2].
1. Long-Read, long reach Bioinformatics Tutorial[^3].
1. Galaxy Training! seuqence analysis resources[^4].

---

## Citations

[^1]: Erin Alison Becker, Anita Schürch, Tracy Teal, Sheldon John McKay, Jessica Elizabeth Mizzi, François Michonneau, *et al.* (2019, June). datacarpentry/shell-genomics: Data Carpentry: Introduction to the shell for genomics data, June 2019 (Version v2019.06.1). Zenodo. [http://doi.org/10.5281/zenodo.3260560](http://doi.org/10.5281/zenodo.3260560).

[^2]: Jian Sheng Boey, Dinindu Senanayake, Michael Hoggard *et al.* (2022). Metagenomics Summer School [https://github.com/GenomicsAotearoa/metagenomics_summer_school](https://github.com/GenomicsAotearoa/metagenomics_summer_school).

[^3]: Tim Kahlke (2021). Long-Read Data Analysis [https://timkahlke.github.io/LongRead_tutorials/](https://timkahlke.github.io/LongRead_tutorials/).

[^4]: Joachim Wolff, Bérénice Batut, Helena Rasche (2023). Sequence Analysis (revision 96e01807afff10d6060ac0691d004f0469676534). [https://training.galaxyproject.org/training-material/topics/sequence-analysis/](https://training.galaxyproject.org/training-material/topics/sequence-analysis/).
