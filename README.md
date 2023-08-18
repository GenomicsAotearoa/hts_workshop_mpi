# Overview

This is the opening page for the HTS workshop. Content is divided according to the teaching stream depending on user experience and profficiency.

---

## Contents

   1. Introduction to NeSI
      1. [Connecting to NeSI](docs/level1/11_nesi_connection.md)
      1. [Navigation on the command line](docs/level1/12_shell_navigation.md)
      1. [Working with files on the command line](docs/level1/13_shell_manipulation.md)
   1. Quality filtering Illumina data
      1. [Inspecting reads](docs/level1/21_illumina_inspection.md)
      1. [Trimming paired-end reads](docs/level1/22_illumina_filtering.md)
   1. [Introduction to slurm](docs/level1/31_slurm_introduction.md)
   1. [Quality filtering Nanopore data](docs/level1/32_quality_filter_nanopore.md)
   1. [Annotating sequences with BLAST](docs/level1/4X_blastn_annotation.md)
1. Level 2 - Advanced
   1. [Shell navigation (advanced)](docs/level2/11_shell_manipulation.md)
   1. [Loops and variables in the command line](docs/level2/12_shell_variables.md)
   1. [Redirection in the command line](docs/level2/13_shell_redirection.mdd)
   1. [*De novo* assembly of sequencing data](docs/level2/21_assembly_de_novo.md)
   1. [Polishing of genome assemblies](docs/level2/22_assembly_polishing.md)
   1. [Mapping reads to a reference](docs/level2/31_coverage_mapping.md)
      1. [Illumina mapping](docs/level2/32_illumina_mapping.md)
      1. [Nanopore mapping](docs/level2/33_nanopore_mapping.md)
      1. [Filtering and sorting mapping files](docs/level2/34_mapping_filters.md)
      1. [Summarising mapping statistics](docs/level2/35_mapping_statistics.md)
   1. [Introduction to workflow managers](docs/level2/41_workflows_introduction.md)
      1. [Building a basic workflow](docs/level2/42_workflow_starting.md)
      1. [Extending an existing workflow](docs/level2/43_workflow_extending.md)

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

This workshop provides a basic introduction to working with the `slurm` scheduling system, and begins working with Illumina MiSeq and Oxford Nanopore Technology sequence data. The data used in this workshop is mostly using simulated reads, produced using `InSilicoSeq`[^2] from the *Mycoplasma bovis* 8790 reference genome [NZ_LAUS01000004.1](https://www.ncbi.nlm.nih.gov/nuccore/NZ_LAUS01000004.1). We also make use of publicly available sequencing data from the studies [PRJNA813586](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA813586), [PRJEB38441](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB38441), and [PRJEB38523](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB38523).

Additional teaching materials were sourced from:

1. Genomics Aoteoroa Metagenomic Summer School workshop[^3].
1. Long-Read, long reach Bioinformatics Tutorial[^4].
1. Galaxy Training! seuqence analysis resources[^5].

---

## Citations

[^1]: Erin Alison Becker, Anita Schürch, Tracy Teal, Sheldon John McKay, Jessica Elizabeth Mizzi, François Michonneau, *et al.* (2019, June). datacarpentry/shell-genomics: Data Carpentry: Introduction to the shell for genomics data, June 2019 (Version v2019.06.1). Zenodo. [http://doi.org/10.5281/zenodo.3260560](http://doi.org/10.5281/zenodo.3260560).

[^2] Hadrien Gourlé, Oskar Karlsson-Lindsjö, Juliette Hayer, Erik Bongcam-Rudloff (2019). Simulating Illumina metagenomic data with InSilicoSeq. Bioinformatics 35(3), 521-522.

[^3]: Jian Sheng Boey, Dinindu Senanayake, Michael Hoggard *et al.* (2022). Metagenomics Summer School [https://github.com/GenomicsAotearoa/metagenomics_summer_school](https://github.com/GenomicsAotearoa/metagenomics_summer_school).

[^4]: Tim Kahlke (2021). Long-Read Data Analysis [https://timkahlke.github.io/LongRead_tutorials/](https://timkahlke.github.io/LongRead_tutorials/).

[^5]: Joachim Wolff, Bérénice Batut, Helena Rasche (2023). Sequence Analysis (revision 96e01807afff10d6060ac0691d004f0469676534). [https://training.galaxyproject.org/training-material/topics/sequence-analysis/](https://training.galaxyproject.org/training-material/topics/sequence-analysis/).
