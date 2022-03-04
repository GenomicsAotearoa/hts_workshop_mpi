# Annotation of sequences

* Teaching: 30 minutes
* Exercises: 30 minutes

#### Objectives

* Understand the differences betweem nucleotide and protein sequence matching
* Be aware of which publicly available databases are appropriate for which data
* Learn how to submit `BLAST` jobs to the NeSI cluster using `slurm`

#### Keypoints

* Annotation is required in order to identify the function and origin of sequences obtaining from HTS analysis
* There are different databases available for annotation and classification
* Interpreting the results of BLAST alignments can be a tricky process
  * Understand the meanings of the identity, coverage, e-value, and bitscore metrics when assessing a BLAST output

---

## Contents

1. [The `BLAST` process](#the-blast-process)
1. [Submitting a nucleotide `BLAST` job on NeSI](#submitting-a-nucleotide-blast-job-on-nesi)
1. [Submitting a protein `BLAST` job on NeSI](#submitting-a-protein-blast-job-on-nesi)
1. [Interpretting the results of BLAST queries](#interpretting-the-results-of-blast-queries)
---

## The `BLAST` process

The fundamental unit of sequence comparison is the 'Blast Local Alignment Search Tool' (BLAST) method, originally developed in 1990 ([Altschul *et al.*, 1990](https://doi.org/10.1016/S0022-2836(05)80360-2)). Without going into excessive detail, the approach is to take a set of novel (query) sequences and compare each sequence against a reference database of sequences with a known annotation of taxonomic origin (target sequences). From sufficiently well-matched pairings between query and target, we infer that the origin of the sequence, and possibly function of the gene, are shared.

`BLAST` works by breaking each query sequence into a set of smaller *seed* sequences, and searching each target in the database for the presence of these seeds. Where matches are found, the tool then extends the ends of the seed and assesses how well the ends of the query beyond the seed match the corresponding regions of the target. The quality of the match between the query and the target are evaluated in terms of how well conserved the sequence composition is between the pair, as well as how many insertion or deletion events need to be introduced to maintain the matching.

<img src='../img/12_blast.png' alt='Toy example of the BLAST process in action' width='600' />

In the example above, we are trying to match a query sequence `TMATO` against a database of four vegetable names, with the following process

a. The seed sequence `ATO` is identified in the query and two target sequences (green).
b. The query/target pairs are aligned and the character upstream from `ATO` is examined. Matches (blue) and mismatches (orange) are recorded.
c. The process repeats at the next position, scoring mismatches for both targets.
d. The final position is examined. For the target sequence `TOMATO` it is possible to insert a gap into the query sequence to preserve positional homology between query and target and improve the final score, but for the `POTATO` target the number of matching positions is never increased from the initial three.

Over the course of this matching, every target in the database will be assessed for it's suitability to match with the query sequence, and results ranked by how well they match. Typically, we restrict our results to only return a certain number of best matches, rather than report everything with any degree of similarity to the query. Regardless of how many results are returned, for a good set of matches we would expect to see a strong consensus in the gene function and taxonomy of the top hits for our query. From this consensus we can make inferences about the role of the sequence and which organism from which it may have originated.

---

## Submitting a nucleotide `BLAST` job on NeSI

...

---

## Submitting a protein `BLAST` job on NeSI

...

---

## Interpretting the results of BLAST queries

...

---
