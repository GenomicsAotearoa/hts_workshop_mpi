# 2.1 - Genome Assembly 

## Overview

!!! clock "time"

    * Teaching: 15 minutes
    
!!! circle-info "Key points"
   
    #### Keypoints
    
    * Assembly is the process of reconstructing long nucletic sequences (contigs) from fragmented read data.
    * There are two major approaches for how assembly is bets performed - Debruijn graph and Overlap Layout Consensus assemblies.
    * There are multiple tools for performing assembly, and specifically different tools for assembling short or long read data.

---

## The assembly process

Genome assembly is the act of organising sequencing data to produce a representation of the complete genome sequence of an organism (ie a plant, fungus, bacterium or virus). 
An assembled genome can then be used to perform down stream analysis such as annotation of coding regions, BLAST analysis or comparasion with known or reference material. 

The success of a genome assembly project depends on a range of factors including;

1. The type of genome being assembled (ie a viral genome is smaller and easier to assemble than a fungal genome) 
1. The qaulity of the original sample and the extracted DNA 
1. The sequencing technology used to generate the data
   1. Read length, read quality and read quantity can all impact genome assembly
1. The software used to assemble the genome 

---

## Assembly methods 

Although there is a wide range of genome assembly tools, two main approaches which genome assemblers use depending on the type of data available. These classes are: Overlap Layout Consensus (OLC) and De Bruijn Graph (DBG).   

Understanding how each assembly type works is generally beyond the scope of this training, but it is important to use an approriate assembly method for your data type. 

We're not going to dive deeply into the differences between these methods, but as a brief differentiation, OLC methods work by aligning the sequence reads against each other, and trying to build consensus sequences where sequence align and extend upon one another. By contrast DBG methods cut the input sequence data into smaller pieces, and maps the transition from one of these sub-reads to the next. Where common sub-reads are identified in different sequences, a map of how these sub-reads relate to each other is made and eventually a consensus sequence can be extracted.

The key important distinction between these methods is that OLC is robust to a degree of sequence variation between the reads it is trying to align, whereas DBG methods require very high quality data to keep the map of how sub-reads relate to each other manageable. Therefore, the key takeaway is that Overlap Layout Consensus methods are best suited to long read sequence data while De Bruijn Graph methods are better suited to short read sequence data. 

In this training we will use short and long read sequence data to familarise ourselves with two popular assemblers **SPAdes** (short read) and **Flye** (long read).

---
