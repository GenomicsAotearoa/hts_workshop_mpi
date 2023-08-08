## Genome Assembly 
---

Genome assembly is the act of organising sequencing data to produce a representation of the complete genome sequence of an organism (ie a plant, fungus, bacterium or virus). 
An assembled genome can then be used to perform down stream analysis such as annotation of coding regions, BLAST analysis or comparasion with known or reference material. 

The success of a genome assembly project depends on a range of factors including 
- The type of genome being assembled (ie a viral genome is smaller and easier to assemble than a fungal genome) 
- The qaulity of the original sample and the extracted DNA 
- The sequencing technology used to generate the data
  - Read length, read quality and read quantity can all impact genome assembly
- The software used to assemble the genome 

---

# Assembly methods 

Although there is a wide range of genome assembly tools, two main approaches which genome assemblers use depending on the type of data available. These classes are: Overlap Layout Consensus (OLC) and Debruijn Graph (DBG).   
Understanding how each assembly type works is generally beyond the scope of this training, but it is important to use an approriate assembly method for your data type. Therefore, it can be useful to understand that Overlap Layout Consensus methods are best suited to long read sequence data while Debruijn Graph methods are better suited to short read sequence data. 

In this training we will use short and long read sequence data to familarise ourselves with two popular assemblers **SPAdes** and **Flye**.

We will start with a SPades assembly of XXX data 

[Next lesson](SPAdes-assembly.md)

