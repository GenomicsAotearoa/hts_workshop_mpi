# Peforming hybrid assembly of Illumina and Oxford Nanopore sequence data (optional)

* Teaching: 30 minutes

#### Objectives

* Understand the different strategies which can be taken when combining short and long read sequence data

#### Keypoints

* ...
* ...

---

## Contents

1. [Options for hybrid assembly](#options-for-hybrid-assembly)
   1. [Polish then assemble](#polish-then-assemble)
   1. [Assemble then polish](#assemble-then-polish)
   1. [Hybrid assembly](#hybrid-assembly)

---

## Options for hybrid assembly

Both short read sequence data, such as Illumina HiSeq and MiSeq data, and long read sequences (Oxford Nanopore, PacBio) have their advantages over the other type of sequencing. However, it is completely reasonable to want to obtain the best of both worlds and combine sequencing data from multiple approaches to obtain a superior result.

You should by now be familiar with the following advantages/disadvantages of long and short read sequencing, namely that short read sequencing yields much higher overall sequence quality, but long reads enable assembly through repeat regions, or low complexity regions, and can result in longer contigs than what high-quality short reads can provide. If we want to combine these methods to get the advantages from both, then there are three broad approaches we can take:

1. Polish the raw long reads with short sequences, then assemble as normal
1. Assemble the long reads as normal, then polish the assembled contigs with the short read data
1. Perform short read assembly, with the long reads provided as reference material for scaffolding contigs

We have already briefly mentioned the third option during our walk through of `SPAdes` commands ([Exercise 6, Input data and reference files](06-assembly-choices.md#input-data-and-reference-files)) so we will start with the first options.

>**Note:** Before continuing, please keep in mind that this is an actively growing area of research and new results are being published all the time. The following results and recommendations may not prove true as sequencing error rates improve, new tools are released, and existing tools are refined.

### Polish then assemble

This approach was benchmarked in a study by [Fu and colleages (2019)](https://doi.org/10.1186/s13059-018-1605-z), who compared 10 polishing tools against a data set consisting of one bacterial, fungal, invertebrate, and plant genome. Profiles were made from raw and corrected sequences and assembly efficiency was evaluated for all analysis methods. They summarised their results in the following figure, noting several notable high-performing options.

![Fu *et al.* (2019). Figure 6](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-018-1605-z/MediaObjects/13059_2018_1605_Fig6_HTML.png?as=webp)

From this figure, the following tools appear to perform the best for polishing long read data:

1. `FMLRC` ([Wang *et al.*, 2018](https://doi.org/10.1186/s12859-018-2051-3))
1. `Jabba` ([Miclotte *et al.*, 2015](https://doi.org/10.1007/978-3-662-48221-6_13))
1. `LoRDEC` ([Salmela & Rivals, 2014](https://doi.org/10.1093/bioinformatics/btu538))

However, it is important to note the caveats of these findings. Most importantly, the plots display several profiles for each tool, reporting how the performance of each tool varies with depth of target coverage. In the case of `FMLRC` this is probably the best performing tool in the study *at high coverage depth*, but as coverage reduces the sensitivity is impacted much more severely than that of `LoRDEC` which boasts a lower maximum sensitivity but is more robust to a decreasing depth of coverage.

It is also worth noting that in their data they report that assemblies were generally improved by performing polishing but again, results varied with coverage. The following table is comprised of data summarised from [Fu *et al.* (2019), Table S7](https://doi.org/10.1186/s13059-018-1605-z):

<table>
  <th>
    <td>Method</td>
    <td>Fraction assembled 5x</td>
    <td>10x</td>
    <td>50x</td>
    <td>Sequence accuracy 5x</td>
    <td>10x</td>
    <td>50x</td>
  </th>
  <tr>
    <td><i>E. coli</i> (ONT)</td>
    <td>Raw</td>
    <td colspan=3>1.02</td>
    <td colspan=3>0.82</td>
  </tr>
  <tr>
    <td></td>
    <td>FMLRC</td>
    <td>1.02</td>
    <td>1.00</td>
    <td>1.00</td>
    <td>0.84</td>
    <td>0.97</td>
    <td>0.96</td>
  </tr>
  <tr>
    <td></td>
    <td>Jabba</td>
    <td>N/A</td>
    <td>N/A</td>
    <td>0.46</td>
    <td>N/A</td>
    <td>N/A</td>
    <td>0.99</td>
  </tr>
  <tr>
    <td></td>
    <td>LoRDEC</td>
    <td>1.01</td>
    <td>1.00</td>
    <td>1.00</td>
    <td>0.88</td>
    <td>0.96</td>
    <td>0.96</td>
  </tr>
  <tr>
    <td><i>E. coli</i> (PacBio)</td>
    <td>Raw</td>
    <td colspan=3>1.05</td>
    <td colspan=3>0.87</td>
  </tr>
  <tr>
    <td></td>
    <td>FMLRC</td>
    <td>1.04</td>
    <td>1.01</td>
    <td>1.00</td>
    <td>0.88</td>
    <td>0.98</td>
    <td>0.99</td>
  </tr>
  <tr>
    <td></td>
    <td>Jabba</td>
    <td>N/A</td>
    <td>N/A</td>
    <td>0.37</td>
    <td>N/A</td>
    <td>N/A</td>
    <td>0.99</td>
  </tr>
  <tr>
    <td></td>
    <td>LoRDEC</td>
    <td>1.02</td>
    <td>1.02</td>
    <td>1.02</td>
    <td>0.94</td>
    <td>0.99</td>
    <td>0.99</td>
  </tr>
  <tr>
    <td><i>S. cerevisae</i> (ONT)</td>
    <td>Raw</td>
    <td colspan=3>0.95</td>
    <td colspan=3>0.83</td>
  </tr>
  <tr>
    <td></td>
    <td>FMLRC</td>
    <td>0.93</td>
    <td>0.97</td>
    <td>0.98</td>
    <td>0.84</td>
    <td>0.95</td>
    <td>0.97</td>
  </tr>
  <tr>
    <td></td>
    <td>Jabba</td>
    <td>N/A</td>
    <td>N/A</td>
    <td>0.00</td>
    <td>N/A</td>
    <td>N/A</td>
    <td>0.91</td>
  </tr>
  <tr>
    <td></td>
    <td>LoRDEC</td>
    <td>0.95</td>
    <td>0.97</td>
    <td>0.98</td>
    <td>0.87</td>
    <td>0.95</td>
    <td>0.95</td>
  </tr>
  <tr>
    <td><i>S. cerevisae</i> (PacBio)</td>
    <td>Raw</td>
    <td colspan=3>1.01</td>
    <td colspan=3>0.85</td>
  </tr>
  <tr>
    <td></td>
    <td>FMLRC</td>
    <td>1.01</td>
    <td>1.01</td>
    <td>1.01</td>
    <td>0.85</td>
    <td>0.85</td>
    <td>0.85</td>
  </tr>
  <tr>
    <td></td>
    <td>Jabba</td>
    <td>0.00</td>
    <td>N/A</td>
    <td>0.00</td>
    <td>1.00</td>
    <td>N/A</td>
    <td>0.99</td>
  </tr>
  <tr>
    <td></td>
    <td>LoRDEC</td>
    <td>0.98</td>
    <td>1.02</td>
    <td>1.01</td>
    <td>0.90</td>
    <td>0.96</td>
    <td>0.96</td>
  </tr>
</table>

### Assemble then polish

As you might note, the order of polshing prior to assembly is the opposite of what we performed in the last exercise, and the `racon` and `medaka` documentation specifically detail polishing an assembly not the raw reads. This is generally the approach tested in the development of polishing tools. For example, along with the two tools we have worked with, the following polish software describes a polish/assemble workflow:

1. `POLCA` ([Zimin *et al.*, 2013](https://doi.org/10.1093/bioinformatics/btt476))
1. `Pilon` ([Walker *et al.*, 2014](https://doi.org/10.1371/journal.pone.0112963))
1. `NanoPolish` ([Loman *et al.*, 2015](https://doi.org/10.1038/nmeth.3444))
1. `HomoPolish` ([Huang *et al.*, 2021](https://doi.org/10.1186/s13059-021-02282-6))

The advantage of polishing post-assembly is that it allows multiple tools to be tested on a single draft assembly to compare their performance. This approach also carries the benefit that the benchmarking of assemblers is performed on uncorrected reads. If correction is performed prior to assembly the results may be excellent (as seen in the data above) but the corrected reads will violate some of the assumptions that assembers make about their input data. Some tools are aware of the difference between raw and corrected input reads, for example with `Canu` there is an optional parameter `-corrected` which tells the assembler that the input reads have already been corrected by an unspecified technology and they are treated differently to raw sequences.

### Hybrid assembly

Alternatively, it is possible for many assemblers to take a multitude of data types as input files and combine the data during assembly. This approach has been used to great effect in several publications (for example [Wick *et al.*, 2017](https://doi.org/10.1099/mgen.0.000132), [Chen *et al.*, 2020](https://doi.org/10.1186/s12864-020-07041-8), [Brown *et al.*, 2021](https://doi.org/10.1038/s41598-021-83081-8) [Neal-McKinney *et al.*, 2021](https://doi.org/10.1038/s41598-021-84956-6). Hybrid assembly is supported by a number of assemblers, including;

1. `SPAdes` ([Antipov *et al.*, 2016](https://doi.org/10.1093/bioinformatics/btv688))
1. `MaSuRCA` ([Zimin *et al.*, 2013](https://doi.org/10.1093/bioinformatics/btt476))
1. `Unicycler` / `Trycycler` ([Wick *et al.*, 2017](https://doi.org/10.1371/journal.pcbi.1005595), [Wick *et al.*, preprint](https://doi.org/10.1101/2021.07.04.451066))

>**Note:**The `Unicycler` assembler does make use of `SPAdes` under the hood, so many not be a complete separate approach from `hybridSPAdes`.

---

