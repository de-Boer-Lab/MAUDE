# MAUDE: Mean Alterations Using Discrete Expression

[![DOI](https://zenodo.org/badge/135627989.svg)](https://zenodo.org/badge/latestdoi/135627989)

MAUDE is an R package for finding differences in means of normally distributed (or nearly so) data, via measuring abundances in discrete bins. For example, a pooled CRISPRi screen with expression readout by FACS sorting into discrete bins and sequencing the abundances of the guides in each bin.  Most of the documentation and examples are written with a CRISPRi-type sorting screen in mind, but there is no reason why it can't be used for any experiment where normally distributed expression values are read out via abundances in discrete expression bins. For example, MAUDE can also be used for [CRISPR base editor screens](https://de-boer-lab.github.io/MAUDE/doc/BACH2_base_editor_screen.html) where the readout is expression of a target gene, and reporter assays with expression readouts (e.g. [Rafi et al](https://www.biorxiv.org/content/10.1101/2023.04.26.538471v2)). 

See 'Usage' below for more information.


<img src="images/logo2.png" alt="Maude Flanders" width="400"/>

# Table of contents
<!--ts-->
   * [R Installation](#r-installation)
   * [Requirements](#requirements)
   * [Usage](#usage)
   * [Citation](#citation)
<!--te-->

# R Installation

## Option 1: Install directly from GitHub

If you don't already have `devtools`, install it:
```
install.packages("devtools")
```

Load `devtools` and install from the GitHub page:

```
devtools::install_github("de-Boer-Lab/MAUDE")
```

## Option 2: Install from download

Download the latest MAUDE release (Under "Releases" on the right hand side of this page).

Then in R, run:
```
install.packages("C:\\Path\\To\\Download\\MAUDE-1.0.1.zip")
```

# Requirements
Right now we have three main requirements: 
1. Negative control guides are included in the experiment; (these are used for calibrating Z-scores and P-values, and so are not strictly needed if only the expression means are desired).
2. The abundance of the guides must have been measured somehow (usually by sequencing the guide DNA of unsorted cells; though there are ways to estimate this post-sort if the bins cover the majority of the distribution)
3. The fractions of cells sorted into each expression bin was quantified (typically the cell counts/fractions read off of the cell sorter)


# Usage

## Tutorials
We provide two tutorials on how to run a MAUDE analysis in R here:
1. [Re-analysis of CD69 screen data](https://de-boer-lab.github.io/MAUDE/doc/CD69_tutorial.html)
2. [Analysis of a simulated screen](https://de-boer-lab.github.io/MAUDE/doc/simulated_data_tutorial.html)
3. [Analysis of a CRISPR base editor non-coding mutation screen](https://de-boer-lab.github.io/MAUDE/doc/BACH2_base_editor_screen.html)

For additional examples, see the [script for evaluating and comparing sorting-based CRISPR screen analysis methods.](https://de-boer-lab.github.io/MAUDE/Evaluation/method_evaluation.html)

## Quantifying guide DNA abundance
After sequencing, you get fastqs, one per sorting bin and experiment.  The first step for a MAUDE analysis is to quantify the number of guides residing in each bin.  Here, we provide some guidance as to how to do this.

We have previously used the aligner [`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

To make the bowtie2 reference `guide_seq_reference`:
```bash
bowtie2-build guide_seqs.fa guide_seq_reference

```
where `guide_seqs.fa` is a fasta file including the sequences you are mapping against, which will include the guide DNA sequence and any flanking constant regions as well. The amount of constant sequence you include in the reference should be at least as much as what was sequenced.

For example, with 20bp guides with constant flanking `GTTTAAGAGCTATGCTGGAAACAGCATAG`:
```
>guide1
GTCGCATATCGCGATAGCGAGTTTAAGAGCTATGCTGGAAACAGCATAG
>guide2
GTCGTGAAAGTGCTGTTGAGGTTTAAGAGCTATGCTGGAAACAGCATAG
...
```

The following command is an example of how to quantify guide abundance into a format that can easily be input into R for MAUDE analysis:
```bash
bowtie2 --no-head -x guide_seq_reference -U $sample.fastq.gz -S $sample.mapped.sam
#here, we include all mapped reads, but by using Samtools, you can filter out reads that map to the wrong strand, have indels, etc.
cat $sample.mapped.sam | awk '{print $3}' | sort | uniq -c | sort > $sample.counts
```
Here, `$sample` is the sample name, with `$sample.fastq.gz` the corresponding fastq file, and `guide_seq_reference` is the `bowtie2` reference.  The file `$sample.counts` will contain guide counts that can be input into R. 

To turn this into a format that can easily be used for a MAUDE analysis, you can input the data using something like the following:
```R
#here, allSamples is a data.frame containing one sample per row, with columns including ID, expt, and Bin.  There should be one file for every row in allSamples
allData = data.frame();
for (i in 1:nrow(allSamples)){
  curData = read.table(file=sprintf("%s/%s.counts",inDir,allSamples$ID[i]), quote="", header = F, row.names = NULL, stringsAsFactors = F)
  names(curData) = c("count","guideID");
  curData = curData[curData$gID!="*",] # remove unmapped counts
  curData$ID = allSamples$ID[i];
  curData$expt = allSamples$expt[i];
  curData$Bin = allSamples$Bin[i];
  allData = rbind(allData, curData)
}
#now you have the data in a data.frame that can be reshaped to a MAUDE-compatible format:
library(reshape)
allDataCounts = as.data.frame(cast(allData, expt + guideID ~ Bin, value="count"));
allDataCounts[is.na(allDataCounts)]=0; # fill in 0s for guides not observed at all
#now you just need to label the non-targeting guides and this will be in the correct format
```

## Encountering problems
Should you encounter a problem using MAUDE:
1. [Consult the Common Problems](CommonProblems.md)
2. [Submit an Issue](https://github.com/Carldeboer/MAUDE/issues)
3. Contact the authors.


# Citation
Please cite:

Carl G de Boer*, John P Ray*, Nir Hacohen, Aviv Regev. [_MAUDE: Inferring Expression Changes in Sorting-Based CRISPR Screens_.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02046-8) 2020 Jun 3;21(1):134. doi: 10.1186/s13059-020-02046-8. PMID: [32493396](https://pubmed.ncbi.nlm.nih.gov/32493396/).
