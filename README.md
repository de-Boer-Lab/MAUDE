# MAUDE: Mean Alterations Using Discrete Expression

MAUDE is an R package for finding differences in means of normally distributed (or nearly so) data, via measuring abundances in discrete bins. For example, a pooled CRISPRi screen with expression readout by FACS sorting into discrete bins and sequencing the abundances of the guides in each bin. The documentation and examples are written with this usage in mind, but there is no reason why it can't be used for any experiment where normally distributed expression values are read out via abundances in discrete expression bins. See 'Usage' below for more information.


<img src="images/logo2.png" alt="Maude Flanders" width="400"/>

# Table of contents
<!--ts-->
   * [R Installation](#r-installation)
   * [Requirements](#requirements)
   * [Usage](#usage)
   * [Citation](#citation)
<!--te-->

# R Installation

If you don't already have `devtools`, install it:
```
install.packages("devtools")
```

Load `devtools` and install from the GitHub page:

```
library(devtools)
install_github("Carldeboer/MAUDE")
```
# Requirements
Right now we have three main requirements: 
1. Non-targeting guides are included in the experiment to serve as negative controls and are used for calibrating Z-scores and P-values.
2. The abundance of the guides must have been measured somehow (usually by sequencing the guide DNA of unsorted cells)
3. The fractions of cells sorted into each expression bin was quantified


# Usage

## Tutorials
We provide two tutorials on how to run a MAUDE analysis in R here:
1. [Re-analysis of CD69 screen data](Tutorial/Tutorial.md)
2. [Analysis of a simulated screen](Simulation.md)

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
> Carl G de Boer*, John Ray*, Nir Hacohen, and Aviv Regev, MAUDE, bioRxiv. https://www.biorxiv.org/content/10.1101/819649v1
