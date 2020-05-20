# Tutorial: Re-analysis of CD69 CRISPRa tiling screen

See Simeonov et al. Discovery of stimulation-responsive immune enhancers with CRISPR activation. Nature. 2017 Sep 7;549(7670):111-115. doi: 10.1038/nature23875. Epub 2017 Aug 30. [PMC5675716](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5675716/) for the original publication associated with this dataset.

This tutorial is designed to get you quickly up and running with MAUDE. The only requirements are that you have R and the following R libraries installed:
`ggplot2, openxlsx, MAUDE, reshape`

## Load required libraries into R
```R
#load required libraries
library(openxlsx)
library(reshape)
library(ggplot2)
library(MAUDE)
library(GenomicRanges)
library(ggbio)
library(Homo.sapiens)
```

## Load CD69 screen count data from Simeonov Supplementary Table 1
```R
#read in the CD69 screen data
CD69Data = read.xlsx('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5675716/bin/NIHMS913084-supplement-supplementary_table_1.xlsx')

#identify non-targeting guides
CD69Data$isNontargeting = grepl("negative_control", CD69Data$gRNA_systematic_name)

CD69Data = unique(CD69Data) # for some reason there were duplicated rows in this table - remove duplicates

#reshape the count data so we can label the experimental replicates and bins, and remove all the non-count data
cd69CountData = melt(CD69Data, id.vars = c("PAM_3primeEnd_coord","isNontargeting","gRNA_systematic_name"))
cd69CountData = cd69CountData[grepl(".count$",cd69CountData$variable),] # keep only read count columns
cd69CountData$Bin = gsub("CD69(.*)([12]).count","\\1",cd69CountData$variable)
cd69CountData$expt = gsub("CD69(.*)([12]).count","\\2",cd69CountData$variable)
cd69CountData$reads= as.numeric(cd69CountData$value); cd69CountData$value=NULL;
cd69CountData$Bin = gsub("_","",cd69CountData$Bin) # remove extra underscores

#reshape into a matrix
binReadMat = data.frame(cast(cd69CountData[!is.na(cd69CountData$PAM_3primeEnd_coord) | cd69CountData$isNontargeting,], PAM_3primeEnd_coord+gRNA_systematic_name+isNontargeting+expt ~ Bin, value="reads"))
#binReadMat now contains a matrix in the proper format for MAUDE analysis

```
##read in the Encode DHS peaks at the locus
We can use this as a reference point with which to view our results, and to annotate guides with regulatory elements. These are all the DHS peaks at CD69 in Jurkat cells, merged across two replicates.

```R
dhsPeakBED = read.table('https://raw.githubusercontent.com/Carldeboer/MAUDE/master/Tutorial/Encode_Jurkat_DHS_both.merged.bed', stringsAsFactors=F, row.names=NULL, sep="\t", header=F)
names(dhsPeakBED) = c("chrom","start","end");
#add a column to include peak names
dhsPeakBED$name = paste(dhsPeakBED$chrom, paste(dhsPeakBED$start, dhsPeakBED$end, sep="-"), sep=":")
```

## Read in and inspect bin fractions

```R
#read in the bin fractions derived from Simeonov et al Extended Data Fig 1a and the "digitize" R package
#Ideally, you derive this from the FACS sort data. 
binStats = read.table('https://raw.githubusercontent.com/Carldeboer/MAUDE/master/Tutorial/CD69_bin_percentiles.txt', stringsAsFactors=F, row.names=NULL, sep="\t", header=T)
binStats$fraction = binStats$binEndQ - binStats$binStartQ; #the fraction of cells captured is the difference in bin start and end percentiles

#plot the bins as the percentiles of the distribution captured by each bin
p = ggplot(binStats, aes(colour=Bin)) + geom_segment(aes(x=binStartQ, xend=binEndQ, y=fraction, yend=fraction)) + xlab("Bin bounds as percentiles") + ylab("Fraction of the distribution captured") +theme_classic()+scale_y_continuous(expand=c(0,0))+coord_cartesian(ylim=c(0,0.7)); print(p)
```
![Bin percentiles for the CD69 screen](CD69_bin_percentiles.png "Bin percentiles")
This shows the sizes of the bins in terms of the fractions of the overall distribution covered by each bin. In terms of expression space, we need to convert these to Z scores with the reverse normal CDF function `qnorm`
```R
#convert bin fractions to Z scores
binStats$binStartZ = qnorm(binStats$binStartQ)
binStats$binEndZ = qnorm(binStats$binEndQ)
```
Now let's plot the bins as the Z score bounds of a normal distribution, with an actual normal distribution overlaid in gray.
```R
p = ggplot(binStats, aes(colour=Bin))  + geom_density(data=data.frame(x=rnorm(100000)), aes(x=x), fill="gray", colour=NA)+ geom_segment(aes(x=binStartZ, xend=binEndZ, y=fraction, yend=fraction)) + xlab("Bin bounds as expression Z-scores") + ylab("Fraction of the distribution captured") +theme_classic()+scale_y_continuous(expand=c(0,0))+coord_cartesian(ylim=c(0,0.7)); print(p)
```
![Bin expression bounds for the CD69 screen](CD69_bin_expression_bounds.png "Bin bounds")

Here, we are forced to use the same distribution for both experimental replicates since we didn't have the underlying data. You should derive this data from the FACS sorting data separately for each replicate (although it doesn't affect things that much if the bins were drawn similarly).
```R
binStats = rbind(binStats, binStats) #duplicate data
binStats$expt = c(rep("1",4),rep("2",4)); #name the first duplicate expt "1" and the next expt "2";
```

## 1) MAUDE: Calculate guide level statistics
Now we've finally gotten to the part where we're running MAUDE.  To compute guide-level stats, we run `findGuideHitsAllScreens`
This step takes about a minute to run.
```R
guideLevelStats = findGuideHitsAllScreens(experiments = unique(binReadMat["expt"]), countDataFrame = binReadMat, binStats = binStats, sortBins = c("baseline","high","low","medium"), unsortedBin = "back", negativeControl = "isNontargeting")
```
Here, we calculated the mean expression of each guide
```R
# Plot the guide-level mus
p = ggplot(guideLevelStats, aes(x=mean, colour=isNontargeting, linetype=expt)) + geom_density()+theme_classic()+scale_y_continuous(expand=c(0,0)) + geom_vline(xintercept = 0)+xlab("Learned mean guide expression"); print(p);
```
![Guide mean distributions](CD69_guide_mu_dist.png "Guide mean distributions")
But notice how the distributions for the non-targeting guides are off-centre.  This is why we also calculate a re-centred Z-score using the non-targeting guides.

```R
# Plot the guide-level Zs
p = ggplot(guideLevelStats, aes(x=Z, colour=isNontargeting, linetype=expt)) + geom_density()+theme_classic()+scale_y_continuous(expand=c(0,0)) + geom_vline(xintercept = 0)+xlab("Learned guide expression Z score"); print(p)
```
![Guide Z-score distributions](CD69_guide_Z_dist.png "Guide Z distributions")
Now the non-targeting guides are centred at 0, as expected since they have no effect on expression.  

Now that we have guide-level statistics, we can inspect them. We can see there is a high correlation between replicates:
```R
guideEffectsByRep = cast(guideLevelStats, gRNA_systematic_name + isNontargeting +PAM_3primeEnd_coord ~ expt, value="Z")

p = ggplot(guideEffectsByRep[!guideEffectsByRep$isNontargeting,], aes(x=`1`, y=`2`)) + geom_point(size=0.3) + xlab("Replicate 1 Z score") + ylab("Replicate 2 Z score") + ggtitle(sprintf("r = %f",cor(guideEffectsByRep$`1`[!guideEffectsByRep$isNontargeting],guideEffectsByRep$`2`[!guideEffectsByRep$isNontargeting])))+theme_classic(); print(p)
```
![Guide Z correlations](CD69_guide_rep1_vs_rep2_Z.png "Guide Z correlations")
They are pretty highly correlated at the effect size.  Now let's map the results across the CD69 locus.

```R
dhsPos = min(guideLevelStats$Z)*1.05;
p=ggplot(guideLevelStats, aes(x=PAM_3primeEnd_coord, y=Z)) +geom_point(size=0.5)+facet_grid(expt ~.)+ geom_segment(data = dhsPeakBED, aes(x=start, xend=end,y=dhsPos, yend=dhsPos), colour="red") + theme_classic() + xlab("Genomic position") + ylab("Guide Z score"); print(p)
```
![Guide locus view](CD69_guide-level_locus_view.png "Guide locus view")
Here, the DHS peaks are shown in red. Clearly some regions contain many active guides, the majority of which have high Z-scores, indicating CD69 activation.

## 2) MAUDE: Identify active elements
There are two ways to get element-level statistics, depending on how the screen was done.  If you have regulatory element annotations, you can use `getElementwiseStats` to identify active elements.  If you did a tiling screen, you could use the same, but you can also identify active regions in an unbiased way with `getTilingElementwiseStats`.
### 2a) Get element level stats with sliding window
Now, we can combine adjacent guides in an unbiased way, using a sliding window across the locus to identify regions with more active guides than expected by chance. Here we use a sliding window of 200 bp, and the default minimum guide number (5). Any 200 bp with fewer than 5 guides will not be tested.

```R
guideLevelStats$chrom = "chr12"; # we need to tell it what chromosome our guides are on - they're all on chr12
slidingWindowElements = getTilingElementwiseStats(experiments = unique(binReadMat["expt"]), normNBSummaries = guideLevelStats, tails="both", window = 200, location = "PAM_3primeEnd_coord",chr="chrom",negativeControl = "isNontargeting")
names(slidingWindowElements)[names(slidingWindowElements)=="chr"]="chrom" #override the default chromosome field 'chr' with the GRanges compatible 'chrom'
```

Now that we have element-level stats, let's inspect them! First, let's look at the whole locus.
```R
dhsPos = min(slidingWindowElements$meanZ)*1.05;
p=ggplot(slidingWindowElements, aes(x=start, xend=end, y=meanZ,yend=meanZ, colour=FDR<0.01)) +geom_segment(size=1)+facet_grid(expt ~.) + theme_classic() + xlab("Genomic position") + ylab("Element Z score") + geom_hline(yintercept = 0) + geom_segment(data = dhsPeakBED, aes(x=start, xend=end,y=dhsPos, yend=dhsPos), colour="black"); print(p)
```
![Sliding element locus view](CD69_sliding_element-level_locus_view.png "Sliding element locus view")
In the above, we see that there are some regions significantly upregulated or down-regulated, which are often shared between the two replicates. The regions that upregulate CD69 preferentially lie in open chromatin (DHS - black, in the above).

How correlated are the effect size estimates for the two replicates at the element level?
```R
slidingWindowElementsByRep = cast(slidingWindowElements, chrom + start + end +numGuides ~ expt, value="meanZ")
p = ggplot(slidingWindowElementsByRep, aes(x=`1`, y=`2`)) + geom_point(size=0.5) + xlab("Replicate 1 element effect Z score") + ylab("Replicate 2 element effect Z score") + ggtitle(sprintf("r = %f",cor(slidingWindowElementsByRep$`1`,slidingWindowElementsByRep$`2`)))+theme_classic(); print(p)
```
![Sliding element Z correlations](CD69_sliding_element_rep1_vs_rep2_Z.png "Sliding element Z correlations")
Quite highly correlated!

### 2b) Get element level stats with annotated elements
Since we have annotated elements anyway in `dhsPeakBED`, it might help us to identify which guides should be combined.  This would be the only way to do it if we had only targeted these regions to begin with.

```R
#the next command annotates our guides with any DHS peak they lie in.
annotatedGuides = findOverlappingElements(guides = unique(guideLevelStats[!guideLevelStats$isNontargeting, c("PAM_3primeEnd_coord","gRNA_systematic_name","chrom")]), elements = dhsPeakBED, elements.start = "start", elements.end = "end", elements.chr = "chrom", guides.pos = "PAM_3primeEnd_coord", guides.chr = "chrom")

#merge regulatory element annotations back onto guideLevelStats
guideLevelStats = merge(guideLevelStats, annotatedGuides[c("gRNA_systematic_name", "name")], by="gRNA_systematic_name", all.x=T)

#this is where we are actually running MAUDE to find element-level stats
dhsPeakStats = getElementwiseStats(experiments = unique(binReadMat["expt"]), normNBSummaries = guideLevelStats, negativeControl = "isNontargeting", elementIDs = "name") # "name" is the peak IDs from the DHS BED file

#merge peak info back into dhsPeakStats
dhsPeakStats = merge(dhsPeakStats, dhsPeakBED, by="name");
```

Now we can again look at the activity of the elements across the locus:
```R
p=ggplot(dhsPeakStats, aes(x=start, xend=end, y=meanZ,yend=meanZ, colour=FDR<0.01)) +geom_segment(size=1)+facet_grid(expt ~.) + theme_classic() + xlab("Genomic position") + ylab("Element Z score") + geom_hline(yintercept = 0); print(p)
```
![DHS element locus view](CD69_DHS_element-level_locus_view.png "DHS element locus view")

And the correlation between replicates.
```R
dhsPeakStatsByRep = cast(dhsPeakStats, name ~ expt, value="meanZ")

p = ggplot(dhsPeakStatsByRep, aes(x=`1`, y=`2`)) + geom_point() + xlab("Replicate 1 DHS effect Z score") + ylab("Replicate 2 DHS effect Z score") + ggtitle(sprintf("r = %f",cor(dhsPeakStatsByRep$`1`,dhsPeakStatsByRep$`2`)))+theme_classic(); print(p)
```
![DHS element Z correlations](CD69_DHS_element_rep1_vs_rep2_Z.png "DHS element Z correlations")
Again, these are highly correlated.
In this case, it didn't look like the regulatory element annotations helped our analysis, but that is not always true. 

Let's see investicate why this might be by looking at the guide effect sizes for the guides within each DHS element.
```R
p=ggplot(guideLevelStats, aes(x=Z, group=name, colour=name == "chr12:9912678-9915275")) + stat_ecdf(alpha=0.3)+ 
  stat_ecdf(data=guideLevelStats[!is.na(guideLevelStats$name) & guideLevelStats$name=="chr12:9912678-9915275",], size=1)+
  facet_grid(expt ~.) + theme_classic() + xlab("Guide Z score")+scale_y_continuous(expand=c(0,0)) + 
  scale_x_continuous(expand=c(0,0)) + scale_colour_manual(values=c("black","red")) + labs(colour = "CD69 promoter?")+ylab("Cumulative fraction"); print(p)
```
![Guide Z score CDFs for elements](CD69_guide_Z_dist_per_element.png "Guide Z score CDFs for elements")
Here, the promoter DHS peak is in red and all the other DHS peaks are in black. You can appreciate that the promoter tends to have among the most influencial guides (rightward shift in the CDF curves). Even for the promoter, there is substantial variability in the estimated guide effect sizes, ranging from not doing anything (near 0) to pretty sizable effects (>3).  This is probably a combination of factors, including experimental noise, how effective each guide targets the region, and the effect of actually targeting each region.

Let's zoom in on the promoter to see what exactly is happening here.
```R
p=ggplot(guideLevelStats[!is.na(guideLevelStats$name) & guideLevelStats$name=="chr12:9912678-9915275",], aes(x=PAM_3primeEnd_coord, y=Z, colour=expt)) +
  geom_point(size=1)+ geom_line()+theme_classic() + xlab("Genomic position") + ylab("Guide Z score")+ geom_vline(xintercept = 9913497, colour="black"); print(p)
```
![Guide Z scores in the promoter](CD69_guide_Zs_in_promoter.png "Guide Z scores in the promoter")

Here, I've marked the transcription start site with a vertical black line. CD69 is an antisense gene, so it is transcribed to the left of the black line. The guides targeting the 5' UTR of CD69 (left of black line) appear to be less effective than those in the promoter (right of black line), which is expected. The high correlation in the guide effect sizes between the two experiments indicates that many of the low-effect size guides cannot be attributed to experimental noise. So probably most of the signal diversity within the promoter DHS is attributable to how effective the guide is at targeting the DNA, and how effective the activation domain of CRISPRa is where it is bound.

Now let's look more closely at the guides and select a few that work especially well.
```R
p=ggplot(guideEffectsByRep[guideEffectsByRep$PAM_3primeEnd_coord < 9915275 & guideEffectsByRep$PAM_3primeEnd_coord > 9912678,], aes(x=`2`, y=`1`))+geom_point()+geom_text(data=guideEffectsByRep[guideEffectsByRep$PAM_3primeEnd_coord < 9915275 & guideEffectsByRep$PAM_3primeEnd_coord > 9912678 & guideEffectsByRep$`1`>2 & guideEffectsByRep$`2`>2,], aes(label=gRNA_systematic_name))+theme_classic()+xlab("Replicate 2 guide Z score") + ylab("Replicate 1 guide Z score"); print(p) 

```
![Best guides in promoter](CD69_best_guides.png "Best promoter guides")

We can see that, although there are four guides that have a large effect in one replicate, they don't have as great an effect in the other. Here, I've labeled the IDs of three guides that have a Z score over two in both replicates. If I were going to follow up targeting the CD69 promoter with more guides, these are the ones I would use.

`GenomicRanges` is a great package to work with data associated with genomic coordinates.  First, let's restructure our tiling element data.frame so that there is only one row per region and there is one of each data field for each replicate. 

```R
slidingWindowElementsByReplicate = cast(melt(slidingWindowElements, id.vars=c("expt","numGuides","chrom","start","end")), numGuides+chrom+start+end ~ variable+expt, value="value")
head(slidingWindowElementsByReplicate)
```
```
  numGuides chrom   start     end  stoufferZ_1 stoufferZ_2      meanZ_1     meanZ_2 significanceZ_1 significanceZ_2  p.value_1
1         5 chr12 9882993 9883182 -0.009472366 -0.23953341 -0.004236171 -0.10712260     -0.06237373      -1.1481439 0.47513261
2         5 chr12 9888047 9888212 -0.213765598 -0.35832774 -0.095598882 -0.16024904     -1.40760595      -1.7175550 0.07962389
3         5 chr12 9893958 9894158  0.147176396  0.13554497  0.065819285  0.06061755      0.96912867       0.6497012 0.16624050
4         5 chr12 9898697 9898894  0.112898517 -0.03747552  0.050489752 -0.01675956      0.74341533      -0.1796296 0.22861513
5         5 chr12 9898702 9898898  0.135955105 -0.03142427  0.060800971 -0.01405336      0.89523860      -0.1506244 0.18532978
6         5 chr12 9901267 9901451 -0.021920411 -0.08749781 -0.009803106 -0.03913021     -0.14434175      -0.4193990 0.44261531
  p.value_2     FDR_1     FDR_2
1 0.1254546 0.9771498 0.4344497
2 0.0429389 0.3389303 0.2207419
3 0.2579426 0.5546233 0.6732716
4 0.4287217 0.6691416 0.9149289
5 0.4401360 0.5908959 0.9263483
6 0.3374623 0.9456464 0.7964048
```

We can convert our data structures to `GRanges` objects as easily as this:

```R
slidingWindowElementsByReplicateGR = GRanges(as.data.frame(slidingWindowElementsByReplicate)) #casting to data.frame is only needed if using cast
```

Now let's label those ranges that were significant in both replicates, and merge overlapping significant ranges to get a set of expanded active regions.

```R
#require that both replicates are significant at an FDR of 0.1 and that the signs agree
slidingWindowElementsByReplicateGR$significantUp   = slidingWindowElementsByReplicateGR$FDR_1< 0.01 & slidingWindowElementsByReplicateGR$FDR_2 < 0.01 & slidingWindowElementsByReplicateGR$meanZ_1 > 0 &  slidingWindowElementsByReplicateGR$meanZ_2 > 0; 
slidingWindowElementsByReplicateGR$significantDown = slidingWindowElementsByReplicateGR$FDR_1< 0.01 & slidingWindowElementsByReplicateGR$FDR_2 < 0.01 & slidingWindowElementsByReplicateGR$meanZ_1 < 0 &  slidingWindowElementsByReplicateGR$meanZ_2 < 0;

#merge overlapping regions in each set
overlappingSlidingWindowElementsUp =reduce(slidingWindowElementsByReplicateGR[slidingWindowElementsByReplicateGR$significantUp])
overlappingSlidingWindowElementsDown =reduce(slidingWindowElementsByReplicateGR[slidingWindowElementsByReplicateGR$significantDown])
```

Now let's plot it all together as a genome browser track!

```R

#which gene models do I want to plot?
data(genesymbol, package = "biovizBase")
wh <- genesymbol[c("CD69", "CLECL1", "KLRF1", "CLEC2D","CLEC2B")]
wh <- range(wh, ignore.strand = TRUE)

#make the genome tracks
tracks(autoplot(Homo.sapiens, which = wh, gap.geom="chevron"), autoplot(overlappingSlidingWindowElementsUp, fill="red"), autoplot(overlappingSlidingWindowElementsDown, fill="blue"), heights=c(5,2,2)) + theme_classic()
```

![Genome track](CD69_ggbio_locus.png "Genome track")

Here, the red are CD69-activating regions, blue are CD69-repressing regions.

