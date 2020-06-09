## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load libraries, results="hide", message = FALSE, warning=FALSE-----------
#load required libraries
library(openxlsx)
library(reshape)
library(ggplot2)
library(MAUDE)
library(GenomicRanges)
library(ggbio)
library(Homo.sapiens)

## ----input data---------------------------------------------------------------
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
binReadMat = data.frame(cast(cd69CountData[!is.na(cd69CountData$PAM_3primeEnd_coord) | cd69CountData$isNontargeting,], 
  PAM_3primeEnd_coord+gRNA_systematic_name+isNontargeting+expt ~ Bin, value="reads"))
#binReadMat now contains a matrix in the proper format for MAUDE analysis


## ----input set of DHS peaks---------------------------------------------------
dhsPeakBED = read.table(system.file("extdata", "Encode_Jurkat_DHS_both.merged.bed", package = "MAUDE", mustWork = TRUE), 
  stringsAsFactors=FALSE, row.names=NULL, sep="\t", header=FALSE)
names(dhsPeakBED) = c("chrom","start","end");
#add a column to include peak names
dhsPeakBED$name = paste(dhsPeakBED$chrom, paste(dhsPeakBED$start, dhsPeakBED$end, sep="-"), sep=":")

## ----read in bin fractions----------------------------------------------------
#read in the bin fractions derived from Simeonov et al Extended Data Fig 1a and the "digitize" R package
#Ideally, you derive this from the FACS sort data. 
binStats = read.table(system.file("extdata", "CD69_bin_percentiles.txt", package = "MAUDE", mustWork = TRUE), 
  stringsAsFactors=FALSE, row.names=NULL, sep="\t", header=TRUE)
binStats$fraction = binStats$binEndQ - binStats$binStartQ; #the fraction of cells captured is the difference in bin start and end percentiles

#plot the bins as the percentiles of the distribution captured by each bin
p = ggplot(binStats, aes(colour=Bin)) + ggplot2::geom_segment(aes(x=binStartQ, xend=binEndQ, y=fraction, yend=fraction)) + 
  xlab("Bin bounds as percentiles") + ylab("Fraction of the distribution captured") +theme_classic() + 
  scale_y_continuous(expand=c(0,0))+coord_cartesian(ylim=c(0,0.7)); print(p)

## ----convert bin percentiles to Z scores--------------------------------------
#convert bin fractions to Z scores
binStats$binStartZ = qnorm(binStats$binStartQ)
binStats$binEndZ = qnorm(binStats$binEndQ)

## ----plot bins----------------------------------------------------------------
p = ggplot(binStats, aes(colour=Bin))  + 
  geom_density(data=data.frame(x=rnorm(100000)), aes(x=x), fill="gray", colour=NA)+ 
  ggplot2::geom_segment(aes(x=binStartZ, xend=binEndZ, y=fraction, yend=fraction)) + 
  xlab("Bin bounds as expression Z-scores") + 
  ylab("Fraction of the distribution captured") +theme_classic()+scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0,0.7)); print(p)

## ----duplicate bins for second replicate--------------------------------------
binStats = rbind(binStats, binStats) #duplicate data
binStats$expt = c(rep("1",4),rep("2",4)); #name the first duplicate expt "1" and the next expt "2";

## ----find guide effects-------------------------------------------------------
guideLevelStats = findGuideHitsAllScreens(experiments = unique(binReadMat["expt"]), 
  countDataFrame = binReadMat, binStats = binStats, 
  sortBins = c("baseline","high","low","medium"), 
  unsortedBin = "back", negativeControl = "isNontargeting")

## ----plot guide effects-------------------------------------------------------
# Plot the guide-level mus
p = ggplot(guideLevelStats, aes(x=mean, colour=isNontargeting, linetype=expt)) + geom_density()+
  theme_classic()+scale_y_continuous(expand=c(0,0)) + geom_vline(xintercept = 0)+
  xlab("Learned mean guide expression"); print(p);

## ----plot guide level Zs------------------------------------------------------
# Plot the guide-level Zs
p = ggplot(guideLevelStats, aes(x=Z, colour=isNontargeting, linetype=expt)) + geom_density()+
  theme_classic()+scale_y_continuous(expand=c(0,0)) + geom_vline(xintercept = 0)+
  xlab("Learned guide expression Z score"); 
print(p)

## ----plot replicate scatter---------------------------------------------------
guideEffectsByRep = cast(guideLevelStats, 
  gRNA_systematic_name + isNontargeting + PAM_3primeEnd_coord ~ expt, value="Z")

p = ggplot(guideEffectsByRep[!guideEffectsByRep$isNontargeting,], aes(x=`1`, y=`2`)) + 
  geom_point(size=0.3) + xlab("Replicate 1 Z score") + ylab("Replicate 2 Z score") + 
  ggtitle(sprintf("r = %f",cor(guideEffectsByRep$`1`[!guideEffectsByRep$isNontargeting],
    guideEffectsByRep$`2`[!guideEffectsByRep$isNontargeting])))+theme_classic(); 
print(p)

## ----plot locus---------------------------------------------------------------
dhsPos = min(guideLevelStats$Z)*1.05;
p=ggplot(guideLevelStats, aes(x=PAM_3primeEnd_coord, y=Z)) +geom_point(size=0.5)+facet_grid(expt ~.)+ 
  ggplot2::geom_segment(data = dhsPeakBED, aes(x=start, xend=end,y=dhsPos, yend=dhsPos), colour="red") + 
  theme_classic() + xlab("Genomic position") + ylab("Guide Z score"); 
print(p)

## ----infer sliding window effects---------------------------------------------
guideLevelStats$chrom = "chr12"; # we need to tell it what chromosome our guides are on - they're all on chr12
slidingWindowElements = getTilingElementwiseStats(experiments = unique(binReadMat["expt"]), 
  normNBSummaries = guideLevelStats, tails="both", window = 200, location = "PAM_3primeEnd_coord",
  chr="chrom",negativeControl = "isNontargeting")
#override the default chromosome field 'chr' with the GRanges compatible 'chrom'
names(slidingWindowElements)[names(slidingWindowElements)=="chr"]="chrom" 

## ----tiles locus effects------------------------------------------------------
dhsPos = min(slidingWindowElements$meanZ)*1.05;
p=ggplot(slidingWindowElements, aes(x=start, xend=end, y=meanZ,yend=meanZ, colour=FDR<0.01)) +
  ggplot2::geom_segment(size=1)+facet_grid(expt ~.) + theme_classic() + xlab("Genomic position") + 
  ylab("Element Z score") + geom_hline(yintercept = 0) + 
  ggplot2::geom_segment(data = dhsPeakBED, aes(x=start, xend=end,y=dhsPos, yend=dhsPos), colour="black");
print(p)

## ----tiled replicate scatter--------------------------------------------------
slidingWindowElementsByRep = cast(slidingWindowElements, chrom + start + end +numGuides ~ expt, 
  value="meanZ")
p = ggplot(slidingWindowElementsByRep, aes(x=`1`, y=`2`)) + geom_point(size=0.5) + 
  xlab("Replicate 1 element effect Z score") + ylab("Replicate 2 element effect Z score") + 
  ggtitle(sprintf("r = %f",cor(slidingWindowElementsByRep$`1`,slidingWindowElementsByRep$`2`)))+
  theme_classic(); 
print(p)

## ----element-level stats------------------------------------------------------
#the next command annotates our guides with any DHS peak they lie in.
annotatedGuides = findOverlappingElements(guides = unique(guideLevelStats[!guideLevelStats$isNontargeting,
  c("PAM_3primeEnd_coord","gRNA_systematic_name","chrom")]), elements = dhsPeakBED, 
  elements.start = "start", elements.end = "end", elements.chr = "chrom", 
  guides.pos = "PAM_3primeEnd_coord", guides.chr = "chrom")

#merge regulatory element annotations back onto guideLevelStats
guideLevelStats = merge(guideLevelStats, annotatedGuides[c("gRNA_systematic_name", "name")],
                        by="gRNA_systematic_name", all.x=TRUE)

#this is where we are actually running MAUDE to find element-level stats
dhsPeakStats = getElementwiseStats(experiments = unique(binReadMat["expt"]), 
  normNBSummaries = guideLevelStats, negativeControl = "isNontargeting", 
  elementIDs = "name") # "name" is the peak IDs from the DHS BED file

#merge peak info back into dhsPeakStats
dhsPeakStats = merge(dhsPeakStats, dhsPeakBED, by="name");

## ----element locus effect view------------------------------------------------
p=ggplot(dhsPeakStats, aes(x=start, xend=end, y=meanZ,yend=meanZ, colour=FDR<0.01)) +
  ggplot2::geom_segment(size=1)+facet_grid(expt ~.) + theme_classic() + xlab("Genomic position") + 
  ylab("Element Z score") + geom_hline(yintercept = 0); 
print(p)

## ----element replicate scatter------------------------------------------------
dhsPeakStatsByRep = cast(dhsPeakStats, name ~ expt, value="meanZ")

p = ggplot(dhsPeakStatsByRep, aes(x=`1`, y=`2`)) + geom_point() + 
  xlab("Replicate 1 DHS effect Z score") + ylab("Replicate 2 DHS effect Z score") + 
  ggtitle(sprintf("r = %f",cor(dhsPeakStatsByRep$`1`,dhsPeakStatsByRep$`2`)))+theme_classic(); 
print(p)

## ----guide effects per element------------------------------------------------
p=ggplot(guideLevelStats, aes(x=Z, group=name, colour=name == "chr12:9912678-9915275")) + stat_ecdf(alpha=0.3)+ 
  stat_ecdf(data=guideLevelStats[!is.na(guideLevelStats$name) &
                                   guideLevelStats$name=="chr12:9912678-9915275",], size=1)+
  facet_grid(expt ~.) + theme_classic() + xlab("Guide Z score")+scale_y_continuous(expand=c(0,0)) + 
  scale_x_continuous(expand=c(0,0)) + scale_colour_manual(values=c("black","red")) + 
  labs(colour = "CD69 promoter?")+ylab("Cumulative fraction"); 
print(p)

## ----promoter view------------------------------------------------------------
p=ggplot(guideLevelStats[!is.na(guideLevelStats$name) & guideLevelStats$name=="chr12:9912678-9915275",],
         aes(x=PAM_3primeEnd_coord, y=Z, colour=expt)) +
  geom_point(size=1)+ geom_line()+theme_classic() + xlab("Genomic position") + ylab("Guide Z score")+
  geom_vline(xintercept = 9913497, colour="black"); 
print(p)

## ----promoter guide zoom------------------------------------------------------
p=ggplot(guideEffectsByRep[guideEffectsByRep$PAM_3primeEnd_coord < 9915275 & 
                             guideEffectsByRep$PAM_3primeEnd_coord > 9912678,], 
         aes(x=`2`, y=`1`)) + geom_point() + 
  geom_text(data=guideEffectsByRep[guideEffectsByRep$PAM_3primeEnd_coord < 9915275 & 
                                     guideEffectsByRep$PAM_3primeEnd_coord > 9912678 & 
                                     guideEffectsByRep$`1`>2 & guideEffectsByRep$`2`>2,],
            aes(label=gRNA_systematic_name)) + 
  theme_classic()+xlab("Replicate 2 guide Z score") + ylab("Replicate 1 guide Z score"); 
print(p) 


## ----reshaping window effects-------------------------------------------------
slidingWindowElementsByReplicate = cast(melt(slidingWindowElements, 
  id.vars=c("expt","numGuides","chrom","start","end")), 
  numGuides+chrom+start+end ~ variable+expt, value="value")
head(slidingWindowElementsByReplicate)

## ----cast to GRanges----------------------------------------------------------
#casting to data.frame is only needed if using cast
slidingWindowElementsByReplicateGR = GRanges(as.data.frame(slidingWindowElementsByReplicate)) 

## ----find doubly significant tiles--------------------------------------------
#require that both replicates are significant at an FDR of 0.1 and that the signs agree
slidingWindowElementsByReplicateGR$significantUp   = slidingWindowElementsByReplicateGR$FDR_1< 0.01 & 
  slidingWindowElementsByReplicateGR$FDR_2 < 0.01 & slidingWindowElementsByReplicateGR$meanZ_1 > 0 &
  slidingWindowElementsByReplicateGR$meanZ_2 > 0; 
slidingWindowElementsByReplicateGR$significantDown = slidingWindowElementsByReplicateGR$FDR_1< 0.01 &
  slidingWindowElementsByReplicateGR$FDR_2 < 0.01 & slidingWindowElementsByReplicateGR$meanZ_1 < 0 &
  slidingWindowElementsByReplicateGR$meanZ_2 < 0;

#merge overlapping regions in each set
overlappingSlidingWindowElementsUp =
  reduce(slidingWindowElementsByReplicateGR[slidingWindowElementsByReplicateGR$significantUp])
overlappingSlidingWindowElementsDown =
  reduce(slidingWindowElementsByReplicateGR[slidingWindowElementsByReplicateGR$significantDown])

## ----genome browser view, fig.width=10, fig.height=5--------------------------

#which gene models do I want to plot?
data(genesymbol, package = "biovizBase")
wh <- genesymbol[c("CD69", "CLECL1", "KLRF1", "CLEC2D","CLEC2B")]
wh <- range(wh, ignore.strand = TRUE)

#make the genome tracks
tracks(autoplot(Homo.sapiens, which = wh, gap.geom="chevron"), 
  autoplot(overlappingSlidingWindowElementsUp, fill="red"), 
  autoplot(overlappingSlidingWindowElementsDown, fill="blue"), heights=c(5,2,2)) + theme_classic()

