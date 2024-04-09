



#' Find overlaps between guides and annotated elements
#'
#' Finds guides that overlap the elements of a BED-like data.frame (e.g. open chromatin regions) and returns a new data.frame containing those overlaps
#' @param guides a data.frame containing guide information including the guides genomic position in a column named guides.pos
#' @param elements a data.frame containing element information, as in a BED file, including the element's genomic start, end, and chromosome in columns named elements.start, elements.end, and elements.chr
#' @param guides.pos the name of the column in guides that contains the genomic position targeted by the guide (defaults to "pos")
#' @param guides.chr the name of the column in guides that contains the genomic chromosome targeted by the guide  (defaults to "chr")
#' @param elements.start the name of the column in elements that contains the start coordinate of the element (defaults to "st")
#' @param elements.end the name of the column in elements that contains the start coordinate of the element (defaults to "en")
#' @param elements.chr the name of the column in elements that contains the start coordinate of the element (defaults to "chr")
#' @return Returns a new data.frame containing the intersection of elements and guides
#' @export
#' @examples
#' set1 = data.frame(gid=1:10, chr=c(rep("chr1",5), rep("chr5",5)),
#'   pos= c(1:5,1:5)*10, stringsAsFactors = FALSE)
#' set2 = data.frame(eid=1:4, chr=c("chr1","chr1","chr4","chr5"), st=c(5,25,1,45),
#'     en=c(15,50,50,55), stringsAsFactors = FALSE)
#' findOverlappingElements(set1, set2)
findOverlappingElements = function(guides,elements,guides.pos="pos", guides.chr="chr",elements.start="st",elements.end="en", elements.chr="chr"){
  if (! all(c(guides.pos, guides.chr) %in% names(guides))){
    stop(sprintf("Cannot find all of %s and %s in guides data.frame!  Did you forget to set guides.pos or guides.chr?", guides.pos, guides.chr))
  }
  if (! all(c(elements.start, elements.end, elements.chr) %in% names(elements))){
    stop(sprintf("Cannot find all of %s, %s, and %s in elements data.frame!  Did you forget to set elements.start or elements.end?", elements.start, elements.end, elements.chr))
  }
  if(any(is.na(c(elements[[elements.start]], elements[[elements.end]], elements[[elements.chr]])))){
    stop(sprintf("Elements data.frame contains NAs in either %s, %s, or %s.",elements.start,elements.end,elements.chr))
  }
  if(any(is.na(c(guides[[guides.pos]], guides[[guides.chr]])))){
    stop(sprintf("Guides data.frame contains NAs in either %s, or %s.",guides.pos, guides.chr))
  }
  intersection = data.frame()
  for (xi in 1:nrow(elements)){
    curGuides = guides[elements[[elements.start]][xi] <= guides[[guides.pos]] & elements[[elements.end]][xi] >= guides[[guides.pos]] & elements[[elements.chr]][xi]==guides[[guides.chr]],];
    if (nrow(curGuides)>0){
      curGuides[[guides.chr]]=NULL; #remove chr column 
      intersection = rbind(intersection, cbind(curGuides,elements[xi,]))
    }
  }
  return(intersection)
}



#' Combines Z-scores using Stouffer's method
#'
#' This function takes a vector of Z-scores and combines them into a single Z-score using Stouffer's method.
#' @param x a vector of Z-scores to be combined
#' @return Returns a single Z-score.
#' @export
#' @examples
#' combineZStouffer(rnorm(10))
combineZStouffer = function(x){sum(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))}

#' Calculate the log likelihood of observed read counts
#'
#' Uses a normal distribution (N(mu,sigma)) to estimate how many reads are expected per bin under nullModel, and calculates the log likelihood under a negative binomial model. This function is usually not used directly.
#' @param x a vector of guide counts per bin
#' @param mu the mean for the normal expression distribution
#' @param k the vector of total counts per bin
#' @param sigma for the normal expression distribution (defaults to 1)
#' @param nullModel the bin bounds for the null model (for no change in expression)
#' @param libFract the fraction of the unsorted library this guide comprises (e.g. from unsorted cells, or sequencing the vector)
#' @return the log likelihood
#' @export
#' @examples
#' #usually not used directly
#' #make a bin sorting model with 6 10% bins
#' curSortBins = makeBinModel(data.frame(Bin = c("A","B","C","D","E","F"), fraction = rep(0.1,6)))
#' readsForGuideX =c(10,20,30,100,200,100); #the reads for this guide
#' getNBGaussianLikelihood(x=readsForGuideX, mu=1, k=rep(1E6,6), sigma=1, nullModel=curSortBins, 
#'   libFract = 50/1E6)
#' getNBGaussianLikelihood(x=readsForGuideX, mu=-1, k=rep(1E6,6), sigma=1, nullModel=curSortBins, 
#'   libFract = 50/1E6)
#' #mu=1 is far more likely (closer to 0) than mu=-1 for this distribution of reads
#' @import stats
getNBGaussianLikelihood = function(x, mu, k, sigma=1, nullModel, libFract){
  #mu = mean of distribution
  #k = vector of total counts per bin
  #sigma = sigma for normal distribution
  #nullModel = bin bounds for null model
  #calculate the probabilities of being within each bin
  binFractions = pnorm(nullModel$binEndZ-mu) - pnorm(nullModel$binStartZ-mu)
  binFractions = binFractions*libFract/(nullModel$binEndQ - nullModel$binStartQ)
  likelihood =0;
  for (i in 1:nrow(nullModel)){
    #dnbinom(x = number of reads for this guide, size = number of reads total, prob= probability of getting a read at each drawing)
    suppressWarnings({likelihood = likelihood  + dnbinom(x=x[i], size=k[i], prob =1- binFractions[i], log=TRUE)})
  }
  return(likelihood)
}

#' Create a bin model for a single experiment
#'
#' Provided with the fractions captured by each bin, creates a bin model for use with MAUDE analysis, assuming 3 contiguous bins on the tails of the distribution. You can easily remove or rename bins after they have been created with this function. An example is provided in the BACH2 Vignette.
#' @param curBinBounds  a data.frame containing two columns: Bin (must be {A,B,C,D,E,F}), and fraction (the fractions of the total captured by each bin)
#' @param tailP the fraction of the tails of the distribution not captured in any bin (defaults to 0.001)
#' @return returns a data.frame with additional columns including the bin starts and ends in Z-score space, and in quantile space.
#' @export
#' @examples
#' #generally, the bin bounds are retrieved from the FACS data
#' binBounds = makeBinModel(data.frame(Bin=c("A","B","C","D","E","F"), fraction=rep(0.1,6))) 
#' if(require("ggplot2")){
#'   p = ggplot() + 
#'     geom_vline(xintercept= sort(unique(c(binBounds$binStartZ, binBounds$binEndZ))),colour="gray")+
#'     theme_classic() + xlab("Target expression") + 
#'     geom_segment(data=binBounds, aes(x=binStartZ, xend=binEndZ, colour=Bin, y=0, yend=0), 
#'       size=5, inherit.aes = FALSE); 
#'   print(p)
#' }
#' @import stats
makeBinModel = function(curBinBounds,tailP=0.001){
  if (!all(curBinBounds$Bin %in% c("A","B","C","D","E","F"))){
    stop("Bin IDs must be A-F inclusive")
  }
  if (sum(curBinBounds$Bin %in% c("A","B","C","D","E","F"))<2){
    stop("Must have at least two bins.")
  }
  curBinBounds = curBinBounds[c("Bin","fraction")] # only want these two columns;
  
  #add any missing bins (null bins with no area)
  if (!all(c("A","B","C","D","E","F") %in% curBinBounds$Bin)){
    curBinBounds = rbind(curBinBounds, data.frame(Bin = c("A","B","C","D","E","F")[! c("A","B","C","D","E","F") %in% curBinBounds$Bin], fraction=0.00000000001))
  }
  #make bin cumulative probabilities (Q) 
  lastP = tailP;
  for (b in c("A","B","C")){
    curBinBounds$binStartQ[curBinBounds$Bin==b] = lastP
    curBinBounds$binEndQ[curBinBounds$Bin==b] = curBinBounds$binStartQ[curBinBounds$Bin==b] + curBinBounds$fraction[curBinBounds$Bin==b];
    lastP= curBinBounds$binEndQ[curBinBounds$Bin==b];
  }
  lastP = 1-tailP;
  for (b in c("F","E","D")){
    curBinBounds$binEndQ[curBinBounds$Bin==b] = lastP
    curBinBounds$binStartQ[curBinBounds$Bin==b] = curBinBounds$binEndQ[curBinBounds$Bin==b] - curBinBounds$fraction[curBinBounds$Bin==b];
    lastP= curBinBounds$binStartQ[curBinBounds$Bin==b];
  }
  curBinBounds$binStartZ = qnorm(curBinBounds$binStartQ)
  curBinBounds$binEndZ = qnorm(curBinBounds$binEndQ)
  return(curBinBounds)
}

#' Calculate guide-level statistics for a single screen
#'
#' Given a table of counts per guide/bin and a bin model for an experiment, calculate the optimal mean expression for each guide
#' @param countTable a table containing one column for each bin (A-F) and another column for non-targeting guide (logical-"NT"), and unsorted abundance (NS)
#' @param curBinBounds a bin model as created by makeBinModel
#' @param pseudocount the count to be added to each bin count, per 1e6 reads/bin total (default=10 pseudo reads per 1e6 reads total)
#' @param meanFunction how to calculate the mean of the non-targeting guides for centering Z-scores.  Defaults to 'mean'
#' @param sortBins the names in countTable of the sorting bins.  Defaults to c("A","B","C","D","E","F")
#' @param unsortedBin the name in countTable of the unsorted bin.  Defaults to "NS"
#' @param negativeControl the name in countTable containing a logical representing whether or not the guide is non-Targeting (i.e. a negative control guide).  Defaults to "NT"
#' @param limits the limits to the mu optimization. Defaults to c(-4,4)
#' @return a data.frame containing the guide-level statistics, including the Z score 'Z', log likelihood ratio 'llRatio', and estimated mean expression 'mean'.
#' @export
#' @examples
#' curSortBins = makeBinModel(data.frame(Bin = c("A","B","C","D","E","F"), fraction = rep(0.1,6)))
#' fakeReadData = data.frame(id=1:1000, A=rpois(1000, lambda = 100), B=rpois(1000, lambda = 100),
#'                           C=rpois(1000, lambda = 100), D=rpois(1000, lambda = 100),
#'                           E=rpois(1000, lambda = 100), F=rpois(1000, lambda = 100),
#'                           NotSorted=rpois(1000, lambda = 100), negControl = rnorm(1000)>0)
#' guideHits = findGuideHits(fakeReadData, curSortBins, unsortedBin = "NotSorted", 
#'   negativeControl="negControl")
#' if(require("ggplot2")){
#'   p=ggplot(guideHits, aes(x=Z, colour=negControl))+geom_density(); print(p)
#' }
#' @import stats
findGuideHits = function(countTable, curBinBounds, pseudocount=10, meanFunction = mean, sortBins = c("A","B","C","D","E","F"), unsortedBin = "NS", negativeControl="NT", limits=c(-4,4)){
  if(nrow(curBinBounds) > nrow(unique(curBinBounds))){ stop("Duplicate rows in curBinBounds table!")}
  if(!negativeControl %in% names(countTable)){
    stop(sprintf("Column '%s' is missing from countTable; did you specify 'negativeControl='?", negativeControl))
  }
  if(!unsortedBin %in% names(countTable)){
    stop(sprintf("Column '%s' is missing from countTable; did you specify 'unsortedBin='?", unsortedBin))
  }
  if(!all(sortBins %in% names(countTable))){
    stop(sprintf("Not all sort bins '%s' are present as columns in countTable; did you specify 'sortBins='?", paste(sortBins, collapse=", ")))
  }
  if(!all(curBinBounds$Bin %in% sortBins)){
    message("Warning: curBinBounds contained bins that are not present in sortBins: removing these.")
    curBinBounds = curBinBounds[curBinBounds$Bin %in% sortBins,]
  }
  if(any((curBinBounds$binEndQ-curBinBounds$binStartQ)<=0)){
    stop("Some sort bins cover <=0% of the distribution")
  }
  #the next line resorts curBinBounds to be in the same bin order as sortBins; this is assumed later in the code because sortBins is used to extract the count matrix
  curBinBounds = curBinBounds[order(as.numeric(factor(as.character(curBinBounds$Bin), levels=sortBins))),]
  if (pseudocount>0){
    for(b in sortBins){
      countTable[b]=countTable[b]+max(1,round(pseudocount*(sum(countTable[b])/1E6)));#add a pseudocount in proportion to depth
    }
    #countTable[sortBins]=countTable[sortBins]+pseudocount;#add a pseudocount
    countTable[unsortedBin]=countTable[unsortedBin]+pseudocount;
  }
  curNormNBSummaries = countTable
  countTable$libFraction = countTable[[unsortedBin]]/sum(countTable[[unsortedBin]],na.rm=TRUE)
  
  curNormNBSummaries$libFraction = curNormNBSummaries[[unsortedBin]]/sum(curNormNBSummaries[[unsortedBin]],na.rm=TRUE)
  binCounts = apply(curNormNBSummaries[sortBins],2,function(x){sum(x, na.rm = TRUE)})
  
  #for each guide, find the optimal mu given the count data and bin percentages
  for (i in 1:nrow(curNormNBSummaries)){
    curOptim=list()
    #interval: The probability of observing a guide outside of this interval in one of the non-terminal bins is very unlikely, and so estimating a true mean for these is too difficult anyway. Besides, we get some local optima at the extremes for sparsely sampled data.
    tryCatch(
      {curOptim = optim(0, fn = function(mu) { -getNBGaussianLikelihood(x = as.numeric(curNormNBSummaries[sortBins][i, ]), mu = mu, k = binCounts, nullModel = curBinBounds, libFract = curNormNBSummaries$libFraction[i]) }, method = "L-BFGS-B", lower = limits[1], upper = limits[2])},
      error =function(e){
        message("Mu optimization returned NaN objective: restricting search space")
        uniformObjectiveEval = data.frame(mu = ((0:1000)/1000 - 1) * (limits[2] - limits[1]) + limits[2], ll = NA)
        for (k in 1:nrow(uniformObjectiveEval)) {
          uniformObjectiveEval$ll[k] = getNBGaussianLikelihood(x = as.numeric(curNormNBSummaries[sortBins][i, ]), mu = uniformObjectiveEval$mu[k], k = binCounts, nullModel = curBinBounds, libFract = curNormNBSummaries$libFraction[i])
        }
        uniformObjectiveEval = uniformObjectiveEval[!is.infinite(uniformObjectiveEval$ll), ]
        uniformObjectiveEval = uniformObjectiveEval[!is.nan(uniformObjectiveEval$ll), ]
        stopifnot(nrow(uniformObjectiveEval) >= 1) # no points in the tested space have valid LLs
        message(sprintf("New interval = c(%f, %f)", min(uniformObjectiveEval$mu), max(uniformObjectiveEval$mu)))
        curOptim <<- optim(0, fn = function(mu) { -getNBGaussianLikelihood(x = as.numeric(curNormNBSummaries[sortBins][i, ]), mu = mu, k = binCounts, nullModel = curBinBounds, libFract = curNormNBSummaries$libFraction[i]) }, method = "L-BFGS-B", lower = min(uniformObjectiveEval$mu), upper = max(uniformObjectiveEval$mu))
      }
    )
    curNormNBSummaries$mean[i]=  curOptim$par
    curNormNBSummaries$ll[i]  = -curOptim$value
  }
  if (sum(curNormNBSummaries[[negativeControl]]) == 0){
    warning(sprintf("Cannot calculate logliklihood ratio or Z score without non-targeting guides and %s indicates there are no such guides",negativeControl))
  }
  #recalculate LL ratio and calculate a Z score for the mean WRT the observed mean expression of the non-targeting (NT) guides
  muNT = meanFunction(curNormNBSummaries$mean[curNormNBSummaries[[negativeControl]]]) # mean of the non-targeting guides mean expressions
  for (i in 1:nrow(curNormNBSummaries)){
    curNormNBSummaries$llRatio[i]=curNormNBSummaries$ll[i] -getNBGaussianLikelihood(x=as.numeric(curNormNBSummaries[sortBins][i,]), mu=muNT, k=binCounts, nullModel=curBinBounds, libFract = curNormNBSummaries$libFraction[i])
    if (!is.finite(curNormNBSummaries$llRatio[i])){
      warning(sprintf("Got non-finite llRatio: %f; %f - %f", curNormNBSummaries$llRatio[i],curNormNBSummaries$ll[i] , getNBGaussianLikelihood(x=as.numeric(curNormNBSummaries[sortBins][i,]), mu=muNT, k=binCounts, nullModel=curBinBounds, libFract = curNormNBSummaries$libFraction[i])))
    }
    curNormNBSummaries$Z[i]=curNormNBSummaries$mean[i]-muNT
  }
  if (any(curNormNBSummaries$llRatio[!is.na(curNormNBSummaries$llRatio)]<0)){
    warning("Some log-likelihood ratios are negative! (i.e. optimized mu is less likely to mu=0)")
  }
  return(curNormNBSummaries)
} 



#' Calculate Z-score scaling factors using non-targeting guides
#'
#' Calculates scaling factors to calibrate  element-wise Z-scores by repeatedly calculating a set of "null" Z-scores by repeatedly sampling the given numbers of non-targeting guides per element. This function is not normally used directly.
#' @param ntData data.frame containing the data for the non-targeting guides
#' @param uGuidesPerElement a unique vector of guide counts per element
#' @param mergeBy a character vector containing the header(s) that demarcate the screen/experiment/replicate ID(s)
#' @param ntSampleFold how many times to sample each non-targeting guide to make the Z score scale (defaults to 10)
#' @return a data.frame containing a Z-score scaling factor, one for every number of guides and unique entry in mergeBy
#' @export
#' @examples
#' fakeReadData = data.frame(id=rep(1:1000,2), expt=c(rep("e1",1000), rep("e2",1000)), 
#'                           A=rpois(2000, lambda = 100), B=rpois(2000, lambda = 100),
#'                           C=rpois(2000, lambda = 100), D=rpois(2000, lambda = 100),
#'                           E=rpois(2000, lambda = 100), F=rpois(2000, lambda = 100),
#'                           NotSorted=rpois(2000, lambda = 100), 
#'                           negControl = rep(rnorm(1000)>0,2), stringsAsFactors = FALSE)
#' expts = unique(fakeReadData["expt"]);
#' curSortBins = makeBinModel(data.frame(Bin = c("A","B","C","D","E","F"), fraction = rep(0.1,6)))
#' curSortBins = rbind(curSortBins, curSortBins) # duplicate and use same for both expts
#' curSortBins$expt = c(rep(expts$expt[1],6),rep(expts$expt[2],6))
#' guideHits = findGuideHitsAllScreens(experiments = expts, countDataFrame=fakeReadData,
#'                                     binStats = curSortBins, unsortedBin = "NotSorted",
#'                                     negativeControl="negControl")
#' guideZScales = getZScalesWithNTGuides(guideHits[guideHits$negControl,], uGuidesPerElement=1:10, 
#'   mergeBy=names(expts))
#' if(require("ggplot2")){
#'   p=ggplot(guideZScales, aes(x=numGuides, y=Zscale, colour=expt))+geom_point()+geom_line();
#'   print(p)
#' }
#' @importFrom reshape cast
getZScalesWithNTGuides = function(ntData, uGuidesPerElement, mergeBy, ntSampleFold=10){
  message(sprintf("Building background with %i non-targeting guides", nrow(ntData)))
  ntData = ntData[sample(1:nrow(ntData), nrow(ntData)*ntSampleFold, replace=TRUE),]
  zScales = data.frame();
  for(i in uGuidesPerElement){
    ntData = ntData[order(runif(nrow(ntData))),]
    for(sortBy in mergeBy){ ntData = ntData[order(ntData[[sortBy]]),]} #sort by screen, then by random
    ntData$groupID = floor((0:(nrow(ntData)-1))/i)
    #message(sprintf("Unique groups for %i guides per locus: %i", i, length(unique(ntData$groupID))))
    ntStats = as.data.frame(cast(ntData, as.formula(sprintf("%s + groupID ~ .", paste(mergeBy, collapse = " + "))), value="Z", fun.aggregate = function(x){return(list(numGuides = length(x), stoufferZ=combineZStouffer(x)))}))
    ntStats = ntStats[ntStats$numGuides==i,]
    ntStats = as.data.frame(cast(ntStats, as.formula(sprintf("%s ~ .", paste(mergeBy, collapse = " + "))), value="stoufferZ",fun.aggregate = function(x){sd(x,na.rm=TRUE)}))
    names(ntStats)[ncol(ntStats)]="Zscale"
    ntStats$numGuides=i;
    zScales = rbind(zScales, ntStats)
  }
  return(zScales)
}

#' Find active elements by sliding window
#'
#' Tests guides for activity by considering a sliding window across the tested region and including all guides within the window for the test.
#' @param experiments a data.frame containing the headers that demarcate the screen ID, which are all also present in normNBSummaries
#' @param normNBSummaries data.frame of guide-level statistics as generated by findGuideHits()
#' @param tails whether to test for increased expression ("upper"), decreased ("lower"), or both ("both"); (defaults to "both")
#' @param location the name of the column in normNBSummaries containing the chromosomal location (defaults to "pos")
#' @param chr the name of the column in normNBSummaries containing the chromosome name (defaults to "chr")
#' @param window the window width in base pairs (defaults to 500)
#' @param minGuides the minimum number of guides in a window required for a test (defaults to 5)
#' @param negativeControl the name in normNBSummaries containing a logical representing whether or not the guide is non-Targeting (i.e. a negative control guide).  Defaults to "NT"
#' @param ... other parameters for getZScalesWithNTGuides
#' @return a data.frame containing the statistics for all windows tested for activity
#' @export
#' @examples
#' fakeReadData = data.frame(id=rep(1:10000,2), expt=c(rep("e1",10000), rep("e2",10000)), 
#'                           A=rpois(20000, lambda = 100), B=rpois(20000, lambda = 100),
#'                           C=rpois(20000, lambda = 100), D=rpois(20000, lambda = 100),
#'                           E=rpois(20000, lambda = 100), F=rpois(20000, lambda = 100),
#'                           NotSorted=rpois(20000, lambda = 100), 
#'                           position = rep(c(rep(NA, 1000), (1:9000)*10 + 5E7),2), 
#'                           chr=rep(c(rep(NA, 1000), rep("chr1", 9000)),2), 
#'                           negControl = rep(c(rep(TRUE,1000),rep(FALSE,9000)),2), 
#'                           stringsAsFactors = FALSE)
#' #make one region an "enhancer" and "repressor" by skewing the reads 
#' enhancers = data.frame(name = c("enh","repr"), start = c(40000, 70000) + 5E7, 
#'   end = c(40500, 70500) + 5E7, chr="chr1")
#' enhancerData = findOverlappingElements(fakeReadData[!is.na(fakeReadData$position),], enhancers, 
#'   guides.pos = "position",elements.start = "start", elements.end = "end")
#' readSkew=1.2 # we will scale up/down the reads in ABC and DEF by this amount
#' enhancerData[enhancerData$name=="enh", c("D","E","F")] = 
#'   floor(readSkew* enhancerData[enhancerData$name=="enh", c("D","E","F")]);
#' enhancerData[enhancerData$name=="repr", c("D","E","F")] = 
#'   floor(enhancerData[enhancerData$name=="repr", c("D","E","F")]/readSkew);
#' enhancerData[enhancerData$name=="repr", c("A","B","C")] = 
#'   floor(readSkew* enhancerData[enhancerData$name=="repr", c("A","B","C")]);
#' enhancerData[enhancerData$name=="enh", c("A","B","C")] = 
#'   floor(enhancerData[enhancerData$name=="enh", c("A","B","C")]/readSkew);
#' #replace the original data for these elements
#' fakeReadData = rbind(fakeReadData[!(fakeReadData$id %in% enhancerData$id), ],
#'   enhancerData[names(fakeReadData)])
#' #make experiments and sorting strategy
#' expts = unique(fakeReadData["expt"]);
#' curSortBins = makeBinModel(data.frame(Bin = c("A","B","C","D","E","F"), fraction = rep(0.1,6)))
#' curSortBins = rbind(curSortBins, curSortBins) # duplicate and use same for both expts
#' curSortBins$expt = c(rep(expts$expt[1],6),rep(expts$expt[2],6))
#' guideHits = findGuideHitsAllScreens(experiments = expts, countDataFrame=fakeReadData,
#'                                     binStats = curSortBins, unsortedBin = "NotSorted",
#'                                     negativeControl="negControl")
#' if(require("ggplot2")){
#'   p=ggplot(guideHits, aes(x=position, y=Z))+geom_point() ; print(p)
#' }
#' tilingElementStats = getTilingElementwiseStats(experiments = expts, normNBSummaries = guideHits, 
#'   tails = "both", chr = "chr", location="position", negativeControl = "negControl")
#' if(require("ggplot2")){
#'   p=ggplot(tilingElementStats, aes(x=start, xend=end, y=significanceZ, yend=significanceZ, 
#'     colour=FDR<0.01))+geom_segment() ; print(p)
#' }
#' @import stats
getTilingElementwiseStats = function(experiments, normNBSummaries, tails="both", location="pos", chr="chr", window=500, minGuides=5, negativeControl="NT", ...){
  if(!location %in% names(normNBSummaries)){
    stop(sprintf("Column '%s' is missing from normNBSummaries; did you specify 'location='?", location))
  }
  if(!negativeControl %in% names(normNBSummaries)){
    stop(sprintf("Column '%s' is missing from normNBSummaries; did you specify 'negativeControl='?", negativeControl))
  }
  if(!"Z" %in% names(normNBSummaries)){
    stop(sprintf("Column 'Z' is missing from normNBSummaries"))
  }
  if(!chr %in% names(normNBSummaries)){
    stop(sprintf("Column '%s' is missing from normNBSummaries", chr))
  }
  experiments = unique(experiments);
  mergeBy = names(experiments);
  ntData = normNBSummaries[normNBSummaries[[negativeControl]],]
  if(nrow(ntData)<10){
    stop(sprintf("You have too few negative controls (only %i). This can indicate a mistake in the input labels for %s", nrow(ntData), negativeControl))
  }
  normNBSummaries = normNBSummaries[!normNBSummaries[[negativeControl]],]
  if(nrow(normNBSummaries)<1){
    stop(sprintf("You have no non-negative controls. Usually this is a mistake in the input labels for %s", negativeControl))
  }
  elementwiseStats = data.frame();
  for (i in 1:nrow(experiments)){
    for (curChr in unique(normNBSummaries[[chr]])){
      curData = merge(experiments[i,, drop=FALSE],normNBSummaries[normNBSummaries[[chr]]==curChr,], by=mergeBy)
      curData = curData[order(curData[[location]]),]
      lagging=1
      leading=1
      while (leading <=nrow(curData)){
        while (curData[[location]][lagging] + window < curData[[location]][leading]){
          lagging = lagging +1;
        }
        while (leading+1 <=nrow(curData) && curData[[location]][lagging] + window >= curData[[location]][leading+1]){
          leading = leading +1;
        }
        #message(sprintf("%i:%i %s:%i-%i",lagging,leading, curChr, curData[[location]][lagging], curData[[location]][leading]))
        if (leading-lagging +1 >= minGuides){
          curRes = data.frame(chr = curChr, start = curData[[location]][lagging], end = curData[[location]][leading], numGuides = leading-lagging+1, stoufferZ=combineZStouffer(curData$Z[lagging:leading]), meanZ=mean(curData$Z[lagging:leading]))
          for (k in names(experiments)){
            curRes[k]=experiments[i,k];
          }
          elementwiseStats = rbind(elementwiseStats, curRes)
        }
        leading = leading + 1;
      }
    }
  }
  elementwiseStats = elementwiseStats[order(elementwiseStats$stoufferZ),]
  #head(elementwiseStats)
  #calibrate Z scores
  uGuidesPerElement = sort(unique(elementwiseStats$numGuides))
  zScales = data.frame();
  #build background
  message("Building null model")
  zScales = getZScalesWithNTGuides(ntData,uGuidesPerElement, mergeBy, ...)
  if (!all(c(mergeBy,"numGuides") %in% names(elementwiseStats))){
    stop(sprintf("{%s} not all present in names of elementwiseStats: {%s}", paste(c(mergeBy,"numGuides"), collapse=", "), paste(names(elementwiseStats), collapse=", ")));
  }
  if (!all(c(mergeBy,"numGuides") %in% names(zScales))){
    stop(sprintf("{%s} not all present in names of zScales: {%s}", paste(c(mergeBy,"numGuides"), collapse=", "), paste(names(zScales), collapse=", ")));
  }
  elementwiseStats = merge(elementwiseStats, zScales, by=c(mergeBy,"numGuides"))
  message("Calculating P-values")
  elementwiseStats$significanceZ = elementwiseStats$stoufferZ/elementwiseStats$Zscale;
  if (tails=="both" || tails =="lower"){
    elementwiseStats$p.value = pnorm(elementwiseStats$significanceZ)
  }
  if (tails=="both"){
    elementwiseStats$p.value[elementwiseStats$significanceZ>0] = pnorm(-elementwiseStats$significanceZ[elementwiseStats$significanceZ>0])
  }
  if (tails=="upper"){
    elementwiseStats$p.value = pnorm(-elementwiseStats$significanceZ)
  }
  
  elementwiseStats$FDR = calcFDRByExperiment(experiments, elementwiseStats, tails)
  elementwiseStats$Zscale=NULL
  return(elementwiseStats);
}

#' FDR correction per experiment
#'
#' Perform Benjamini-Hochberg FDR correction of p-values within each experiment. The returned values correspond to Q values. This is run automatically within getTilingElementwiseStats and getElementwiseStats, and doesn't generally need to be used directly.
#' @param experiments a data.frame containing the headers that demarcate the screen ID, which are all also present in normNBSummaries
#' @param x data.frame of element-level statistics, including columns for every column in 'experiments' and a column named 'p.value'
#' @param tails whether to test for increased expression ("upper"), decreased ("lower"), or both ("both"); (defaults to "both")
#' @return a numerical vector containing B-H Q values corrected separately for every experiment.
#' @export
#' @examples
#' fakeReadData = data.frame(id=rep(1:10000,2), expt=c(rep("e1",10000), rep("e2",10000)), 
#'                           A=rpois(20000, lambda = 100), B=rpois(20000, lambda = 100),
#'                           C=rpois(20000, lambda = 100), D=rpois(20000, lambda = 100),
#'                           E=rpois(20000, lambda = 100), F=rpois(20000, lambda = 100),
#'                           NotSorted=rpois(20000, lambda = 100), 
#'                           position = rep(c(rep(NA, 1000), (1:9000)*10 + 5E7),2), 
#'                           chr=rep(c(rep(NA, 1000), rep("chr1", 9000)),2), 
#'                           negControl = rep(c(rep(TRUE,1000),rep(FALSE,9000)),2), 
#'                           stringsAsFactors = FALSE)
#' #make one region an "enhancer" and "repressor" by skewing the reads 
#' enhancers = data.frame(name = c("enh","repr"), start = c(40000, 70000) + 5E7, 
#'   end = c(40500, 70500) + 5E7, chr="chr1")
#' enhancerData = findOverlappingElements(fakeReadData[!is.na(fakeReadData$position),], enhancers, 
#'   guides.pos = "position",elements.start = "start", elements.end = "end")
#' readSkew=1.2 # we will scale up/down the reads in ABC and DEF by this amount
#' enhancerData[enhancerData$name=="enh", c("D","E","F")] = 
#'   floor(readSkew* enhancerData[enhancerData$name=="enh", c("D","E","F")]);
#' enhancerData[enhancerData$name=="repr", c("D","E","F")] = 
#'   floor(enhancerData[enhancerData$name=="repr", c("D","E","F")]/readSkew);
#' enhancerData[enhancerData$name=="repr", c("A","B","C")] = 
#'   floor(readSkew* enhancerData[enhancerData$name=="repr", c("A","B","C")]);
#' enhancerData[enhancerData$name=="enh", c("A","B","C")] = 
#'   floor(enhancerData[enhancerData$name=="enh", c("A","B","C")]/readSkew);
#' #replace the original data for these elements
#' fakeReadData = rbind(fakeReadData[!(fakeReadData$id %in% enhancerData$id), ],
#'   enhancerData[names(fakeReadData)])
#' #make experiments and sorting strategy
#' expts = unique(fakeReadData["expt"]);
#' curSortBins = makeBinModel(data.frame(Bin = c("A","B","C","D","E","F"), fraction = rep(0.1,6)))
#' curSortBins = rbind(curSortBins, curSortBins) # duplicate and use same for both expts
#' curSortBins$expt = c(rep(expts$expt[1],6),rep(expts$expt[2],6))
#' guideHits = findGuideHitsAllScreens(experiments = expts, countDataFrame=fakeReadData,
#'                                     binStats = curSortBins, unsortedBin = "NotSorted",
#'                                     negativeControl="negControl")
#' tilingElementStats = getTilingElementwiseStats(experiments = expts, normNBSummaries = guideHits, 
#'   tails = "both", chr = "chr", location="position", negativeControl = "negControl")
#' tilingElementStats$Q = calcFDRByExperiment(expts, tilingElementStats,"both") 
#' if(require("ggplot2")){
#'   p=ggplot(tilingElementStats, aes(x=FDR, y=Q)) +geom_point()+geom_abline(intercept=0, slope=1)+
#'     scale_y_log10()+scale_x_log10(); print(p)
#' }
#' @import stats
calcFDRByExperiment = function(experiments, x, tails){
  allFDRs = rep(NA, nrow(x));
  for (i in 1:nrow(experiments)){
    curRows = rep(TRUE, nrow(x));
    for(n in names(experiments)){
      curRows = curRows & x[n]==experiments[i,n]
    }
    allFDRs[curRows] = p.adjust(x$p.value[curRows], method = "BH", n = sum(curRows) * (1 + ((tails=="both")) ))
  }
  return(allFDRs);
}


#' Find active elements by annotation
#'
#' Tests guides for activity by considering a set of provided regulatory elements within the region and considering all guides within each region for the test.
#' @param experiments a data.frame containing the headers that demarcate the screen ID, which are all also present in normNBSummaries
#' @param normNBSummaries data.frame of guide-level statistics as generated by findGuideHits()
#' @param elementIDs the names of one or more columns within guideLevelStats that contain the element annotations.
#' @param tails whether to test for increased expression ("upper"), decreased ("lower"), or both ("both"); (defaults to "both")
#' @param negativeControl the name in normNBSummaries containing a logical representing whether or not the guide is non-Targeting (i.e. a negative control guide).  Defaults to "NT"
#' @param ... other parameters for getZScalesWithNTGuides
#' @return a data.frame containing the statistics for all elements
#' @export
#' @examples
#' fakeReadData = data.frame(id=rep(1:10000,2), expt=c(rep("e1",10000), rep("e2",10000)), 
#'                           A=rpois(20000, lambda = 100), B=rpois(20000, lambda = 100),
#'                           C=rpois(20000, lambda = 100), D=rpois(20000, lambda = 100),
#'                           E=rpois(20000, lambda = 100), F=rpois(20000, lambda = 100),
#'                           NotSorted=rpois(20000, lambda = 100), 
#'                           position = rep(c(rep(NA, 1000), (1:9000)*10 + 5E7),2), 
#'                           chr=rep(c(rep(NA, 1000), rep("chr1", 9000)),2), 
#'                           negControl = rep(c(rep(TRUE,1000),rep(FALSE,9000)),2), 
#'                           stringsAsFactors = FALSE)
#' #make one region an "enhancer" and "repressor" by skewing the reads 
#' enhancers = data.frame(name = c("enh","repr"), start = c(40000, 70000) + 5E7, 
#'   end = c(40500, 70500) + 5E7, chr="chr1")
#' enhancerData = findOverlappingElements(fakeReadData[!is.na(fakeReadData$position),], 
#'   enhancers, guides.pos = "position",elements.start = "start", elements.end = "end")
#' readSkew=1.2 # we will scale up/down the reads in ABC and DEF by this amount
#' enhancerData[enhancerData$name=="enh", c("D","E","F")] = 
#'   floor(readSkew* enhancerData[enhancerData$name=="enh", c("D","E","F")]);
#' enhancerData[enhancerData$name=="repr", c("D","E","F")] = 
#'   floor(enhancerData[enhancerData$name=="repr", c("D","E","F")]/readSkew);
#' enhancerData[enhancerData$name=="repr", c("A","B","C")] = 
#'   floor(readSkew* enhancerData[enhancerData$name=="repr", c("A","B","C")]);
#' enhancerData[enhancerData$name=="enh", c("A","B","C")] = 
#'   floor(enhancerData[enhancerData$name=="enh", c("A","B","C")]/readSkew);
#' #replace the original data for these elements
#' fakeReadData = rbind(fakeReadData[!(fakeReadData$id %in% enhancerData$id), ],
#'   enhancerData[names(fakeReadData)])
#' #make experiments and sorting strategy
#' expts = unique(fakeReadData["expt"]);
#' curSortBins = makeBinModel(data.frame(Bin = c("A","B","C","D","E","F"), fraction = rep(0.1,6)))
#' curSortBins = rbind(curSortBins, curSortBins) # duplicate and use same for both expts
#' curSortBins$expt = c(rep(expts$expt[1],6),rep(expts$expt[2],6))
#' guideHits = findGuideHitsAllScreens(experiments = expts, countDataFrame=fakeReadData,
#'                                     binStats = curSortBins, unsortedBin = "NotSorted",
#'                                     negativeControl="negControl")
#' 
#' potentialEnhancers = rbind(enhancers, data.frame(name = sprintf("EC%i", 1:16), 
#'   start = (5:20)*6000 + 5E7, end = (5:20)*6000 + 5E7+500, chr="chr1"))
#' #annotate guides with elements, but first get all non-targeting guides
#' guideHitsAnnotated = guideHits[is.na(guideHits$position),];
#' guideHitsAnnotated$name=NA; guideHitsAnnotated$start=NA; 
#' guideHitsAnnotated$end=NA; guideHitsAnnotated$chr=NA;
#' guideHitsAnnotated = rbind(guideHitsAnnotated, 
#'   findOverlappingElements(guideHits[!is.na(guideHits$position),], potentialEnhancers, 
#'     guides.pos = "position",elements.start = "start", elements.end = "end") )
#' allElementHits = getElementwiseStats(experiments = expts, normNBSummaries = guideHitsAnnotated, 
#'   elementIDs = "name", tails = "both", negativeControl = "negControl")
#' allElementHits = merge(allElementHits, potentialEnhancers, by="name")
#' if(require("ggplot2")){
#'   p=ggplot(allElementHits, aes(x=start, xend=end, y=significanceZ, yend=significanceZ, 
#'     colour=FDR<0.01, label=name)) + geom_segment()+
#'     geom_text(data= allElementHits[allElementHits$FDR<0.01,], colour="black") + 
#'     facet_grid(expt ~ .); print(p)
#' }
#' @importFrom reshape cast
#' @import stats
getElementwiseStats = function(experiments, normNBSummaries, elementIDs, tails="both", negativeControl="NT",...){
  if(!negativeControl %in% names(normNBSummaries)){
    stop(sprintf("Column '%s' is missing from normNBSummaries; did you specify 'negativeControl='?", negativeControl))
  }
  experiments = unique(experiments);
  mergeBy = names(experiments);
  ntData = normNBSummaries[normNBSummaries[[negativeControl]],]
  normNBSummaries = normNBSummaries[!normNBSummaries[[negativeControl]],]
  elementwiseStats = cast(normNBSummaries[!apply(is.na(normNBSummaries[elementIDs]), 1, any),], as.formula(sprintf("%s ~ .", paste(c(elementIDs, mergeBy), collapse = " + "))), value="Z", fun.aggregate = function(x){return(list(numGuides = length(x), stoufferZ=combineZStouffer(x), meanZ=mean(x)))})
  elementwiseStats = elementwiseStats[order(elementwiseStats$stoufferZ),]
  #head(elementwiseStats)
  #calibrate Z scores
  uGuidesPerElement = sort(unique(elementwiseStats$numGuides))
  #build background
  zScales = getZScalesWithNTGuides(ntData,uGuidesPerElement, mergeBy, ...)
  elementwiseStats = merge(elementwiseStats, zScales, by=c(mergeBy,"numGuides"))
  elementwiseStats$significanceZ = elementwiseStats$stoufferZ/elementwiseStats$Zscale;
  if (tails=="both" || tails =="lower"){
    elementwiseStats$p.value = pnorm(elementwiseStats$significanceZ)
  }
  if (tails=="both"){
    elementwiseStats$p.value[elementwiseStats$significanceZ>0] = pnorm(-elementwiseStats$significanceZ[elementwiseStats$significanceZ>0])
  }
  if (tails=="upper"){
    elementwiseStats$p.value = pnorm(-elementwiseStats$significanceZ)
  }
  elementwiseStats$FDR = calcFDRByExperiment(experiments, elementwiseStats, tails)
  elementwiseStats$Zscale=NULL
  return(elementwiseStats);
}

#' Calculate guide-level stats for multiple experiments
#'
#' Uses findGuideHits to find guide-level stats for each unique entry in 'experiments'.
#' @param experiments a data.frame containing the headers that demarcate the screen ID, which are all also present in countDataFrame and binStats
#' @param countDataFrame a table containing one column for each bin (A-F) and another column for non-targeting guide (logical-"NT"), and unsorted abundance (NS), as well as columns corresponding to those in  'experiments' 
#' @param binStats a bin model as created by makeBinModel, as well as columns corresponding to those in  'experiments'
#' @param ... other parameters for findGuideHits
#' @return guide-level stats for all experiments
#' @export
#' @examples
#'  fakeReadData = data.frame(id=rep(1:1000,2), expt=c(rep("e1",1000), rep("e2",1000)), 
#'                           A=rpois(2000, lambda = 100), B=rpois(2000, lambda = 100),
#'                           C=rpois(2000, lambda = 100), D=rpois(2000, lambda = 100),
#'                           E=rpois(2000, lambda = 100), F=rpois(2000, lambda = 100),
#'                           NotSorted=rpois(2000, lambda = 100), 
#'                           negControl = rep(rnorm(1000)>0,2), stringsAsFactors = FALSE)
#' expts = unique(fakeReadData["expt"]);
#' curSortBins = makeBinModel(data.frame(Bin = c("A","B","C","D","E","F"), fraction = rep(0.1,6)))
#' curSortBins = rbind(curSortBins, curSortBins) # duplicate and use same for both expts
#' curSortBins$expt = c(rep(expts$expt[1],6),rep(expts$expt[2],6))
#' guideHits = findGuideHitsAllScreens(experiments = expts, countDataFrame=fakeReadData,
#'                                     binStats = curSortBins, unsortedBin = "NotSorted",
#'                                     negativeControl="negControl")
findGuideHitsAllScreens = function(experiments, countDataFrame, binStats, ...){
  if (!"Bin" %in% names(binStats)){
    if (!"bin" %in% names(binStats)){
      warning("'Bin' column not found in binStats; using 'bin' instead")
      binStats$Bin = binStats$bin;
    }else{
      stop("No 'Bin' column in binStats!")
    }
  }
  if (!"fraction" %in% names(binStats)){
    stop("No 'fraction' column in binStats!")
  }
  experiments = unique(experiments);
  mergeBy = names(experiments);
  if (!all(mergeBy %in% names(countDataFrame))){
    message(names(experiments))
    message(names(countDataFrame))
    stop("Columns from 'experiments' missing from 'countDataFrame'");
  }
  if (!all(mergeBy %in% names(binStats))){
    stop("Columns from 'experiments' missing from 'binStats'");
  }
  experiments[ncol(experiments)+1]=1:nrow(experiments);
  idCol = names(experiments)[ncol(experiments)];
  countDataFrame = merge(countDataFrame, experiments, by=mergeBy);
  binStats = merge(binStats, experiments, by=mergeBy);
  
  allSummaries = data.frame()
  for(j in 1:nrow(experiments)){
    normNBSummaries = findGuideHits(
      countDataFrame[countDataFrame[[idCol]]==experiments[[idCol]][j],], 
      #makeBinModel(binStats[binStats[[idCol]]==experiments[[idCol]][j] , c("Bin","fraction")]), ...)
      binStats[binStats[[idCol]]==experiments[[idCol]][j] , ], ...)
    allSummaries = rbind(allSummaries,normNBSummaries)
  }
  allSummaries[[idCol]]=NULL;
  return(allSummaries);
}


