
library(codetools)

getNBGaussianLikelihood = function(x, mu, k, sigma=1, nullModel, libFract){
  #x= vector of guide counts per bin
  #mu = mean of distribution
  #k = vector of total counts per bin
  #sigma = sigma for normal distribution
  #nullModel = bin bounds for null model
  #calculate the probabilities of being within each bin
  binFractions = pnorm(nullModel$binEndZ-mu) - pnorm(nullModel$binStartZ-mu)
  #message(sprintf("binFractions = %s",paste(as.character((binFractions)),collapse=", ")))
  binFractions = binFractions*libFract/(nullModel$binEndQ - nullModel$binStartQ)
  #message(sprintf("scaled binFractions = %s",paste(as.character((binFractions)),collapse=", ")))
  likelihood =0;
  #message(sprintf("observed fractions = %s",paste(as.character(((x/k))),collapse=", ")))
  for (i in 1:6){
    #dnbinom(x = number of reads for this guide, size = number of reads total, prob= probability of getting a read at each drawing)
    likelihood = likelihood  + dnbinom(x=x[i], size=k[i], prob =1- binFractions[i], log=T)
  }
  return(likelihood)
}
checkUsage(getNBGaussianLikelihood)

makeBinModel = function(curBinBounds,tailP=0.001){
  curBinBounds$binStartQ[curBinBounds$Bin=="A"]=tailP #A
  curBinBounds$binEndQ[curBinBounds$Bin=="A"] = curBinBounds$binStartQ[curBinBounds$Bin=="A"] + curBinBounds$fraction[curBinBounds$Bin=="A"];
  curBinBounds$binStartQ[curBinBounds$Bin=="B"] = curBinBounds$binEndQ[curBinBounds$Bin=="A"]; #B
  curBinBounds$binEndQ[curBinBounds$Bin=="B"] = curBinBounds$binStartQ[curBinBounds$Bin=="B"] + curBinBounds$fraction[curBinBounds$Bin=="B"];
  curBinBounds$binStartQ[curBinBounds$Bin=="C"] = curBinBounds$binEndQ[curBinBounds$Bin=="B"]; #C
  curBinBounds$binEndQ[curBinBounds$Bin=="C"] = curBinBounds$binStartQ[curBinBounds$Bin=="C"] + curBinBounds$fraction[curBinBounds$Bin=="C"];
  curBinBounds$binEndQ[curBinBounds$Bin=="F"]=1-tailP # F
  curBinBounds$binStartQ[curBinBounds$Bin=="F"] = curBinBounds$binEndQ[curBinBounds$Bin=="F"] - curBinBounds$fraction[curBinBounds$Bin=="F"];
  curBinBounds$binEndQ[curBinBounds$Bin=="E"] = curBinBounds$binStartQ[curBinBounds$Bin=="F"]; #E
  curBinBounds$binStartQ[curBinBounds$Bin=="E"] = curBinBounds$binEndQ[curBinBounds$Bin=="E"] - curBinBounds$fraction[curBinBounds$Bin=="E"];
  curBinBounds$binEndQ[curBinBounds$Bin=="D"] = curBinBounds$binStartQ[curBinBounds$Bin=="E"]; #D
  curBinBounds$binStartQ[curBinBounds$Bin=="D"] = curBinBounds$binEndQ[curBinBounds$Bin=="D"] - curBinBounds$fraction[curBinBounds$Bin=="D"];
  curBinBounds$binStartZ = qnorm(curBinBounds$binStartQ)
  curBinBounds$binEndZ = qnorm(curBinBounds$binEndQ)
  return(curBinBounds)
}
checkUsage(makeBinModel)

findGuideHits = function(countTable, curBinBounds, pseudocount=10, meanFunction = mean){
  allBins = c("A","B","C","D","E","F")
  #for(b in allBins){
  #  countTable[b]=countTable[b]+round(pseudocount*sum(countTable);#add a pseudocount
  #}
  countTable[allBins]=countTable[allBins]+pseudocount;#add a pseudocount
  countTable["NS"]=countTable["NS"]+pseudocount;
  curNormNBSummaries = countTable
  allDataCounts$libFraction = allDataCounts$NS/sum(allDataCounts$NS,na.rm=T)
  
  curNormNBSummaries$libFraction = curNormNBSummaries$NS/sum(curNormNBSummaries$NS,na.rm=T)
  binCounts = apply(curNormNBSummaries[allBins],2,function(x){sum(x, na.rm = T)})
  
  #for each guide, find the optimal mu given the count data and bin percentages
  for (i in 1:nrow(curNormNBSummaries)){
    #interval: The probability of observing a guide outside of this interval in one of the non-terminal bins is very unlikely, and so estimating a true mean for these is too difficult anyway. Besides, we get some local optima at the extremes for sparsely sampled data.
    temp = optimize(f=function(mu){getNBGaussianLikelihood(x=as.numeric(curNormNBSummaries[allBins][i,]), mu=mu, k=binCounts, nullModel=curBinBounds, libFract = curNormNBSummaries$libFraction[i])}, interval=c(-4,4), maximum = T)
    #TODO: This function can get stuck in local minima. Pseudocounts help prevent this, but it can still happen, resulting in a negative logliklihood ratio (i.e. Null is more likely than optimized alternate).  Usually this happens close to an effect size of 0.  I should still explore other optimization functions (e.g. optim)
    curNormNBSummaries$mean[i]=temp$maximum
    curNormNBSummaries$ll[i]=temp$objective
  }
  #recalculate LL ratio and calculate a Z score for the mean WRT the observed mean expression of the non-targeting (NT) guides
  muNT = meanFunction(curNormNBSummaries$mean[curNormNBSummaries$NT]) # mean of the non-targeting guides mean expressions
  for (i in 1:nrow(curNormNBSummaries)){
    curNormNBSummaries$llRatio[i]=curNormNBSummaries$ll[i] -getNBGaussianLikelihood(x=as.numeric(curNormNBSummaries[allBins][i,]), mu=muNT, k=binCounts, nullModel=curBinBounds, libFract = curNormNBSummaries$libFraction[i])
    curNormNBSummaries$Z[i]=curNormNBSummaries$mean[i]-muNT
  }
  return(curNormNBSummaries)
}
checkUsage(findGuideHits)

getElementwiseStats = function(screens, normNBSummaries, elementIDs, ntSampleFold = 10, tails="both"){
  screens = unique(screens);
  mergeBy = names(screens);
  elementwiseStats = cast(normNBSummaries[!apply(is.na(normNBSummaries[elementIDs]), 1, any),], as.formula(sprintf("%s ~ .", paste(c(elementIDs, mergeBy), collapse = " + "))), value="Z", fun.aggregate = function(x){return(list(numGuides = length(x), stoufferZ=combineZStouffer(x), meanZ=mean(x)))})
  elementwiseStats = elementwiseStats[order(elementwiseStats$stoufferZ),]
  #head(elementwiseStats)
  #calibrate Z scores
  uGuideLens = sort(unique(elementwiseStats$numGuides))
  zScales = data.frame();
  #build background
  message(sprintf("Building background with %i non-targeting guides", sum(normNBSummaries$NT)))
  ntData = normNBSummaries[normNBSummaries$NT,]
  ntData = ntData[sample(1:nrow(ntData), nrow(ntData)*ntSampleFold, replace=T),]
  for(i in uGuideLens){
    ntData = ntData[order(runif(nrow(ntData))),]
    for(sortBy in mergeBy){ ntData = ntData[order(ntData[sortBy]),]} #sort by screen, then by random
    ntData$groupID = floor((0:(nrow(ntData)-1))/i)
    message(sprintf("Unique groups for %i guides per locus: %i", i, length(unique(ntData$groupID))))
    #message(str(ntData))
    ntStats = as.data.frame(cast(ntData, as.formula(sprintf("%s + groupID ~ .", paste(mergeBy, collapse = " + "))), value="Z", fun.aggregate = function(x){return(list(numGuides = length(x), stoufferZ=combineZStouffer(x)))}))
    #message(str(ntStats))
    ntStats = ntStats[ntStats$numGuides==i,]
    #message(str(ntStats))
    ntStats = as.data.frame(cast(ntStats, as.formula(sprintf("%s ~ .", paste(mergeBy, collapse = " + "))), value="stoufferZ",fun.aggregate = function(x){sd(x,na.rm=T)}))
    names(ntStats)[ncol(ntStats)]="Zscale"
    ntStats$numGuides=i;
    #message(str(ntStats))
    zScales = rbind(zScales, ntStats)
  }
  elementwiseStats = merge(elementwiseStats, zScales, by=c(mergeBy,"numGuides"))
  elementwiseStats$rescaledZ = elementwiseStats$stoufferZ/elementwiseStats$Zscale;
  if (tails=="both" || tails =="lower"){
    elementwiseStats$p.value = pnorm(elementwiseStats$rescaledZ)
  }
  if (tails=="both"){
    elementwiseStats$p.value[elementwiseStats$rescaledZ>0] = pnorm(-elementwiseStats$rescaledZ[elementwiseStats$rescaledZ>0])
  }
  if (tails=="upper"){
    elementwiseStats$p.value = pnorm(-elementwiseStats$rescaledZ)
  }
  elementwiseStats$FDR = p.adjust(elementwiseStats$p.value, method = "BH", n = nrow(elementwiseStats) + ((tails=="both") * nrow(elementwiseStats)))
  elementwiseStats$Zscale=NULL
  return(elementwiseStats);
}
checkUsage(getElementwiseStats)

findGuideHitsAllScreens = function(screens, countDataFrame, binStats, ...){
  if (!"Bin" %in% names(binStats)){
    if (!"bin" %in% names(binStats)){
      warning("'Bin' column not found in binStats; using 'bin' instead")
      binStats$Bin = binStats$bin;
    }else{
      stop("No 'Bin' column in binStats!")
    }
  }
  if (!"frequency" %in% names(binStats)){
    stop("No 'frequency' column in binStats!")
  }
  if (!all(c("A","B","C","D","E","F","NS") %in% names(countDataFrame))){
    stop("Not all bins (A-F, NS) present in countDataFrame!")
  }
  screens = unique(screens);
  mergeBy = names(screens);
  if (!all(mergeBy %in% names(countDataFrame))){
    stop("Columns from 'screens' missing from 'countDataFrame'");
  }
  if (!all(mergeBy %in% names(binStats))){
    stop("Columns from 'screens' missing from 'binStats'");
  }
  screens[ncol(screens)+1]=1:nrow(screens);
  idCol = names(screens)[ncol(screens)];
  countDataFrame = merge(countDataFrame, screens, by=mergeBy);
  binStats = merge(binStats, screens, by=mergeBy);
  
  allSummaries = data.frame()
  for(j in 1:nrow(screens)){
    normNBSummaries = findGuideHits(
      countDataFrame[countDataFrame[[idCol]]==screens[[idCol]][j],], 
      makeBinModel(binStats[binStats[[idCol]]==screens[[idCol]][j] , c("Bin","fraction")]), ...)
    allSummaries = rbind(allSummaries,normNBSummaries)
  }
  allSummaries[[idCol]]=NULL;
  return(allSummaries);
}
checkUsage(findGuideHitsAllScreens)

