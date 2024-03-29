% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MAUDE.R
\name{makeBinModel}
\alias{makeBinModel}
\title{Create a bin model for a single experiment}
\usage{
makeBinModel(curBinBounds, tailP = 0.001)
}
\arguments{
\item{curBinBounds}{a data.frame containing two columns: Bin (must be {A,B,C,D,E,F}), and fraction (the fractions of the total captured by each bin)}

\item{tailP}{the fraction of the tails of the distribution not captured in any bin (defaults to 0.001)}
}
\value{
returns a data.frame with additional columns including the bin starts and ends in Z-score space, and in quantile space.
}
\description{
Provided with the fractions captured by each bin, creates a bin model for use with MAUDE analysis, assuming 3 contiguous bins on the tails of the distribution. You can easily remove or rename bins after they have been created with this function. An example is provided in the BACH2 Vignette.
}
\examples{
#generally, the bin bounds are retrieved from the FACS data
binBounds = makeBinModel(data.frame(Bin=c("A","B","C","D","E","F"), fraction=rep(0.1,6))) 
if(require("ggplot2")){
  p = ggplot() + 
    geom_vline(xintercept= sort(unique(c(binBounds$binStartZ, binBounds$binEndZ))),colour="gray")+
    theme_classic() + xlab("Target expression") + 
    geom_segment(data=binBounds, aes(x=binStartZ, xend=binEndZ, colour=Bin, y=0, yend=0), 
      size=5, inherit.aes = FALSE); 
  print(p)
}
}
