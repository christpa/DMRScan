% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmrscan_foo.R
\name{dmrscan}
\alias{dmrscan}
\title{DMR Scan function}
\usage{
dmrscan(observations, windowSize, windowThreshold = NULL, chr = NULL,
  pos = NULL, maxGap = 500, ...)
}
\arguments{
\item{observations}{An object of either;  \code{\link{GRangesList}} made by 
\code{\link{makeCpGregions}}, a vector of the test statistic, a \code{\link{GRanges}} object, 
or a "minfi" object (soon to be supported).}

\item{windowSize}{A sequence of windowSizes for the slidingWindow. Must be an 
integer vector, with equal length as the number of windows.}

\item{windowThreshold}{Optional argument with corresponding cut-off for 
each window. Will be estimated if not supplied.}

\item{chr}{A vector of chromosomal position. Only used when the observations
vector is a matrix of test statistic.}

\item{pos}{A vector of genomic coordinates for the CpGs to match the chr argument}

\item{maxGap}{The maximum allowed gap between two CpGs within the same region.}

\item{...}{Optional arguments to be passed to \code{\link{estimateThreshold}}, 
if no grid is specified.}
}
\value{
An object of type \code{\link{GRanges}} with significantly differentially
}
\description{
DMR Scan function
}
\examples{
## methylation data from chromosome 22 
data(DMRScan.methylationData)
## phenotype (end-point for methylation data)

data(DMRScan.phenotypes) 

## Test for an association between phenotype and methylation
test.statistics <- apply(DMRScan.methylationData,1,function(x,y)
  summary(glm(y ~ x, family = binomial(link = "logit")))$coefficients[2,3],
                                                   y = DMRScan.phenotypes)
## Set chromosomal position to each test-statistic
positions <- data.frame(matrix(as.integer(unlist(strsplit(names(test.statistics), 
                               split="chr|[.]"))), ncol = 3, byrow = TRUE))[,-1]
## Set clustering features 
min.cpg <- 4  ## Minimum number of CpGs in a tested cluster
## Maximum distance (in base-pairs) within a cluster 
## before it is broken up into two separate cluster 
max.gap <- 750  

## Identify all clusters, and generate a list for each cluster
regions <- makeCpGregions(observations = test.statistics, 
                          chr = positions[,1], pos = positions[,2], 
                          maxGap = max.gap, minCpG = min.cpg)
## Number of CpGs in the slidingWindows, can be either a single number 
## or a sequence of windowSizes
windowSizes <- 3:7 
nCpG        <- sum(sapply(regions,length)) ## Number of CpGs to be tested

# Estimate the windowThreshold, based on the number of CpGs and windowSizes
windowThresholds <- estimateWindowThreshold(nProbe = nCpG, 
               windowSize = windowSizes, method = "sampling", mcmc = 10000)
## Run the slidingWindow
DMRScanResults   <- dmrscan(observations = regions, 
                            windowSize = windowSizes, 
                            windowThreshold = windowThresholds)
## Print the result
print(DMRScanResults)

}
\keyword{DMRScan}
