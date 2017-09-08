#' DMR Scan function
#' 
#' @name dmrscan
#' @rdname DMRScan_slidingWindow
#' @param observations An object of type RegionList 
#' @param windowSize A sequence of windowSizes for the slidingWindow, 
#' must be an integer 
#' @param windowThreshold Optional argument with corresponding cut-off for 
#' each window. Will be estimated if not supplied.
#' @param ... Optional arguments to be pased to estimate_windowThreshold(), 
#' if no grid is specified.
#' @return An object of type RegionList with signficantly differentially 
# methylated regions
#' @keywords DMRScan
#' @importFrom GenomicRanges GRanges
#' @importFrom stats acf median pnorm na.omit
#' @export
#' @examples
#' ## nProbeoad methylation data from chromosome 22
#' data(DMRScan.methylationData) 
#' ## nProbeoad phenotype (end-point for methylation data)
#' data(DMRScan.phenotypes) 
#' 
#' ## Test for an association between phenotype and Methylation
#' test.statistics <- apply(DMRScan.methylationData,1,function(x,y)
#'   summary(glm(y ~ x, family = binomial(link = "logit")))$coefficients[2,3],
#'                                                    y = DMRScan.phenotypes)
#' ## Set chromosomal position to each test-statistic
#' positions <- data.frame(matrix(as.integer(unlist(strsplit(names(test.statistics), split="chr|[.]"))), ncol = 3, byrow = TRUE))[,-1]
#' ## Set clustering features 
#' min.cpg <- 4  ## Minimum number of CpGs in a tested cluster
#' ## Maxium distance (in base-pairs) within a cluster 
#' ## before it is broken up into two seperate cluster 
#' max.gap <- 750  
#' 
# 
#' ## Identify all clusters, and generate a list for each cluster
#' regions <- makeCpGregions(observations = test.statistics, 
#'                           chr = positions[,1], pos = positions[,2], 
#'                           maxGap = max.gap, minCpG = min.cpg)
#' ## Number of CpGs in the slidingWindows, can be either a single number 
#' ## or a sequence of windowSizes
#' windowSizes <- 3:7 
#' nCpG        <- nCpG(regions) ## Number of CpGs to be tested
#' 
#' # Estimate the windowThreshold, based on the number of CpGs and windowSizes
#' windowThresholds <- estimateWindowThreshold(nProbe = nCpG, 
#'                windowSize = windowSizes, method = "sampling", mcmc = 10000)
#' ## Run the slidingWindow
#' DMRScanResults   <- dmrscan(observations = regions, 
#'                             windowSize = windowSizes, 
#'                             windowThreshold = windowThresholds)
#' ## Print the result
#' print(DMRScanResults)
#' 
dmrscan <- function(observations,windowSize,windowThreshold=NULL,...){

    nProbe      <- nCpG(observations)
    alpha       <- 0.05

    windowSize  <- sort(windowSize)
    if(is.null(windowThreshold)){
        cat("Constructing t-grid\n")
        windowThreshold  <- estimateWindowThreshold(nProbe,windowSize,...)
        print(windowThreshold)
    }

    if(length(windowSize) == 1){
        slidingWindow  <- oneWindowSizeScanner(observations, windowThreshold = windowThreshold, 
                             windowSize = windowSize)
    }else{
        if(!(length(windowSize)==length(windowThreshold)))
            stop("Error; windowSize and windowThreshold MUST be of equal lenght\n")

        slidingWindow  <- manyWindowSizeScanner(observations,windowThreshold = windowThreshold,
                                    windowSize = windowSize)
    }
    zhangWindow    <- slidingWindow[[1]]
    slidingValues  <- slidingWindow[[2]]
    zhangWhichK   <- slidingWindow[[3]]

    lowerBound     <- which(zhangWindow)[c(1,which(diff(which(zhangWindow))>1)+1)]
    upperBound     <- c(which(zhangWindow)[which(diff(which(zhangWindow))>1)],
                                                    rev(which(zhangWindow))[1])
    nregions        <- length(lowerBound)


    ### calculate Region-wise p-values
    if(nregions >= 1 & !anyNA(lowerBound)){
        regionLengths  <- do.call(c,lapply(getRegions(observations),length))
        regionIndex    <- rep(seq_along(regionLengths),regionLengths)[lowerBound]
        
        regions         <- data.frame(start=lowerBound,end=upperBound,length=upperBound-lowerBound+1,bump=regionIndex)
        index           <-    apply(regions, 1,function(x)seq(x[1],x[2]))
        if(is.matrix(index))index <- as.list(data.frame(index))

        kIndex     <- integer(max(windowSize));kIndex[windowSize] <- seq_len(length(windowSize))
        signRegion <- lapply(index,function(x,val,which.k){
                                   vv  <- val[,x]
                                   if(class(vv) == "numeric"){
                                       use <- !is.na(vv)
                                       out <- rbind(vv[use],which.k[x][use])
                                   }else{
                                       use <- !is.na(colSums(vv))
                                       out <- rbind(vv[,use],which.k[x][use])
                                   }
                                   return(out)
                            },
                                        val     = slidingValues,
                                        which.k = zhangWhichK)
        
        tVal       <- sapply(signRegion,function(x,kIndex){
                                      ll    <- ncol(x)
                                      win   <- nrow(x)
                                      j     <- 1
                                      t.new <- 0
                                      k.sum <- 0
                                      while(j <= ll){
                                       k       <- x[win,j]
                                       t.new   <- t.new + x[kIndex[k],j]*k
                                       k.sum   <- k.sum + k
                                       j       <- j + k 
                                    }
                                     return(c(t.new/k.sum,k.sum))
                                  },kIndex=kIndex)

   slidingValues.no.zero  <- stats::na.omit(as.numeric(slidingValues))
   ll                      <- length(slidingValues.no.zero)

   phi                     <- stats::acf(tVal[1,],plot = FALSE)$acf[2:3]
   p.val.empirical         <- apply(tVal,2,function(x,y,ll){(sum(abs(x[1]) <= y)+1)/ll},y = slidingValues.no.zero,ll = ll)
   sd                      <- sapply(tVal[2,],varAR, phi = phi, se = 1,order = 2)
   mean                    <- stats::median(apply(slidingValues, 1, mean, na.rm = TRUE))
   log.p.val.normal        <- -1*stats::pnorm(abs(tVal[1,]), mean=mean, sd =sd,lower.tail=FALSE,log.p=TRUE)/log(10)
   
   position                <- pos(observations)
   tVal.orig               <- tVal(observations)
   CpGnames                <- names(observations)
   
   ## Roll back regions into object region

        signRegions          <- RegionList(nRegions = nregions)

        for(i in seq_along(signRegion)){
            motherRegion        <- observations[[regionIndex[i]]]
            signIndex           <- index[[i]]
            signRegions         <- setRegion(signRegions,
                                             i      = i,
                                        Region(
                                          tValues   = tVal.orig[signIndex],
                                          position  = position[signIndex],
                                          chromosome= motherRegion@chromosome,
                                          pVal      = p.val.empirical[i],
                                          id        = CpGnames[signIndex]
                                        )
                              )
        }
        if(nregions > 1){
            signRegions <- as.GRanges(sort(signRegions))
        }


    }else{
        signRegions <- as.GRanges(RegionList(0L))
    }
          return(signRegions)
}

