#' DMR Scan function
#' 
#' @name dmrscan
#' @rdname DMRScan_slidingWindow
#' @param observations An object of either;  \code{\link{RegionList}} made by 
#' \code{\link{makeCpGregions}}, a vector of the test statistic, a \code{\link[GenomicRanges]{GRanges}} object, 
#' or a "minfi" object (soon to be supported). 
#' @param windowSize A sequence of windowSizes for the slidingWindow. Must be an 
#' integer vector, with equal length as the number of windows. 
#' @param windowThreshold Optional argument with corresponding cut-off for 
#' each window. Will be estimated if not supplied.
#' @param chr A vector of chromosomal position. Only used when the observations
#' vector is a matrix of test statistic.
#' @param pos A vector of genomic coordinates for the CpGs to match the chr argument
#' @param maxGap The maximum allowed gap between two CpGs within the same region. 
#' @param ... Optional arguments to be passed to \code{\link{estimateThreshold}}, 
#' if no grid is specified.
#' @return An object of type \code{\link[GenomicRanges]{GRanges}} with significantly differentially 
# methylated regions
#' @keywords DMRScan
#' @importFrom GenomicRanges GRanges
#' @importFrom stats acf median pnorm na.omit
#' @export
#' @examples
#' ## methylation data from chromosome 22 
#' data(DMRScan.methylationData)
#' ## phenotype (end-point for methylation data)
#' 
#' data(DMRScan.phenotypes) 
#' 
#' ## Test for an association between phenotype and methylation
#' test.statistics <- apply(DMRScan.methylationData,1,function(x,y)
#'   summary(glm(y ~ x, family = binomial(link = "logit")))$coefficients[2,3],
#'                                                    y = DMRScan.phenotypes)
#' ## Set chromosomal position to each test-statistic
#' positions <- data.frame(matrix(as.integer(unlist(strsplit(names(test.statistics), 
#'                                split="chr|[.]"))), ncol = 3, byrow = TRUE))[,-1]
#' ## Set clustering features 
#' min.cpg <- 4  ## Minimum number of CpGs in a tested cluster
#' ## Maximum distance (in base-pairs) within a cluster 
#' ## before it is broken up into two separate cluster 
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
#' nCpG        <- sum(sapply(regions,length)) ## Number of CpGs to be tested
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
dmrscan <- function(observations,windowSize,windowThreshold=NULL,chr = NULL, pos = NULL, maxGap = 500,...){

    if(class(observations) == "matrix"){
        observations    <- makeCpGregions(observations, chr = chr, pos = pos, maxGap = maxGap, minCpG = min(windowSize))
    }else if( class(observations) == "minfi"){
        ## Not correct class name
        print("Minfi objects are not supported yet")
        return(1)
    }
<<<<<<< HEAD

=======
>>>>>>> master
    alpha       <- 0.05

    windowSize  <- sort(windowSize)
    if(is.null(windowThreshold)){
    	nProbe      <- sum(sapply(observations,length))
        cat("Constructing t-grid\n")
    	nProbe      <- sum(sapply(observations,length))
        windowThreshold  <- estimateWindowThreshold(nProbe,windowSize,...)
        print(windowThreshold)
    }

    if(length(windowSize) == 1){
        slidingWindow  <- oneWindowSizeScanner(observations, windowThreshold = windowThreshold, 
                             windowSize = windowSize)
    }else{
        if(!(length(windowSize)==length(windowThreshold)))
            stop("Error; window size and window threshold must be of equal length\n")

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
        regionLengths  <- sapply(observations,length)
        regionIndex    <- rep(seq_along(regionLengths),regionLengths)[lowerBound]
        
        regions         <- data.frame(start=lowerBound,end=upperBound,length=upperBound-lowerBound+1,bump=regionIndex)
        index           <-    apply(regions, 1,function(x)seq(x[1],x[2]))
        if(is.matrix(index))index <- as.list(data.frame(index))

        kIndex     <- integer(max(windowSize));kIndex[windowSize] <- seq_len(length(windowSize))
        signRegion <- lapply(index,function(x,val,which.k){
                                   vv  <- val[,x]
                                   if(class(vv) == "numeric"){
                                       use <- !is.na(vv)
                                       fail    <- sum(which.k[x]) == 0 | is.null(use)
                                       if(fail)
                                           out <- NULL
                                        else
                                            out <- rbind(vv[use],which.k[x][use])
                                   }else{
                                       use <- !is.na(colSums(vv))
                                       if(length(use) == 0)
                                           vv <- na.omitt.dmrscan(vv)
                                       use <- !is.na(colSums(vv))
                                       fail    <- sum(which.k[x]) == 0 | length(use) == 0
                                       if(fail)       
                                           out <- NULL
                                       else
                                           out <- rbind(vv[,use],which.k[x][use])
                                   }
                                   return(out)
                            },
                                        val     = slidingValues,
                                        which.k = zhangWhichK)
        signRegion.noZero <- sapply(signRegion,length)>0
        nregions          <- sum(signRegion.noZero)
        signRegion        <- signRegion[signRegion.noZero]

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
#   sd                      <- sapply(tVal[2,],varAR, phi = phi, se = 1,order = 1)
#   mean                    <- stats::median(apply(slidingValues, 1, mean, na.rm = TRUE))
#   log.p.val.normal        <- -1*stats::pnorm(abs(tVal[1,]), mean=mean, sd =sd,lower.tail=FALSE,log.p=TRUE)/log(10)
   
   position                <- pos(observations)
   chr					   <- chr(observations)
   tVal.orig               <- tVal(observations)
   CpGnames                <- id(observations)


   ## Roll back regions into object region
    
        if(any(!signRegion.noZero)){
            regionIndex <- regionIndex[signRegion.noZero]
            index        <- index[signRegion.noZero]
<<<<<<< HEAD
        }
		signRegions	<- as.data.frame(matrix(NA, ncol = 6, nrow = length(regionIndex), dimnames = list(NULL, c("seqnames","start","end","no.cpg","pVal","tVal"))))
		for(i in seq_along(signRegion)){
            signIndex           <- index[[i]]
			signRegions[i,"seqnames"]	<- chr[signIndex[1]]
			signRegions[i,"start"]	<- min(position[signIndex])
			signRegions[i,"end"]	<- max(position[signIndex])
			signRegions[i,"pVal"]	<- p.val.empirical[i]
			signRegions[i,"tVal"]	<- tVal[1,i]
#			signRegions[i,"id"]		<- paste(CpGnames[signIndex], collapse="|")
			signRegions[i,"no.cpg"]	<- length(signIndex)

                                        
        }
			signRegions <- GRanges(signRegions) 
=======
        }else{
        ## Do I need a list for this?? GRanges is sufficient
		signRegions.df <- data.frame("start" = NA, "stop" = NA, seqname = "", "no.cpg" = NA, "tVal" = NA, "pVal" = NA, "id"= NA, row.names = 1:length(index)) 	
		for(i in seq_along(signRegion)){
            signIndex           <- index[[i]]
            signRegions.df$seqnames <- chr[signIndex][1]
			signRegions.df$start	<- min(position[signIndex])
			signRegions.df$end		<- max(position[signIndex])
			signRegions.df$no.cpg	<- length(signIndex)
		    signRegions.df$pval		<- p.val.empirical[i]
			signRegions.df$id		<- paste(names(CpGnames[signIndex]),collapse="|")
			
#			signRegions.df[i,] <- c(
#									seqnames = chr[signIndex],
#									ranges	= IRanges(start=position[signIndex], end = position[signIndex]),
#                                    no.cpg	= length(signIndex),
#									tVal	= tVal.orig[signIndex],
#                                    pVal    = p.val.empirical[i],
#                                    id      = CpGnames[signIndex]
#                                        )
        }
		}
		signRegions <- GRanges(signRegions.df) 
>>>>>>> master
    }else{
        signRegions <- GRanges()
    }
          return(signRegions)
}

