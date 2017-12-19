#' Estimate window thresholds
#' @name estimateThreshold
#' @description Estimate window thresholds for sliding window, 
#' one unique value for each window size
#'
#' @param nProbe The number of probes (CpGs) in the study. 
#' @param windowSize The different window sizes to be tested. 
#' Must be either one, or an ordered sequence of integers.
#' @param method Gives the method by which the threshold is calculated. 
#' Can be either an analytical solution "siegmund", provided by 
#' Siegnumd et.al (2012), or an iterative process; either importance 
#' sampling "sampling", as suggested by Zhang (2012) or a full MCMC model 
#' "mcmc" which can account for any dependency structure,  wich is pass 
#' to arima.sim, with ...
#' @param mcmc The number of MCMC iterations to be used, when using either 
#' Important Sampling ("zhang") or MCMC estimation of the threshold.
#' @param nCPU When calculating the thresholds on a cluster, how many CPUs 
#' should be used. This option is only compatible with the 'mcmc' method.
#' @param submethod A character string indicating if an AR(5) or ARIMA model 
#' should be used. In the AR(5), the index runs from -2 to 2. A regular AR(p) 
#' model can be obtaine using ARIMA(p,0,0) instead. 
#' @param ... Optinal parameters pased on to \code{\link[stats]{arima}}, when simulating data 
#' using the mcmc option, see arima.sim()
#' @return Returns a vector of the threshold for each window size
#' @examples
#' thresholdGrid <- estimateWindowThreshold(nProbe = 1000, 
#'                                     windowSize = 3:8, method = "siegmund")
#' @export
estimateWindowThreshold <- function(nProbe, windowSize, method = "siegmund",
                                mcmc = 1000, nCPU = 1,submethod = "ar",...){
        method              <- match.arg(method,c("siegmund","mcmc","sampling"))
        if(method == "siegmund"){
                estimateThreshold   <- estimateThresholdGrid.siegmund
        }else if(method == "mcmc"){
              estimateThreshold     <- estimateThresholdGrid.mcmc
        }else{
            estimateThreshold     <- estimateThresholdGrid.sampling
        }
 thresholdGrid      <- estimateThreshold(nProbe = nProbe, 
 windowSize = windowSize, mcmc = mcmc, nCPU = nCPU,submethod = submethod,...)
        return(thresholdGrid)
}

#' @importFrom parallel makeCluster parLapplyLB stopCluster
#' @importFrom stats smooth.spline predict
estimateThresholdGrid.mcmc <- function(nProbe, windowSize, mcmc, nCPU = 1,submethod = "ar",...){

    alpha       <- 0.05
    lambdaStar <- log((-log(1-alpha)/sum(windowSize))*nProbe)
    windowSize  <- sort(windowSize)
    thresholdGrid      <- seq(from=0.5,to=9,length=35)
    ## first row is k, second row is t
    iterGrid   <- as.data.frame(rbind(rep(windowSize,
                    each = length(thresholdGrid)),
                    rep(thresholdGrid,by=length(windowSize))))
    
    iterGrid   <- do.call(cbind,replicate(mcmc,iterGrid))
    
    if(nCPU == 1){
        cat("For parallel use, set nCPU > 1\n")
        significantWindow <- apply(iterGrid,2,function(x,nProbe,method,...) 
                            windowTest(threshold = x[2],windowSize = x[1],
                            nProbe = nProbe,method = submethod,...),
                            nProbe = nProbe,method = submethod,...)
    }else{
      iterGrid <- lapply(seq_len(ncol(iterGrid)),function(i,x,nProbe,submethod)
                          list(threshold = x[2,i],windowSize = x[1,i],
                               nProbe = nProbe,submethod = submethod),
                               x = iterGrid, nProbe, submethod)
#        options(MulticoreParam=quote(MulticoreParam(workers=nCPU))) 


       FUN <- function(x,...)windowTest(threshold  = x$threshold,
                                     windowSize = x$windowSize,
                                     nProbe     = x$nProbe,
                                     method     = x$submethod,...)

        cl <- parallel::makeCluster(nCPU,"PSOCK")
        significantWindow <- unlist(parallel::parLapplyLB(cl, iterGrid,FUN,...))
        parallel::stopCluster(cl)
    }

    resultArray  <- array(significantWindow,dim = c(length(thresholdGrid),
                                                length(windowSize),mcmc))
    
    result      <- matrix(0,length(thresholdGrid),length(windowSize))
    for(i in seq_len(nrow(result)))
        for(j in seq_len(ncol(result)))
            result[i,j] <- mean(resultArray[i,j,])
    
    logResult           <- log(t(result))
    rownames(logResult) <- windowSize
    colnames(logResult) <- thresholdGrid
    
    thresholdGridNew <- numeric(length(windowSize))
    if(any(result == 0)){ 
        for (j in seq_along(windowSize)) {
            m = max(which(logResult[j,]> lambdaStar))
            thresholdGridNew[j] = (log(lambdaStar)- logResult[j,m])/(logResult[j,m+1] - logResult[j,m])*(thresholdGrid[m+1]-thresholdGrid[m]) + thresholdGrid[m]
        }
    }else{
        resultSmooth <- apply(logResult,1,function(x,y,df){
                                        return(stats::smooth.spline(x=x,y=y,df=df))},
                                                        y = thresholdGrid,df=2)
        for(i in seq_along(windowSize))
            thresholdGridNew[i]   <- predict(resultSmooth[[i]],
                                            x=log(windowSize[i]*lambdaStar))$y
    }       
            
    return(thresholdGridNew)

}       

#' @importFrom stats smooth.spline predict
estimateThresholdGrid.siegmund <- function(windowSize, nProbe, 
                                                    mcmc = 1, nCPU = 1,...){
    
    alpha              <- 0.05
    thresholdGrid      <- seq(from=.5,to=8,length=100)
    result             <- matrix(NA,length(windowSize),length(thresholdGrid))
    lambdaStar         <- log((-log(1-alpha)/sum(windowSize))*nProbe)
  
    for(i in seq_along(windowSize)){
        for(j in seq_along(thresholdGrid)){
            result[i,j] <- siegmundLambda(threshold = thresholdGrid[j], windowSize = windowSize[i],nProbe = nProbe)
        }
    }
        thresholdGridNew <- numeric(length(windowSize))
        resultSmooth     <- apply(log(result),1,stats::smooth.spline,y=thresholdGrid,df=4)
       
  
  for(i in seq_along(windowSize))
    thresholdGridNew[i]   <- predict(resultSmooth[[i]],x=log(windowSize[i]*lambdaStar))$y

    return(thresholdGridNew)
}

estimateThresholdGrid.sampling <- function(windowSize, nProbe, mcmc = 1000, nCPU = nCPU,...){
  
    alpha       <- 0.05
    lambdaStar  <- log((-log(1-alpha)/sum(windowSize))*nProbe)
    windowSize  <- sort(windowSize)
    thresholdGrid  <- seq(from = 1,to = 4.5,.10)
    
    
    logResult  <- log(t(importantSampling(thresholdGrid = thresholdGrid,windowSize = windowSize,nProbe = nProbe, mcmc = mcmc)))
    thresholdGridNew  <- numeric(length(windowSize))

    rownames(logResult) <- windowSize
    colnames(logResult) <- thresholdGrid

    useSpline  <- any(apply(logResult,1,function(x)length(unique(x))<5))

    if(useSpline){
        for (j in seq_along(windowSize)) {
            m = max(which(logResult[j,]> lambdaStar))
            thresholdGridNew[j] = (log(lambdaStar)-logResult[j,m])/(logResult[j,m+1]-logResult[j,m])*(thresholdGrid[m+1]-thresholdGrid[m]) + thresholdGrid[m]
        }
    }else{
        resultSmooth <- apply(logResult,1,function(x,y,df){u <- is.finite(x);return(stats::smooth.spline(x=x[u],y=y[u],df=df))},y=thresholdGrid,df=3)
        for(i in seq_along(windowSize))
            thresholdGridNew[i]   <- predict(resultSmooth[[i]],x=log(windowSize[i]*lambdaStar))$y
    }

    return(thresholdGridNew)
}
