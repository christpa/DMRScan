#' @rdname Rt-methods
#' @aliases Rt,GRangesList-method
setMethod("oneWindowSizeScanner", "GRangesList", function(region,windowThreshold,windowSize){
    
    slidingWindow  <- sapply(regions, oneWindowSizeScanner, windowThreshold = windowThreshold,
                          windowSize = windowSize)
      
    signProbe  <- t(as.matrix(do.call(c,slidingWindow[1,])))
    valueProbe <- t(as.matrix(do.call(c,slidingWindow[2,])))
    
    whichK     <- integer(length(signProbe))
    whichK[signProbe] <- windowSize 

    out <- list(signProbe=signProbe,valueProbe,whichK=whichK)
    
        return(out) 
    }
)


#' @rdname Rt-methods
#' @aliases Rt,GRanges-method
#' @importFrom RcppRoll roll_mean roll_sum
setMethod("oneWindowSizeScanner", "GRanges", function(region,windowThreshold,windowSize){

## Assumes that the values are ordered!! 
    dat         <- ncols(region)@listData$tVal
    nProbe      <- length(region)

    if(nProbe <= windowSize){
         signProbe <- rep(FALSE,nProbe)
         valueProbe <- rep(NaN, nProbe)
    }else{
        signProbe  <- logical(nProbe)
        window.observations <- roll_mean(dat,windowSize)
        valueProbe <- c(window.observations, rep(NaN,windowSize-1)) 
        nWindows    <- length(window.observations)
        sign_window <- window.observations > windowThreshold
        if(any(sign_window)){
            if(nProbe < 2*(windowSize -1)){
                for(i in seq_along(window.observations))
                    if(sign_window[i]){
## Significant windows, if first window is significant, TRUE, else 
                        if(i == 1){
                            signProbe[1:windowSize] <- TRUE
                        }else if(i == nWindows){
## Last window, Last window will overwrite NaNs
                            signProbe[(nProbe - windowSize):nProbe] <- TRUE
## chek for overalpp
                        }else if(window.observations[i-1] < windowThreshold){
                            signProbe[i:(i+windowSize -1)] <- TRUE
                        }
                    }
                }else{
                    overlapping_significant_windows <- roll_sum(sign_window, windowSize)
                if(any(overlapping_significant_windows > 1)){
## Overlapping significant windows
                 sign_window[which(overlapping_significant_windows > 1)] <- FALSE
                }
## Overlapping significant windows are removed 
                    which.sign <- do.call(c,lapply(which(sign_window),
                            function(x,windowSize){x:(x + windowSize -1)},
                                         windowSize = windowSize))
                    

                    signProbe[which.sign] <- TRUE
                }
        }
    }
    
        names(signProbe)   <- rownames(dat)
        out         <- list(signProbe=signProbe,valueProbe=valueProbe)

       return(out)
    }
)


#' @rdname St-methods                
setMethod("manyWindowSizeScanner", "GRangesList", function(region,windowThreshold,windowSize){ 

    slidingWindow  <- sapply(regions,manyWindowSizeScanner, windowThreshold = windowThreshold, windowSize = windowSize)

    signProbe  <- do.call(c,slidingWindow[1,])
    valueProbe <- do.call(cbind,slidingWindow[2,])
    whichK     <- do.call(c,slidingWindow[3,])

    out     <- list(signProbe,valueProbe,whichK=whichK)

        return(out)
    }   
)

#' @rdname St-methods             
#' @importFrom RcppRoll roll_mean 
setMethod("manyWindowSizeScanner", "GRanges", function(region,windowThreshold,windowSize){ 
                    

## Assumes that the values are ordered!! 
    dat         <- ncols(region)@listData$tVal
    nProbe      <- length(region)

    signProbe  <- logical(nProbe)
    valueProbe <- matrix(0,length(windowSize),nProbe)
    whichK     <- integer(nProbe)
    
    for(runner in seq_along(windowSize)){
        window   <- windowSize[runner]
## Identify any siginficant windows with windowSize "window"
## If region presented is shorted than windowSize; return NaN   
       if(nProbe > window){                                    
        tmp         <- roll_mean(dat, window)
        nWindows    <- length(tmp)
        sign_window <- (tmp > windowThreshold[runner])

        if(any(sign_window)){
## Significant windows
## Since smalles window comes first, any overlapping windows (ie. next 
## itteration will be non-sigificant.
           sign_window_list <- lapply(which(sign_window), 
                                    function(x,window){return(x:(x+window-1))},
                                    window = window)

        ## Special casese: last window significnat
            if(sign_window[nWindows] & !any(signProbe[(nProbe - window+1):nProbe])){
## Last window significant
                signProbe[(nProbe - window+1):nProbe] <- TRUE
            }

            for(i in seq_along(sign_window_list)){
                   if(!any(signProbe[sign_window_list[[i]]])){                       
                          signProbe[sign_window_list[[i]]] <- TRUE 
                          whichK[sign_window_list[[i]]]   <- length(sign_window_list[[i]])
                   }                                
            }
        }
            valueProbe[runner,] <- c(tmp, rep(NaN,window-1))
        }
    }
    
    names(signProbe)       <- rownames(dat)
    rownames(valueProbe)   <- windowSize 
                
    out <- list(signProbe=signProbe,valueProbe=valueProbe,whichK=whichK)

        return(out)
    }
)
