#.test <- function(x,i,k,t){
#
#        xx      <- x[(i-k+1):i]
#        xx      <- xx[!is.na(xx) & !is.infinite(xx)]
#        ll      <- length(xx)
#        if(ll == k)
#            stat    <- sum(xx)/ll
#        else
#            stat    <- 0
#
#            out     <- abs(stat) > t
#    return(list(test=out,value=stat))
#}
#
#
#.Xt <- function(dat,window.threshold,window.size){
#
#    nn          <- rownames(dat)
#    dat         <- abs(dat)
#    L           <- length(dat) - window.size + 1
#
#    sign_probe  <- logical(L)
#    value_probe <- numeric(L)
#
#    if(L<=k){
# #       t           <- ifelse(!mean,yes=t,no=t*length(dat))
# #       value_probe <- rep(sum(dat,na.rm=TRUE),L)
# #       sign_probe  <- value_probe > t
#        # Return FALSE
#    }else{
#        window.observations <- roll_mean(dat, window.size)
#        
#        if(any(window.observations) > window.threshold){
#            ## significant windows
#
#
#        }
#        
#        
#        
#        
#        
#        for(i in (k+1):L){
#            tmp <- .test(x=dat,i=i,k=k,t=t)
#             if(tmp$test){
#                 sign_probe[(i-k):i] <- TRUE
#            }
#            value_probe[i - k] <- tmp$value
#        }
#    }
#        names(sign_probe)  <- nn
#
#        out         <- list(sign_probe=sign_probe,value_probe=value_probe)
#        class(out)  <- "scan"
#
#       return(out)
##}
#
#.Xt.list <- function(dat,t,k){
#
#    sign_probe  <- matrix(NA,sum(sapply(dat,nrow)),length(k))
#    value_probe <- sign_probe
#    
#        for(i in seq_along(k)){
#            tmp     <- sapply(dat,.Xt,t=t[i],k=k[i])
#            sign_probe[,i]  <- do.call(c,lapply(tmp,function(x)x[[1]]))
#            value_probe[,i] <- do.call(c,lapply(tmp,function(x)x[[2]]))
#        }
#    
#        out <- list(sign_probe=sign_probe,value_probe=value_probe)
#        
#
#    return(out)
#}

.Rt <- function(dat,window.threshold,window.size){

    dat         <- abs(dat)
    nProbe      <- length(dat)

    if(nProbe <= window.size){
         sign_probe <- rep(FALSE,nProbe)
         value_probe <- rep(NaN, nProbe)
    }else{
    sign_probe  <- logical(nProbe)
    value_probe <- c(window.observations, rep(NaN,window.size))
        
    window.observations <- roll_mean(dat,window.size)
    nWindows            <- length(window.observations)
        if(any(window.observations) > window.threshold){
            if(nProbe < 2*(window.size -1)){
                for(i in seq_along(window.observations))
                    if((window.observations[i] > window.threshold)){
## Significant windows
                    if(i == 1){
                        sign_probe[1:window.size] <- TRUE
                    }else if(window.observations[i-1] < window.threshold){
                        sign_probe[i:(i+window.size)] <- TRUE
                    }
            }else{
                overlapping_significant_windows <- roll_sum(window.observations, window.size)
                if(any(overlapping_significant_windows > 1)){
## Overlapping significant windows
                    window.observations[which(overlapping_significant_windows > 1) -1 ] <- 0
                }
## Overlapping significant windows are removed 
                    which.sign <- which(window.observations > window.threshold)
                    sign_probe[which.sign:(which.sign + window.size)] <- TRUE
        }
#        for(i in (k+1):L){
#            tmp <- .test(x=dat,i=i,k=k,t=t)
#            if(tmp$test & !any(sign_probe[(i-k+1):(i-1)])){
#                sign_probe[(i-k+1):i] <- TRUE
#            }
#            value_probe[i - k] <- tmp$value
#        }
    }
        names(sign_probe)   <- rownames(dat)
        out         <- list(sign_probe=sign_probe,value_probe=value_probe)
        class(out)  <- "scan"

       return(out)
}

.Rt.list <- function(dat,window.threshold,window.size){
        
    tmp         <- sapply(dat,.Rt, window.threshold = window.threshold,
                          window.size = window.size)
      
    sign_probe  <- t(as.matrix(do.call(c,tmp[1,])))
    value_probe <- t(as.matrix(do.call(c,tmp[2,])))
    
    which_k     <- integer(length(sign_probe))
    which_k[sign_probe] <- k

    out <- list(sign_probe=sign_probe,value_probe,which_k=which_k)
    
    return(out) 
}

.St <- function(dat,window.threshold,window.size){
## window.threshold & window.size are vectors/grid of values
## Assumes window sizes is sorted from smallest to largers

    dat         <- abs(dat)
    nProbe      <- length(dat)
    sign_probe  <- logical(nProbe)
    value_probe <- matrix(0,length(k_grid),nProbe)
    which_k     <- integer(nProbe)

    for(runner in seq_along(window.size)){
        window   <- window.size[runner]
## Identify any siginficant windows with window size "window"
        tmp     <- roll_mean(dat, window)
        
        sign_window <- (tmp > window.threshold[runner])

        if(any(sign_window)){
## Significant windows
## Since smalles window comes first, any overlapping windows (ie. next 
## itteration will be non-sigificant.
           sign_window_list <- lapply(which(sign_window), function(x,window){
                                    return(x:(x+window-1))}, window = window)

            for(i in seq_along(sign_window_list)){
                        if(!any(sign_probe[sign_window_list[[i]]])){                       
                                sign_probe[sign_window_list[[i]]] <- TRUE 
                                which_k[sign_window_list[[i]]]   <- length(sign_window_list[[i]])
                              }                                
            }
    }
                   
        
    names(sign_probe)       <- rownames(dat)
    rownames(value_probe)   <- k_grid
                
    out <- list(sign_probe=sign_probe,value_probe=value_probe,which_k=which_k)

    class(out) <- "scan"
    return(out)
}

.St.list <- function(dat,t,k){

    ll          <- sum(sapply(dat,length))

    tmp         <- sapply(dat,.St,t=t,k=k)

    sign_probe  <- do.call(c,tmp[1,])
    value_probe <- do.call(cbind,tmp[2,])
    which_k     <- do.call(c,tmp[3,])

    out     <- list(sign_probe,value_probe,which_k=which_k)

    return(out)
}   
