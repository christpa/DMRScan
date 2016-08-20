.test <- function(x,i,k,t){

        xx      <- x[(i-k+1):i]
        xx      <- xx[!is.na(xx) & !is.infinite(xx)]
        ll      <- length(xx)
        if(ll == k)
            stat    <- sum(xx)/ll
        else
            stat    <- 0

            out     <- abs(stat) > t
    return(list(test=out,value=stat))
}

.Xt <- function(dat,t,k){

    nn          <- rownames(dat)
#    dat         <- abs(dat)
    L           <- length(dat)
    sign_probe  <- logical(L)
    value_probe <- numeric(L)

    if(L<=k){
        t           <- ifelse(!mean,yes=t,no=t*length(dat))
        value_probe <- rep(sum(dat,na.rm=TRUE),L)
        sign_probe  <- value_probe > t
    }else{
        for(i in (k+1):L){
            tmp <- .test(x=dat,i=i,k=k,t=t)
             if(tmp$test){
                 sign_probe[(i-k):i] <- TRUE
            }
            value_probe[i - k] <- tmp$value
        }
    }
        names(sign_probe)  <- nn

        out         <- list(sign_probe=sign_probe,value_probe=value_probe)
        class(out)  <- "scan"

       return(out)
}

.Xt.list <- function(dat,t,k){

    sign_probe  <- matrix(NA,sum(sapply(dat,nrow)),length(k))
    value_probe <- sign_probe
    
        for(i in seq_along(k)){
            tmp     <- sapply(dat,.Xt,t=t[i],k=k[i])
            sign_probe[,i]  <- do.call(c,lapply(tmp,function(x)x[[1]]))
            value_probe[,i] <- do.call(c,lapply(tmp,function(x)x[[2]]))
        }
    
        out <- list(sign_probe=sign_probe,value_probe=value_probe)
        

    return(out)
}

.Rt <- function(dat,t,k){

    nn          <- rownames(dat)
#    dat         <- abs(dat)
    L           <- length(dat)
    sign_probe  <- logical(L)
    value_probe <- numeric(L)

    if(L<=k){
        t           <- ifelse(!mean,yes=t,no=t*length(dat))
        sign_probe  <- rep(sum(dat,na.rm=TRUE)>t,L)
        value_probe[1] <- sum(abs(dat))
    }else{
        for(i in (k+1):L){
            tmp <- .test(x=dat,i=i,k=k,t=t)
            if(tmp$test & !any(sign_probe[(i-k+1):(i-1)])){
                sign_probe[(i-k+1):i] <- TRUE
            }
            value_probe[i - k] <- tmp$value
        }
    }
        names(sign_probe)   <- nn

        out         <- list(sign_probe=sign_probe,value_probe=value_probe)
        class(out)  <- "scan"

       return(out)
}

.Rt.list <- function(dat,t,k){
        
    tmp         <- sapply(dat,.Rt,t=t,k=k)
      
    sign_probe  <- t(as.matrix(do.call(c,tmp[1,])))
    value_probe <- t(as.matrix(do.call(c,tmp[2,])))
    
    which_k     <- integer(length(sign_probe))
    which_k[sign_probe] <- k

    out <- list(sign_probe=sign_probe,value_probe,which_k=which_k)
    
    return(out) 
}

.St <- function(dat,t_grid,k_grid){

    dat         <- abs(dat)
    nn          <- rownames(dat)
    L           <- length(dat)
    sign_probe  <- logical(L)
    value_probe <- matrix(0,length(k_grid),L)
    which_k     <- integer(L)

    for(runner in seq_along(k_grid)){
        k   <- k_grid[runner]
        for(i in k:L){
            if(k<L){
                tmp <- .test(x=dat,i=i,k=k,t=t_grid[runner])
                    if(tmp$test & !any(sign_probe[(i-k):i])){
                        sign_probe[(i-k+1):i] <- TRUE
                        which_k[i-k+1] <- k
                }
                value_probe[runner,i - k+1] <-  tmp$value
            }
        }       
    }               

    names(sign_probe)       <- nn
    rownames(value_probe)   <- k_grid
                
    out        <- list(sign_probe=sign_probe,value_probe=value_probe,which_k=which_k)

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
