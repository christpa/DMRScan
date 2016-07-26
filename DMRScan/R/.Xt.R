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

