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
