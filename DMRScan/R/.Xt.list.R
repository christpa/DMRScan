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

