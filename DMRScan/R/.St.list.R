.St.list <- function(dat,t,k){

    ll          <- sum(sapply(dat,length))

    tmp         <- sapply(dat,.St,t=t,k=k)

    sign_probe  <- do.call(c,tmp[1,])
    value_probe <- do.call(cbind,tmp[2,])
    which_k     <- do.call(c,tmp[3,])

    out     <- list(sign_probe,value_probe,which_k=which_k)

    return(out)
}

