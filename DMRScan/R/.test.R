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
