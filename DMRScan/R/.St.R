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

