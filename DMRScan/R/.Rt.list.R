.Rt.list <- function(dat,t,k){
        
    tmp         <- sapply(dat,.Rt,t=t,k=k)
      
    sign_probe  <- t(as.matrix(do.call(c,tmp[1,])))
    value_probe <- t(as.matrix(do.call(c,tmp[2,])))
    
    which_k     <- integer(length(sign_probe))
    which_k[sign_probe] <- k

    out <- list(sign_probe=sign_probe,value_probe,which_k=which_k)
    
    return(out) 
}

