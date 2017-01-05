varAR <- function(order,nProbe,phi,se){
   if(order == 1){
        out     <- se^2 * ( phi*(nProbe*phi - nProbe +2) - 2*phi^(1 + nProbe) + nProbe - phi*(nProbe+2))
        norm    <- (nProbe^2)* (1+phi)*(1-phi)^4
        out     <- out/norm
    }else if(order == 2){
        j   <- .gamma(phi)
        p.1 <- phi[1]
        p.2 <- phi[2]
        a.j <- sapply(seq_len(nProbe),.a_j,phi=phi)
        out   <- j[1]*(sum(a.j[-1])^2 + sum(a.j[-(nProbe+1)]*p.2 )^2) + 2*j[2]*(sum(a.j[-1])*sum(a.j[-(nProbe+1)]*p.2)) + sum(cumsum(a.j)^2)
        norm  <- nProbe^2
        out < out/norm*se^2
    }else{
       out <- NA
    }
    return(out)
}

.gamma <- function(phi){
    p.1 <- phi[1]
    p.2 <- phi[2]
    j.0 <- (1-p.2)/( (1+p.2)*(1-p.2)^2 - p.1^2 )
    j.1 <- j.0 * (p.1/(1-p.2) )
    return(c(j.0,j.1))
}   
 
.a_j <- function(j,phi){
     if(j  < 2){
         return( phi[1]^j)
     }else{
         return( phi[1] * .a_j(j-1,phi) +  phi[2]* .a_j(j-2,phi))
     }
} 
