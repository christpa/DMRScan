var_AR <- function(p,K,phi,s){
    
    var_AR_1 <- function(K,phi,s){
        out     <- s^2 * ( phi*(K*phi - K +2) - 2*phi^(1+K) + K - phi*(K+2))
        norm    <- (K^2)* (1+phi)*(1-phi)^4
        return(out/norm)
    }

    var_AR_2 <- function(K,phi,s){
        .gamma <- function(phi){
            p.1 <- phi[1]
            p.2 <- phi[2]

            j.0 <- (1-p.2)/( (1+p.2)*(1-p.2)^2 - p.1^2 )
            j.1 <- j.0 * (p.1/(1-p.2) )
            return(c(j.0,j.1))
        }   
  
    
    j   <- .gamma(phi)
    p.1 <- phi[1]
    p.2 <- phi[2]
    a.j <- sapply(0:K,.a_j,phi=phi)

    out   <- j[1]*(sum(a.j[-1])^2 + sum(a.j[-(K+1)]*p.2 )^2) + 2*j[2]*(sum(a.j[-1])*sum(a.j[-(K+1)]*p.2)) + sum(cumsum(a.j)^2)
    norm  <- K^2
    
    return(out/norm*s^2)
}   
  
    if(p == 1)
        return(var_AR_1(K,phi,s))
    else if(p == 2)
        return(var_AR_2(K,phi,s))
    else
        return(NA)
}

.a_j <- function(j,phi){
     if(j  < 2){
         return( phi[1]^j)
     }else{
         return( phi[1] * .a_j(j-1,phi) +  phi[2]* .a_j(j-2,phi))
     }
} 
