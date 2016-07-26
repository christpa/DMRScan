mean_AR <- function(p=1,K,phi,c){
    mean_AR_1 <- function(K,phi,c){
        const   <- c/K
        k       <- K-1
        tmp     <- sum((K-0:k)*phi^(0:k)) + sum(phi^(0:k))/(1-phi)
        return(tmp*const)
    }       
    mean_AR_2 <- function(K,phi,c){
        const   <- c/K
        k       <- K -1
        a.j     <- sapply(0:K,.a_j,phi=phi)
        mu      <- 1/(1-sum(phi))
        tmp     <- sum(a.j[-1]*mu + phi[1]*a.j[-K]*mu)
        return(tmp*const)
    }   

    if(p == 1)
        return(mean_AR_1(K,phi,c))
    else if(p == 2)
        return(mean_AR_2(K,phi,c))
    else
        return(NA)
}
