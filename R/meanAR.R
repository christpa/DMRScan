meanAR <- function(order = 1,nProbe, phi, c){
    
    meanAR1 <- function(nProbe, phi, c){
        const   <- c/nProbe
        k       <- nProbe-1
        tmp     <- sum((nProbe-(0:k))*phi^(0:k)) + sum(phi^(0:k))/(1-phi)
        return(tmp*const)
    }       
    meanAR2 <- function(nProbe, phi, c){
        const   <- c/nProbe
        k       <- nProbe -1
        a.j     <- sapply(0:k,.a_j,phi = phi)
        mu      <- 1/(1-sum(phi))
        tmp     <- sum(a.j[-1]*mu + phi[1]*a.j[-nProbe]*mu)
        return(tmp*const)
    }   

    if(order == 1)
        return(meanAR1(nProbe,phi,c))
    else if(order == 2)
        return(meanAR2(nProbe,phi,c))
    else
        return(NA)
}
