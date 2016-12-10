meanAR <- function(order = 1,nProbe, phi, c){
    if(order == 1){
        const   <- c/nProbe
        k       <- nProbe-1
        out     <- sum((nProbe-(0:k))*phi^(0:k)) + sum(phi^(0:k))/(1-phi)*const
    }else if(order == 2){
        const   <- c/nProbe
        k       <- nProbe -1
        a.j     <- sapply(0:k,.a_j,phi = phi)
        mu      <- 1/(1-sum(phi))
        out     <- sum(a.j[-1]*mu + phi[1]*a.j[-nProbe]*mu)*const
    }else{
        out     <- NA
    }

    return(out)
}
