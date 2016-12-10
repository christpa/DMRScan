#' @importFrom stats pnorm dnorm
siegmundLambda <- function(threshold,windowSize,nProbe){
    m   <- nProbe
    w   <- windowSize
    z   <- threshold
    y   <- z*sqrt(2/w)/2
    lambda <- m*z*(w^-1)* stats::dnorm(z)*(1/y)*(stats::pnorm(y) - 
                                    0.5)/(y*stats::pnorm(y) + stats::dnorm(y))
    return(lambda)
}
