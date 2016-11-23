#' @importFrom stats pnorm dnorm
siegmundLambda <- function(threshold,windowSize,nProbe){
                    m   <- nProbe
                    w   <- windowSize
                    z   <- threshold

        nu <- function(x){
            y <- x/2
            out <- (1/y)*(stats::pnorm(y)-0.5)/(y*stats::pnorm(y) + stats::dnorm(y))
            return(out) 
        }

        lambda   <- m*z*(w^-1)* stats::dnorm(z)*nu(z*sqrt(2/w))

    return(lambda)
}
