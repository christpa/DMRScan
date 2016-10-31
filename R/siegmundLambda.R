#' @import stats
siegmundLambda <- function(threshold,windowSize,nProbe){
                    m   <- nProbe
                    w   <- windowSize
                    z   <- threshold

        nu <- function(x){
            y <- x/2
            out <- (1/y)*(pnorm(y)-0.5)/(y*pnorm(y) + dnorm(y))
            return(out) 
        }

        lambda   <- m*z*(w^-1)* dnorm(z)*nu(z*sqrt(2/w))

    return(lambda)
}
