.make.data.ou   <- function(theta,mu=0, L,sigma=1, burn.in=300,grid=0.01){

    L.new   <- L*1/grid + burn.in
    dt      <- grid
    dW      <- rnorm(L.new,sd=sigma*sqrt((1-exp(-2*theta*grid))/(2*theta)))
    x       <- numeric(L.new)
    elt     <- exp(-theta*grid)

    for(t in 2:L.new){
        x[t] <- x[t-1]*elt + (1-elt)*mu + dW[t]
    }

    x   <- x[seq(from=burn.in,to=L.new,length=L)]

    return(x)

}

