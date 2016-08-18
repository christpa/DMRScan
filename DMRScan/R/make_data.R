make_data.ar <- function(L,...){
        
    ll      <- L+2*d+1
    x       <- rnorm(n = ll)
    x.new   <- x%*%.ndiag(d = d, n = ll)

    x.out   <- x.new[(d+1):(L+d)]/sqrt(2*d+1)
    return(x.out)
        
}   
make_data.arima <- function(L,...){

    x.out   <- arima.sim(n=L,...)
    
    return(as.numeric(x.out))

}

.ndiag <- function(par,n){
    x <- rep(par[1],n)
    i <- 1:n
    j <- 1:n
        if(length(par)>1)
           for(k in 2:length(par)){
                x <- c(x,rep(par[k],n-k+1))
                i <- c(i,k:n)
                j <- c(j,1:(n -k+1))
          }
    mat <- Matrix::sparseMatrix(i=i,j=j,x=x,dims=c(n,n),symmetric=!TRUE)
    return(mat)
}


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


