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
