make_data.arima <- function(L,...){

    x.out   <- arima.sim(n=L,...)
    
    return(as.numeric(x.out))

}   
