window_test <- function(t,k,L,foo,...){

    foo     <- get(paste("make_data",foo,sep="."))
    x       <- foo(L,...)
    windows <- roll_mean(abs(x),k)
    sign    <- sum(roll_sum(windows > t,k) == 1)

    return(sign)
}

