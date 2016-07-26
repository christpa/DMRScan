make_data.ar <- function(L,...){
        
    ll      <- L+2*d+1
    x       <- rnorm(n = ll)
    x.new   <- x%*%n_diag(n = ll,d = d)

    x.out   <- x.new[(d+1):(L+d)]/sqrt(2*d+1)
    return(x.out)
        
}   
