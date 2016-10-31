#' @import RcppRoll stats 
#' @importFrom Matrix sparseMatrix
windowTest <- function(threshold,windowSize,nProbe,method,...){

    makeData.ar <- function(nProbe,d = 2,...){
        ndiag <- function(par,dim){
        x <- rep(par[1],dim)
        i <- seq_len(dim)
        j <- seq_len(dim)
            if(length(par)>1 & dim > 2)
               for(k in seq_along(par)[-1]){
                    x <- c(x,rep(par[k],dim-k+1))
                    i <- c(i,k:dim)
                    j <- c(j,1:(dim -k+1))
              }
            mat <- Matrix::sparseMatrix(i=i,j=j,x=x,dims=c(dim,dim),symmetric=FALSE)
            return(mat)
        }       

    ll      <- nProbe + 2*d + 1
    x       <- rnorm(n = ll)
    x.new   <- x%*%ndiag(par = rep(1,d+1), dim = ll)

    x.out   <- x.new[seq(d+1,nProbe+d)]/sqrt(2*d+1)
    return(x.out)
        
    }   
    makeData.arima <- function(L,...){
        x.out   <- arima.sim(n=L,...)
        return(as.numeric(x.out))
    }

    submethod   <- get(paste("makeData",method,sep="."))
    x           <- submethod(nProbe,...)
    windows     <- roll_mean(abs(x),windowSize)
    sign        <- sum(roll_sum(windows > threshold,windowSize) == 1)

    return(sign)
}
