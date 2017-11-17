#' @importFrom Matrix sparseMatrix
#' @importFrom RcppRoll roll_mean roll_sum
windowTest <- function(threshold,windowSize,nProbe,method,...){
<<<<<<< HEAD

=======
    
>>>>>>> upstream/master
    method <- match.arg(method,c("ar", "arima"))
    if(method == "arima"){
        submethod <- .makeData.arima
    }else if(method == "ar"){
        submethod <- .makeData.ar
    }
    x           <- submethod(nProbe,...)
    windows     <- RcppRoll::roll_mean(abs(x),windowSize)
    sign        <- sum(RcppRoll::roll_sum(windows > threshold,windowSize) == 1)
    return(sign)
}
<<<<<<< HEAD
=======

>>>>>>> upstream/master
.ndiag <- function(par,dim){
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
<<<<<<< HEAD
=======

>>>>>>> upstream/master
.makeData.ar <- function(nProbe,d = 2,...){
    ll      <- nProbe + 2*d + 1
    x       <- stats::rnorm(n = ll)
    x.new   <- x%*%.ndiag(par = rep(1,d+1), dim = ll)

    x.out   <- x.new[seq(d+1,nProbe+d)]/sqrt(2*d+1)
    return(x.out)
<<<<<<< HEAD
}
.makeData.arima <- function(L,...){
	        x.out   <- stats::arima.sim(n=L,...)
        return(as.numeric(x.out))
}
=======
}   

.makeData.arima <- function(L,...){
        x.out   <- stats::arima.sim(n=L,...)
        return(as.numeric(x.out))
}


>>>>>>> upstream/master
