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

