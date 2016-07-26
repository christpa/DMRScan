n_diag    <- function(n,d){
            
    x <- rep(1,n)
    i <- 1:n
    j <- 1:n

     for(k in 2:(d+1)){
            x <- c(x,rep(1,n-k+1))
            i <- c(i,k:n)
            j <- c(j,1:(n -k+1))
        }
    
    mat <- sparseMatrix(i=i,j=j,x=x,dims=c(n,n),symmetric=TRUE)
    return(mat)
}           
