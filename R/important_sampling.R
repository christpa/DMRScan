#' @importFrom stats rnorm acf 
#' @import MASS
#' @importFrom MASS mvrnorm
#' @import mvtnorm
#' @importFrom mvtnorm dmvnorm
.important_sampling <- function(t_vec,k_vec,L,n=10000){
 ## Writen by LV   
    x   <-  rnorm(L)
    z   <-  numeric(L)
    d   <- 2
    z[1:d]<-  rnorm(d)
    for(i in d:L){z[i]<-sum(x[max(i-d,1):min(i+d,L)])/sqrt(2*d+1)}
    acf <- acf(z,plot=FALSE)$acf

    lamb <- matrix(nrow = length(t_vec), ncol= length(k_vec))

    for (iter in seq_along(k_vec)){
        k   <- k_vec[iter]
        c   <- max(2*d+1,k)
        l   <- 1
        m   <- rep(0,c+k-1)
        m[(c-min(2*d+1,k)+1):(c+min(2*d+1,k)-1)] <- c(acf[min(2*d+1,k):1], acf[2:min(2*d+1,k)])

        Sigma = matrix(0,c+k-1,c+k-1)

        for (i in seq_along(m)){
            for (j in seq_along(m)){
                if (abs(i-j)< (2*d +1))
                    Sigma[j,i] = acf[abs(i-j)+1]
            }
        }
        A <- chol(Sigma)

        Y <- mvrnorm(n=n,mu=rep(0,c+k-1),Sigma=diag(c+k-1))
        #x <- mvrnorm(n=n,mu=m,Sigma=Sigma)
        x <- Y%*%A + m ## kanske t(A)?

        x   <- abs(x) ## To account for abs() of the test statistic in the R and S - test
        m   <- abs(m)

        x2 = matrix(NA, nrow=n, ncol = c)
        for (i in 1:c)
            x2[,i] = apply(x[,i:(i+k-1)],1, function(x,n)sum(x)/n,n=k)

        for (j in seq_along(t_vec)){
            w = 0
            if (c>2) {x3 <- apply((x2[,1:(c-1)]>= t_vec[j]), 1,sum)*(x2[,c] >= t_vec[j])}
            if (c==2){x3 <-(x2[,1]>= t_vec[j])*(x2[,c] >= t_vec[j])}

                for (i in which(x3==1) ) {
                    tmp <- (L-k+1)*dmvnorm(x[i,], mean = rep(0,c+k-1),sigma = Sigma)/dmvnorm(x[i,], mean =m,sigma = Sigma)
                    if(!is.nan(tmp)){
                        w = w + tmp
                    }
                }
            lamb[j,iter] = w/n
        }
    }
    return(lamb)
}

