#' @importFrom stats rnorm 
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm dmvnorm
importantSampling <- function(thresholdGrid,windowSize,nProbe,mcmc = 10000){
 ## Writen by LV   
    x   <-  stats::rnorm(nProbe)
    z   <-  numeric(nProbe)
    d   <- 2
    z[seq_len(d)] <-  stats::rnorm(d)

    for(i in seq_len(nProbe)[-1]){z[i]<-sum(x[max(i-d,1):min(i+d,nProbe)])/sqrt(2*d+1)}
    acf <- stats::acf(z,plot=FALSE)$acf

    lambda <- matrix(nrow = length(thresholdGrid), ncol= length(windowSize))

    for (iter in seq_along(windowSize)){
        k   <- windowSize[iter]
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

        Y <- MASS::mvrnorm(n = mcmc, mu=rep(0,c+k-1),Sigma=diag(c+k-1))
        #x <- mvrnorm(n=n,mu=m,Sigma=Sigma)
        x <- Y%*%A + m ## kanske t(A)?

## To account for abs() of the test statistic in the R and S - test
        x   <- abs(x) 
        m   <- abs(m)

        x2 = matrix(NA, nrow = mcmc, ncol = c)
      for (i in 1:c)
         x2[,i] = apply(x[,i:(i+k-1)],1, function(x,mcmc)sum(x)/mcmc,mcmc = k)

    for (j in seq_along(thresholdGrid)){
      w = 0
      if (c>2) {x3 <- apply((x2[,seq_len(c-1)]>= thresholdGrid[j]), 1,sum)*(x2[,c] >= thresholdGrid[j])}
      if (c==2){x3 <-(x2[,1]>= thresholdGrid[j])*(x2[,c] >= thresholdGrid[j])}

          for (i in which(x3==1) ) {
              tmp <- (nProbe-k+1)*mvtnorm::dmvnorm(x[i,], mean = rep(0,c+k-1),sigma = Sigma)/dmvnorm(x[i,], mean =m,sigma = Sigma)
              if(!is.nan(tmp)){
                  w = w + tmp
              }
          }
      lambda[j,iter] = w/mcmc
   }
    }
    return(lambda)
}

