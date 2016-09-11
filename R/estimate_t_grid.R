#' Estimate window thresholds for sliding window, one unique value for each window size
#'
#' @param L The number of probes (CpGs) in the study. 
#' @param k_grid The different window sizes to be tested. Must be either one, or an ordered sequence of integers.
#' @param method Gives the method by which the threshold is calculated. Can be either an analytical solution 
#'  "siegmund", provided by Siegnumd et.al (2012), or an iterative process; either important sampling "zhang", as suggested 
#'  by Zhang (2012) or a full MCMC model "mcmc" which can account for any dependency structure, wich is pass to arima.sim, with ...
#' @param mcmc The number of MCMC iterations to be used, when using either Important Sampling ("zhang") or MCMC estimation of the threshold.
#' @param n_cpu When calculating the thresholds on a cluster, how many CPUs should be used.
#' @param ... Optinal parameters pased on to arima(), when simulating data using the mcmc option, see arima.sim()
#' @return Returns a vector of the threshold for each window size
#' @importFrom stats smooth.spline predict
#' @examples
#' t_grid <- estimate_t_grid(L = 1000, k_grid = 3:8, method = "siegmund")
#' @export
#'
#'
estimate_t_grid <- function(L, k_grid, method = "siegmund",mcmc = 1000, n_cpu = 1,...){

        method      <- match.arg(method,c("siegmund","mcmc","zhang"))
        estimate_t  <- get(paste("estimate_t_grid",method,sep="."))
        
        t_grid      <- estimate_t(L = L, k_grid = k_grid, mcmc = mcmc, n_cpu = n_cpu,...)

        return(t_grid)

}


estimate_t_grid.mcmc <- function(L, k_grid, mcmc, n_cpu = 1,foo = "ar",...){

    alpha       <- 0.05
    lambda.star <- log((-log(1-alpha)/sum(k_grid))*L)
    k_grid      <- sort(k_grid)
    t_grid      <- seq(from=1,to=4,length=13)
    t_grid      <- seq(from=0.5,to=9,length=35)
    ## first row is k, second row is t
    iter.grid   <- as.data.frame(rbind(rep(k_grid,each=length(t_grid)),rep(t_grid,by=length(k_grid))))
    iter.grid   <- do.call(cbind,replicate(mcmc,iter.grid))
    sign_window <- apply(iter.grid,2,function(x,L,foo,...)window_test(t=x[2],k=x[1],L=L,foo=foo,...),L=L,foo=foo,...)
    
    res.array  <- array(sign_window,dim=c(length(t_grid),length(k_grid),mcmc))
    
    res        <- matrix(0,length(t_grid),length(k_grid))
    for(i in 1:nrow(res))
        for(j in 1:ncol(res))
            res[i,j] <- mean(res.array[i,j,])
    
    l_res       <- log(t(res))
    rownames(l_res) <- k_grid
    colnames(l_res) <- t_grid
    
    t_grid.new <- numeric(length(k_grid))
    if(any(res==0)){ 
        for (j in seq_along(k_grid)) {
            m = max(which(l_res[j,]> lambda.star))
            t_grid.new[j] = (log(lambda.star)-l_res[j,m])/(l_res[j,m+1]-l_res[j,m])*(t_grid[m+1]-t_grid[m]) + t_grid[m]
        }
    }else{
        res.smooth <- apply(l_res,1,function(x,y,df){return(smooth.spline(x=x,y=y,df=df))},y=t_grid,df=2)
        for(i in seq_along(k_grid))
            t_grid.new[i]   <- predict(res.smooth[[i]],x=log(k_grid[i]*lambda.star))$y
    }       
            
    return(t_grid.new)

}       

estimate_t_grid.siegmund <- function(k_grid, L, mcmc = 1, n_cpu = 1){
    
    alpha       <- 0.05
    t_grid      <- seq(from=.5,to=8,length=100)
    res         <- matrix(NA,length(k_grid),length(t_grid))
    lambda.star <- log((-log(1-alpha)/sum(k_grid))*L)


    
    for(i in seq_along(k_grid)){
        for(j in seq_along(t_grid)){
            res[i,j] <- .siegmund_lambda(t=t_grid[j],k=k_grid[i],L=L)
        }
    }
        t_grid.new <- numeric(length(k_grid))
        res.smooth <- apply(log(res),1,smooth.spline,y=t_grid,df=4)
       
        
        for(i in seq_along(k_grid))
            t_grid.new[i]   <- predict(res.smooth[[i]],x=log(k_grid[i]*lambda.star))$y

    return(t_grid.new)
}

estimate_t_grid.zhang <- function(L, k_grid, mcmc = 1000, n_cpu = n_cpu){
  
    alpha       <- 0.05
    lambda.star <- log((-log(1-alpha)/sum(k_grid))*L)
    k_grid      <- sort(k_grid)
    t_grid      <- seq(from=1,to=3.5,.10)
    
    
    l_res         <- log(t(.important_sampling(t_vec=t_grid,k_vec=k_grid,L=L,n=mcmc)))
    t_grid.new  <- numeric(length(k_grid))

    rownames(l_res) <- k_grid
    colnames(l_res) <- t_grid

    use.spline  <- !any(apply(l_res,1,function(x)length(unique(x))<5))

    if(!use.spline){
        for (j in seq_along(k_grid)) {
            m = max(which(l_res[j,]> lambda.star))
            t_grid.new[j] = (log(lambda.star)-l_res[j,m])/(l_res[j,m+1]-l_res[j,m])*(t_grid[m+1]-t_grid[m]) + t_grid[m]
        }
    }else{
        res.smooth <- apply(l_res,1,function(x,y,df){u <- is.finite(x);return(smooth.spline(x=x[u],y=y[u],df=df))},y=t_grid,df=3)
        for(i in seq_along(k_grid))
            t_grid.new[i]   <- predict(res.smooth[[i]],x=log(k_grid[i]*lambda.star))$y
    }

    return(t_grid.new)
}
