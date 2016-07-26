estimate_t_grid.zhang <- function(L,L.star=NULL,k_grid,mcmc=1000){
  
    if(is.null(L.star))L.star <- L

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

    return(list(t_grid.new=t_grid.new,res=l_res,t_grid.orig=t_grid))
}
