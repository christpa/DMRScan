estimate_t_grid.mcmc.ar <- function(L,k_grid,mcmc,n_cpu=1,foo="ar",...){

    alpha       <- 0.05
    lambda.star <- log((-log(1-alpha)/sum(k_grid))*L)
    k_grid      <- sort(k_grid)
    t_grid      <- seq(from=1,to=4,length=13)
    t_grid      <- seq(from=0.5,to=9,length=35)
    ## first row is k, second row is t
    iter.grid   <- as.data.frame(rbind(rep(k_grid,each=length(t_grid)),rep(t_grid,by=length(k_grid))))
    iter.grid   <- do.call(cbind,replicate(mcmc,iter.grid))
    sign_window <- pbapply(iter.grid,2,function(x,L,foo,...)window_test(t=x[2],k=x[1],L=L,foo=foo,...),L=L,foo=foo,...)
    
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
            
    return(list(t_grid.new=t_grid.new,res=l_res,t_grid.orig=t_grid))

}       

