estimate_t_grid.siegmund <- function(k_grid,L){
    
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

    return(list(t_grid.new=t_grid.new,res=res,t_grid.orig=t_grid))
}
