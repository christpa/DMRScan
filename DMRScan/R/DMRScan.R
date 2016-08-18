#' DMR Scan function
#'
#' This function search for DMRs given a orded sequence of test statistics
#' @param obs [double_seq] Sequence of test statistics
#' @param k_grid [int_seq] Sequence of window sizes for the sliding window 
#' @param t_grid [double_seq] Optional argument with corresponding cut-off for each window. Will be estimated if not supplied
#' @keywords DMRScan
#' @import Rcpp
#' @export
#' @examples
#' 
#' N               <- 200
#' pos             <- cumsum(rpois(N,30))
#' chr             <- rep("chr1", N); chr.int <- rep(1,N) ## use plink annotation 1-22, X = 23, Y = 24, XY = 25, MT = 26
#' z_val           <- matrix(rnorm(N),dimnames=list(paste(chr,pos, sep = ".")))
#' 
#' min.cpg <- 3  # minimum number of cps for a region to ba chnsidered
#' max.gap <- 30 # Maximum allowd gap between two observations within one region
#' regions <- make.cpg.region.list(test_statistic=z_val,chr=chr.int,pos=pos,max.gap=max.gap,min.cpg=min.cpg)
#' 
#' nobs    <- sum(sapply(regions,nrow)) ## observations remaining
#' 
#' k_grid  <- c(2,3,4) ## window sizes 
#' mcmc    <- 1000 ## number of simulation for window sizes to be deterimed
#' 
#' t_grid  <- estimate_t_grid.zhang(k_grid=k_grid,L=L,mcmc=mcmc)
#' res     <- DMRScan(obs=regions,k_grid=k_grid,t_grid=t_grid$t_grid.new) ## If no regions are found, the function returns NA
#' 
DMRScan <- function(obs,k_grid,t_grid=NULL,...){

    xx          <- sapply(obs,length)
    L           <- sum(xx)
    alpha       <- 0.05
    lambda.star <- (-log(1-alpha)/sum(k_grid))*L

    k_grid  <- sort(k_grid)
    if(is.null(t_grid)){
        cat("Constructing t-grid\n")
        t_grid  <- estimate_t_grid.zang(L,k_grid)$t_grid.new
        print(t_grid)
    }

    if(length(k_grid) == 1){
        sliding.window  <- .Rt.list(obs,t=t_grid,k=k_grid)
    }else{
        if(!(length(k_grid)==length(t_grid)))
            stop("Error; k_grid and t_grid MUST be of equal lenght\n")

        sliding.window  <- .St.list(obs,t=t_grid,k=k_grid)
    }
    zhang.window    <- sliding.window[[1]]
    sliding.values  <- sliding.window[[2]]
    zhang.which.k   <- sliding.window[[3]]

    lower_bound     <- which(zhang.window)[c(1,which(diff(which(zhang.window))>1)+1)]
    upper_bound     <- c(which(zhang.window)[which(diff(which(zhang.window))>1)],rev(which(zhang.window))[1])
    nregions        <- length(lower_bound)
    region.index    <- rep(1:length(xx),xx)[lower_bound]


    ### calculate Region-wise p-values
    if(nregions > 1){

        regions         <- data.frame(start=lower_bound,end=upper_bound,length=upper_bound-lower_bound+1,bump=region.index)
        index       <-    apply(regions, 1,function(x)x[1]:x[2])
        if(is.matrix(index))index <- as.list(data.frame(index))

        k_index     <- integer(max(k_grid));k_index[k_grid] <- 1:length(k_grid)
        sign_region <- lapply(index,function(x,val,which.k)rbind(val[,x],which.k[x]),val=sliding.values,which.k=zhang.which.k)
        t.val       <- sapply(sign_region,function(x,k_index){
                                                        ll    <- ncol(x)
                                                        win   <- nrow(x)
                                                        j     <- 1
                                                        t.new <- 0
                                                        k.sum <- 0
                                                        while(j < ll){
                                                         k       <- x[win,j]
                                                         t.new   <- t.new + x[k_index[k],j]*k
                                                         k.sum   <- k.sum + k
                                                         j       <- j + k
                                                      }
                                                       return(c(t.new/k.sum,k.sum))
                                                    },k_index=k_index)

         sliding.values.no.zero  <- do.call(cbind,apply(sliding.values,2,function(x){if(!sum(x==0)>0)x}))
         ll                      <- length(sliding.values.no.zero)

    ##
        p.val.empirical         <- apply(t.val,2,function(x,y,ll){sum(abs(x[1]) <= y)/ll},y = sliding.values.no.zero,ll = ll)
        sd                      <- sapply(sapply(t.val[2,],min,16),var_AR_2,phi=c(0.2506930,0.1194123),s=0.8165831)
        #mean                    <- sapply(t.val[2,],mean_AR_1,phi=0.23,c=0.18) 
        mean                    <- median(apply(sliding.values.no.zero,1,mean))
        log.p.val.normal        <- -1*pnorm(abs(t.val[1,]), mean=mean, sd =sd,lower.tail=FALSE,log.p=TRUE)/log(10)
   ###

        if(is.null(rownames(obs[[1]]))){rr <- do.call(c,lapply(obs,names))}else{rr <- do.call(c,lapply(obs,rownames))}
        regions.coord   <- apply(regions,1,function(x,rr){c(rr[x[1]],rr[x[2]])},rr=rr)
        regions.coord   <- paste(regions.coord[1,],matrix(unlist(strsplit(regions.coord[2,],split="[.]")),ncol=2,byrow=TRUE)[,2],sep="-")
        regions.coord   <- gsub(regions.coord,pattern="[.]",replacement=":")

        regions         <- data.frame(regions,"p.val.emp"=p.val.empirical,"logP.normal"=log.p.val.normal,genome.coord=regions.coord,stringsAsFactors = FALSE)
        regions         <- regions[order(regions$p.val.emp),]
   }else{
       regions          <- NA
    }


    return(regions)
}

