#' DMR Scan function
#'
#' This function search for DMRs given a orded sequence of test statistics
#' @param obs A sequence of test statistics for each CpG
#' @param k_grid A sequence of window sizes for the sliding window, must be an integer 
#' @param t_grid Optional argument with corresponding cut-off for each window. Will be estimated if not supplied.
#' @param ... Optional arguments to be pased to estimate_t_grid(), if no grid is specified.
#' @return A list of Differentially Methylated Regions (DMRs). Returns NA if no DMRs were found.
#' @keywords DMRScan
#' @importFrom stats median pnorm
#' @export
#' @examples
#' data(DMRScan.methylationData) ## Load methylation data from chromosome 22
#' data(DMRScan.phenotypes) ## Load phenotype (end-point for methylation data)
#' 
#' ## Test for an association between phenotype and Methylation
#' test.statistics <- apply(DMRScan.methylationData,1,function(x,y)summary(glm(y ~ x, family = binomial(link = "logit")))$coefficients[2,3], y = DMRScan.phenotypes)
#' ## Set chromosomal position to each test-statistic
#' pos<- data.frame(matrix(as.integer(unlist(strsplit(names(test.statistics), split="chr|[.]"))), ncol = 3, byrow = TRUE))[,-1] 
#' ## Set clustering features 
#' min.cpg <- 4  ## Minimum number of CpGs in a tested cluster
#' ## Maxium distance (in base-pairs) within a cluster 
#' ## before it is broken up into two seperate cluster 
#' max.gap <- 750  
#' 
#' window.sizes <- 3:7 ## Number of CpGs in the sliding windows, can be either a single number or a sequence of window sizes
#' n.CpG        <- nrow(DMRScan.methylationData) ## Number of CpGs to be tested
#' 
#' ## Estimate the window threshold, based on the number of CpGs and window sizes
#' window.thresholds <- estimate_t_grid(k_grid = window.sizes, L = n.CpG, method = "zhang", mcmc = 10000)
#' 
#' ## Identify all clusters, and generate a list for each cluster
#' regions <- make.cpg.regions(test_statistic = test.statistics, chr = pos[,1], pos = pos[,2], max.gap = max.gap, min.cpg = min.cpg)
#' ## Run the sliding window
#' dmrscan.results   <- dmr_scan(obs = regions, k_grid = window.sizes, t_grid = window.thresholds)
#' ## Print the result
#' print(dmrscan.results)
#' 
#' 
dmr_scan <- function(obs,k_grid,t_grid=NULL,...){

    xx          <- sapply(obs,length)
    L           <- sum(xx)
    alpha       <- 0.05

    k_grid  <- sort(k_grid)
    if(is.null(t_grid)){
        cat("Constructing t-grid\n")
        t_grid  <- estimate_t_grid(L,k_grid,...)
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
        sign_region <- lapply(index,function(x,val,which.k){
                                        rbind(val[,x],which.k[x])},
                                val=sliding.values,which.k=zhang.which.k)
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

     sliding.values.no.zero  <- do.call(cbind,apply(sliding.values,2,
                                            function(x){if(!sum(x==0)>0)x}))
     ll                      <- length(sliding.values.no.zero)

    ##

   phi                     <- acf(t.val[1,],plot = FALSE)$acf[2:3]
   p.val.empirical         <- apply(t.val,2,function(x,y,ll){(sum(abs(x[1]) <= y)+1)/ll},y = sliding.values.no.zero,ll = ll)
   sd                      <- sapply(sapply(t.val[2,],min,16),var_AR,phi=phi,s=1,p=2)
   #mean                    <- sapply(t.val[2,],mean_AR_1,phi=0.23,c=0.18) 
   mean                    <- median(apply(sliding.values.no.zero,1,mean))
   log.p.val.normal        <- -1*pnorm(abs(t.val[1,]), mean=mean, sd =sd,lower.tail=FALSE,log.p=TRUE)/log(10)

   if(is.null(rownames(obs[[1]]))){
       rr <- do.call(c,lapply(obs,names))
   }else{
       rr <- do.call(c,lapply(obs,rownames))
   }
   regions.coord   <- apply(regions,1,function(x,rr){c(rr[x[1]],rr[x[2]])},rr=rr)
   regions.coord   <- paste(regions.coord[1,],matrix(unlist(strsplit(regions.coord[2,],split="[.]")),ncol=2,byrow=TRUE)[,2],sep="-")
   regions.coord   <- gsub(regions.coord,pattern="[.]",replacement=":")

   regions         <- data.frame(regions,
                       "p.val.emp"    = p.val.empirical,
   #                    "logP.normal"  = log.p.val.normal,
                       "genome.coord" =regions.coord,stringsAsFactors = FALSE)
   regions         <- regions[order(regions$p.val.emp),]
   }else{
       regions          <- NA
    }


    return(regions)
}

