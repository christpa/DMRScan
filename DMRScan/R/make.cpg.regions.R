#' Bundles CpGs into a list of clusters 
#'
#' This function creates a list of regions/clusters of CpGs in close proximity, with no gap larger than max.gap, and 
#' all where all regions are longer than min.cpg. 
#' @param chr [integer] (Sorted) Vector of chromosome location for each CpG (use integer annotation 1-22, X = 23, Y = 24, XY = 25, MT = 26)
#' @param pos [integer] (Sorted) Vector giving base pair position for each CpG 
#' @param test_statistic [numeric] Vector of corresponding observed T-value for each CpG
#' @param max.gap [integer] Maximum allowed base pair gap within a cluster. Default is set to 500.
#' @param min.cpg [integer] Minimum number of CpGs allowed in each region to be considered. Default is set to at least 2 CpGs within each region
#' @keywords CpG Regions
#' @export
#' @examples
#' pos             <- cumsum(rpois(100,30))
#' chr             <- sort(sample(1:3,100,TRUE)) 
#' z_val           <- matrix(rnorm(100),dimnames=list(paste(chr,pos, sep = ".")))
#' 
#' min.cpg <- 3  # minimum number of cps for a region to be chnsidered
#' max.gap <- 30 # Maximum allowd gap between two observations within one region
#' regions <- make.cpg.regions(test_statistic=z_val,chr=chr,pos=pos,max.gap=max.gap,min.cpg=min.cpg)
#' ## Remaining observations 
#' sum(sapply(regions,nrow))
make.cpg.regions <- function(chr,pos,test_statistic,max.gap = 500, min.cpg = 2){
        if(any(is.na(c(chr,pos))))
         ## if any NA or Chr is not integer
         ## stop
        
	if(is.null(dim(test_statistic))){nn <- names(test_statistic);test_statistic <- as.matrix(test_statistic);rownames(test_statistic)<-nn}
        rowNames    <- rownames(test_statistic)

        idx          <- split(seq_along(pos), chr)
        region.idx   <- rep(NA,length(pos))
        end          <- 0
        
        for(i in seq_along(idx)){
             oo      <- order(pos[idx[[i]]])
             idx.reg <- idx[[i]][oo]
             xx      <- pos[idx.reg]
            
             diff    <- as.numeric(diff(xx) > max.gap)
             cluster <- cumsum(c(1,diff))
             region.idx[idx.reg] <- cluster + end
             end     <- max(cluster) + end
         }



        n.reg               <- max(region.idx)
        region.list         <-array(list(),n.reg)
        rr.true             <- !is.null(rownames(test_statistic))

        for(i in seq_along(region.list)){
            reg.id  <- which(region.idx==i)
            xx     <- test_statistic[reg.id,]
            region.list[[i]]  <- as.matrix(xx)
                if(rr.true)
                    rownames(region.list[[i]]) <- rowNames[reg.id]

            cat(i, "/", n.reg,"\r")
        };cat("\n") 
        ll               <- sapply(region.list,length)
        region.list.long <- region.list[ll >= min.cpg]

        return(region.list.long)
}
