#' Bundles CpGs into a list of clusters 
#'
#' This function creates a list of regions/clusters of CpGs in close proximity, with no gap larger than max.gap, and 
#' all where all regions are longer than min.cpg. 
#' @param chr (Sorted) Vector of chromosome location for each CpG (use integer annotation 1-22, X = 23, Y = 24, XY = 25, MT = 26)
#' @param pos (Sorted) Vector giving base pair position for each CpG 
#'  If unsorted, use order(chr,pos) to sort the genomic positions within each chromosome.
#' @param test_statistic Vector of corresponding observed T-value for each CpG, must be ordered in the same way as chr and pos
#' @param max.gap Maximum allowed base pair gap within a cluster. Default is set to 500.
#' @param min.cpg Minimum number of CpGs allowed in each region to be considered. Default is set to at least 2 CpGs within each region.
#' @return The suplied test_statistics ordered into into a list, with one entry for each CpG region.
#' @keywords CpG Regions
#' @export
#' @examples
#' data(DMRScan.methylationData) ## Load methylation data from chromosome 22
#' data(DMRScan.phenotypes) ## Load phenotype (end-point for methylation data)
#' 
#' ## Test for an association between phenotype and Methylation
#' test.statistics <- apply(DMRScan.methylationData,1,function(x,y)summary(glm(y ~ x, family = binomial(link = "logit")))$coefficients[2,3], y = DMRScan.phenotypes)
#' pos<- data.frame(matrix(as.integer(unlist(strsplit(names(test.statistics),split="chr|[.]"))), ncol = 3, byrow = TRUE))[,-1] ## Set chromosomal position to each test-statistic
#' 
#' ## Set clustering features 
#' min.cpg <- 3  ## Minimum number of CpGs in a tested cluster
#' max.gap <- 750  ## Maxium distance (in base-pairs) within a cluster before it is broken up into two seperate cluster 
#' regions <- make.cpg.regions(test_statistic = test.statistics, chr = pos[,1], pos = pos[,2], max.gap = max.gap, min.cpg = min.cpg)
make.cpg.regions <- function(chr,pos,test_statistic,max.gap = 500, min.cpg = 2){
        if(any(is.na(c(chr,pos)))){
         ## if any NA or Chr is not integer
         ## stop
        }

	if(is.null(dim(test_statistic))){
        nn <- names(test_statistic);test_statistic <- as.matrix(test_statistic)
        rownames(test_statistic)<-nn
    }
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
