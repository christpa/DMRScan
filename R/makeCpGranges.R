#' @title Cluster
#' @name makeCpGregions
#' @description Cluster CpGs together in regions based on proximity
#' @param chr Vector of chromosome location for each CpG 
#' @param pos Vector giving base pair position for each CpG If unsorted, 
#' use order(chr,pos) to sort the genomic positions within each chromosome.
#' @param observations Vector of corresponding observed T-value for each CpG, 
#' must be ordered in the same way as chr and pos
#' @param maxGap Maximum allowed base pair gap within a cluster. 
#' Default is set to 500.
#' @param minCpG Minimum number of CpGs allowed in each region to be 
#' considered. Default is set to at least 2 CpGs within each region.
#' @return The suplied observations ordered into into a GRangesList object. 
#' To be parsed further into \code{\link{dmrscan}} 
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom methods as
#' @keywords CpG Regions
#' @export
#' @examples
#' data(DMRScan.methylationData) ## Load methylation data from chromosome 22
#' data(DMRScan.phenotypes) ## Load phenotype (end-point for methylation data)
#' 
#' ## Test for an association between phenotype and Methylation
#' testStatistics <- apply(DMRScan.methylationData,1,function(x,y)
#'  summary(glm(y ~ x, family = binomial(link = "logit")))$coefficients[2,3],
#'  y = DMRScan.phenotypes)
#' 
#' ## Set chromosomal position to each test-statistic
#' pos<- data.frame(matrix(as.integer(unlist(strsplit(names(testStatistics),
#' split="chr|[.]"))), ncol = 3, byrow = TRUE))[,-1] 
#' 
#' ## Set clustering features 
#' minCpG <- 3  ## Minimum number of CpGs in a tested cluster 
#' ## Maxium distance (in base-pairs) within a cluster before it is 
#' ## broken up into two seperate cluster
#' maxGap <- 750  
#' regions <- makeCpGregions(observations = testStatistics, chr = pos[,1], 
#'                             pos = pos[,2], maxGap = maxGap, minCpG = minCpG)
makeCpGregions <- function(observations, chr, pos, maxGap = 500, minCpG = 2){

    if(is.null(dim(observations))){
        rowNames <- names(observations)
        observations <- as.matrix(observations)
        rownames(observations) <- rowNames
    }else{
        rowNames <- rownames(observations)
    }

    idx          <- split(seq_along(pos), chr)
    region.idx   <- rep(NA,length(pos))
    end          <- 0
    
    for(i in seq_along(idx)){
         oo      <- order(pos[idx[[i]]])
         idx.reg <- idx[[i]][oo]
         xx      <- pos[idx.reg]
        
         diff    <- as.numeric(diff(xx) > maxGap)
         cluster <- cumsum(c(1,diff))
         region.idx[idx.reg] <- cluster + end
         end     <- max(cluster) + end
    }

    

## Substetute reginos.index with cluster
    if(minCpG >= 2){
        nRegions    <- sum(table(region.idx) >= minCpG)
    }else{
        nRegions    <- as.integer(max(region.idx))
    }
        regionList  <- array(list(),nRegions)

    j                  <- 1
    for(i in seq_len(max(region.idx))){
        reg.id  <- which(region.idx==i)
        xx     <- observations[reg.id,]
        if(length(reg.id) >= minCpG){
            regionList[[j]] <- GRanges(seqnames = chr[reg.id],
                                  	   ranges = IRanges(start=pos[reg.id], end = pos[reg.id]),
										no.cpgs= length(reg.id),
										tVal = xx,
										id   = rowNames[reg.id])
            #cat(j, "/", nRegions,"\r")
                j   <- j + 1 
        }
    };cat("\n") 
	if(i > (j-1))
		regionList <- regionList[1:(j-1)]

	regionList <- as(unlist(regionList), "GRangesList")
    return(regionList)
}


#' @title Cluster
#' @name makeCpGgenes
#' @description
#' Cluster CpGs together in genes based on annotation
#' @param chr Vector of chromosome location for each CpG 
#' @param pos Vector giving base pair position for each CpG If unsorted, 
#' use order(chr,pos) to sort the genomic positions within each chromosome.
#' @param observations Vector of corresponding observed T-value for each CpG, 
#' must be ordered in the same way as chr and pos
#' @param gene A vector asigning each probe to a gene.
#' @param minCpG Minimum number of CpGs allowed in each region to be 
#' considered. Default is set to at least 2 CpGs within each region.
#' @return The suplied observations ordered into into a list, 
#' with one entry for each CpG region.
#' @keywords CpG Regions
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom methods as
#' @keywords CpG Regions
#' @export
#' @examples
#' data(DMRScan.methylationData) ## Load methylation data from chromosome 22
#' data(DMRScan.phenotypes) ## Load phenotype (end-point for methylation data)
#' 
#' ## Test for an association between phenotype and Methylation
#' testStatistics <- apply(DMRScan.methylationData,1,function(x,y)
#'  summary(glm(y ~ x, family = binomial(link = "logit")))$coefficients[2,3],
#'   y = DMRScan.phenotypes)
#' 
#' ## Set chromosomal position to each test-statistic
#' pos <- data.frame(matrix(as.integer(unlist(strsplit(names(testStatistics),
#'  split="chr|[.]"))), ncol = 3, byrow = TRUE))[,-1]
#' 
#' ## Set clustering features 
#' minCpG   <- 3  ## Minimum number of CpGs in a tested cluster
#' gene     <- sample(paste("Gene",1:100,sep=""), 
#'                            length(testStatistics),replace=TRUE)
#' regions  <- makeCpGgenes(observations = testStatistics, 
#'                          chr = pos[,1], pos = pos[,2], 
#'                          gene = gene, minCpG = minCpG)
#'
makeCpGgenes <- function(observations, chr, pos, gene, minCpG = 2){

        if(anyNA(c(chr,pos))){
         ## if any NA or Chr is not integer
         ## stop
        }
        if(is.factor(gene))
            gene    <- as.character(gene)
        
        if(is.null(dim(observations))){
            rowNames <- names(observations)
            observations <- as.matrix(observations)
            rownames(observations) <- rowNames
        }else{
            rowNames <- rownames(observations)
        }
            if(is.numeric(gene) | is.factor(gene))
                gene <- as.character(gene)

## Clean gene names
        gene[grep(" ", gene)] <- NA
        gene[gene == ""] <- NA

        
      geneNames  <- unique(na.omit(gene))
      nGene      <- sum(table(gene) >= minCpG)
      
      regionList  <- array(list(),nGene)

      i          <- 1
      for(geneIter in geneNames){
          idx           <- gene %in% geneIter
          xx            <- observations[idx,]
        
        if(length(idx) >= minCpG){
            regionList[[i]] <- GRanges(seqnames = chr[idx],
                                  	   ranges = IRanges(start=pos[idx], end = pos[idx]),
										no.cpgs = length(idx),
										tVal = xx,
										id   = rowNames[idx],
										gene = geneIter)
            #cat(j, "/", nRegions,"\r")
                i   <- i + 1 
        }
    };cat("\n") 
	if(length(geneNames) > (i-1))
		regionList <- regionList[1:(i-1)]

	regionList <- as(unlist(regionList), "GRangesList")
    return(regionList)

}
