########################
#name : "Comparing the DMRScan package to other DMR calling algorithms" 
#author: "Christian M Page"                                                      
#abstract: Comparison of the three methods; DMRScan, DMRCate and bumphunter. 
#########################

## Load Packages and data
library(DMRScan)                                                                
data(DMRScan.methylationData) ## Load methylation data from chromosome 22, with 52018 CpGs measured
data(DMRScan.phenotypes) ## Load phenotype (end-point for methylation data)     
library(bumphunter)
library(DMRcate)


getPos <- function(x,id)if{(length(x) > 1){lapply(sapply(apply(as.data.frame(ranges(x)),1,function(x)seq(x[1],x[2])),function(x)paste("chr22.",x,sep="")),intersect,id)}else{
                                                list(intersect(paste("chr22.",ranges(x)[[1]],sep=""),id))}}


## Simulate true positive clusters of DMRs   

rr              <- rownames(DMRScan.methylationData)
rr.int          <- as.integer(matrix(unlist(strsplit(rr,"[.]")),ncol=2,byrow=TRUE)[,2])

obs             <- rep(-9,length(rr));names(obs)<-rr
cgIsland        <- makeCpGregions(obs, pos = rr.int, chr = rep("22",length(rr.int)), minCpG = 15)

## Pick regions to add effect to 
region.index    <- order(sapply(cgIsland,length),decreasing=TRUE)[1:100]
region.length   <- c(100,80,60,60,50,50,45,45,rep(c(40,35,30,25,20,15,10,5),c(4,4,6,8,20,20,20,10)))
region.id.list  <- sapply(1:100,function(x,r,ll)names(r[[x]])[31:(ll[x] + 30)],r = cgIsland[region.index],ll = region.length)
region.id       <- unlist(region.id.list) 
region.index    <- which(rr %in% region.id)

beta            <- DMRScan.methylationData
names(rr)       <- rr

window.sizes    <- 3:9

res             <- matrix(0,21,2,dimnames=list(seq(0,0.2,0.01),c("TP","FP")))
regions.observed <- list(DMRScanSigmund = res, DMRScanImpSamp = res, DMRScanMCMC = res, bumphunter = res, DMRCate = res);rm(res) 
cpg.observed     <- regions.observed

test.statistics <- apply(beta, 1, function(x,y)summary(glm(y ~ x, family = binomial()))$coef[2,3],y = DMRScan.phenotypes)
for(i in 1:21){
## Do glm/lm -> DMRScan
print(i)

test.statistics[region.id] <- apply(beta[region.id,], 1, function(x,y)summary(glm(y ~ x, family=binomial()))$coef[2,3],y = DMRScan.phenotypes)
test.statistics.regions     <- makeCpGregions(test.statistics, pos = rr.int, chr = rep("22",length(rr)), maxGap = 500,minCpG = 5)

DMRScan.res                 <- lapply(window.thresholds,function(x)dmrscan(observations = test.statistics.regions, windowSize = window.sizes, windowThreshold = x))


#  Do bumphunter
bumphunter.res  <- bumphunterEngine(mat = beta, design = cbind(1,DMRScan.phenotypes), chr = rep(22,nrow(beta)), pos = rr.int, pickCutoff = TRUE, B = 1000,  maxGap = 500, verbose=FALSE)$table
bumphunter.res  <- bumphunter.res[bumphunter.res$fwer < 0.05,]
## Do DMRCate (for sequencing)
wgbs.diff   <- rowSums(beta[,DMRScan.phenotypes == 1],na.rm=TRUE) - rowSums(beta[,DMRScan.phenotypes == 0], na.rm=TRUE)
wgbsannot   <- cpg.annotate("sequencing", object = data.frame("stat" = test.statistics, "chr" = "22",  "pos" = rr.int, "diff" = wgbs.diff, "fdr" = p.adjust(2*pnorm(abs(test.statistics), lower.tail = FALSE), method = "fdr")))
DMRCate.res <- dmrcate(wgbsannot, lambda = 1000, C = 50, pcutoff = 0.05, mc.cores = 1)$results



## True/false observations in list:
# DMRScan
if(length(DMRScan.res$siegmund) > 0){
       DMRScan.pos.list    <- getPos(DMRScan.res$siegmund,rr)
       DMRScan.pos     <- unlist(DMRScan.pos.list) 
    ## True positivies
    regions.observed$DMRScanSigmund[i,1] <- sum(sapply(region.id.list,function(x,pos)any(x%in%pos),pos = rr[DMRScan.pos]))
    cpg.observed$DMRScanSigmund[i,1]     <- 100*sum(sapply(region.id.list,function(x,pos)sum(x%in%pos),pos = rr[DMRScan.pos]))/length(region.index)
    ## False positive 
    regions.observed$DMRScanSigmund[i,2] <- sum(sapply(DMRScan.pos.list,function(x,pos,rr)!any(rr[x] %in% pos), pos = region.id, rr = rr))
    cpg.observed$DMRScanSigmund[i,2]     <- sum(sapply(DMRScan.pos.list,function(x,pos,rr)sum(! rr[x] %in% pos), pos = region.id, rr = rr))
}


if(length(DMRScan.res$importancSampling) > 0){

    DMRScan.pos.list    <- getPos(DMRScan.res$importancSampling,rr)
    DMRScan.pos         <- unlist(DMRScan.pos.list)
             
    ## True positivies
    regions.observed$DMRScanImpSamp[i,1] <- sum(sapply(region.id.list,function(x,pos)any(x%in%pos),pos = rr[DMRScan.pos]))
    cpg.observed$DMRScanImpSamp[i,1]     <- 100*sum(sapply(region.id.list,function(x,pos)sum(x%in%pos),pos = rr[DMRScan.pos]))/length(region.index)
    ## False positive 
    regions.observed$DMRScanImpSamp[i,2] <- sum(sapply(DMRScan.pos.list,function(x,pos,rr)!any(rr[x] %in% pos), pos = region.id, rr = rr))
    cpg.observed$DMRScanImpSamp[i,2]     <- sum(sapply(DMRScan.pos.list,function(x,pos,rr)sum(! rr[x] %in% pos), pos = region.id, rr = rr))
}

if(length(DMRScan.res$mcmc) > 0){
        
        DMRScan.pos.list    <- getPos(DMRScan.res$mcmc,rr)
        DMRScan.pos         <- unlist(DMRScan.pos.list)
    
    ## True positivies
    regions.observed$DMRScanMCMC[i,1] <- sum(sapply(region.id.list,function(x,pos)any(x%in%pos),pos = rr[DMRScan.pos]))
    cpg.observed$DMRScanMCMC[i,1]     <- 100*sum(sapply(region.id.list,function(x,pos)sum(x%in%pos),pos = rr[DMRScan.pos]))/length(region.index)
    ## False positive 
    regions.observed$DMRScanMCMC[i,2] <- sum(sapply(DMRScan.pos.list,function(x,pos,rr)!any(rr[x] %in% pos), pos = region.id, rr = rr))
    cpg.observed$DMRScanMCMC[i,2]     <- sum(sapply(DMRScan.pos.list,function(x,pos,rr)sum(! rr[x] %in% pos), pos = region.id, rr = rr))
}

if(nrow(bumphunter.res) > 0){
    cat("Bumpts = ", nrow(bumphunter.res),"\n")
    if(nrow(bumphunter.res) == 1){
         bumphunter.index.list <- list(seq(bumphunter.res$indexStart, bumphunter.res$indexEnd)) 
    }else{
        bumphunter.index.list <- apply(cbind(bumphunter.res$indexStart, bumphunter.res$indexEnd),1,function(x)seq(from = x[1], to = x[2]))
    }
    bumphunter.index      <- sort(unlist(bumphunter.index.list))

    regions.observed$bumphunter[i,1] <- sum(sapply(region.id.list,function(x,pos)any(x%in%pos),pos = rr[bumphunter.index]))
    cpg.observed$bumphunter[i,1]     <- 100*sum(sapply(region.id.list,function(x,pos)sum(x%in%pos),pos = rr[bumphunter.index]))/length(region.index)
    ## False positive
    regions.observed$bumphunter[i,2] <- sum(sapply(bumphunter.index.list,function(x,pos,rr)!any(rr[x] %in% pos)>0, pos = region.id, rr = rr))
    cpg.observed$bumphunter[i,2]     <- sum(sapply(bumphunter.index.list,function(x,pos,rr)sum(! rr[x] %in% pos), pos = region.id, rr = rr))


}

# DMRCate 
if(nrow(DMRCate.res) > 0){
    
    if(nrow(DMRCate.res) == 1){
        DMRCate.pos.list <- list(t(apply(as.matrix(which(rr %in%  matrix(paste("chr22.",unlist(strsplit(DMRCate.res$coord,c(":|-"))), sep=""),ncol=3, byrow=TRUE)[,-1])), 2,function(x)seq(x[1],x[2],1))))    
    }else{
        DMRCate.pos.list <- apply(apply(matrix(paste("chr22.",unlist(strsplit(DMRCate.res$coord,c(":|-"))), sep=""),ncol=3, byrow=TRUE)[,-1],1,function(x,rr)which(rr %in% x),rr = rr), 2,function(x)seq(x[1],x[2],1))
    }
    DMRCate.pos <- unlist(DMRCate.pos.list)

    regions.observed$DMRCate[i,1] <- sum(sapply(region.id.list,function(x,pos)any(x%in%pos),pos = rr[DMRCate.pos]))
    cpg.observed$DMRCate[i,1]             <- 100*sum(sapply(region.id.list,function(x,pos)sum(x%in%pos),pos = rr[DMRCate.pos]))/length(region.index)
    ## False positive 
    regions.observed$DMRCate[i,2] <- sum(sapply(DMRCate.pos.list,function(x,pos,rr) !any(rr[x] %in% pos) , pos = region.id, rr = rr))
    cpg.observed$DMRCate[i,2]     <- sum(sapply(DMRCate.pos.list,function(x,pos,rr)sum(! rr[x] %in% pos), pos = region.id, rr = rr))
}


 ## Updating beta last, make sure that the "first update" is zero!
up <- rowMeans(beta[region.id,],na.rm=TRUE)<0.5
beta[region.id[up], DMRScan.phenotypes == 1] <- pmin(beta[region.id[up], DMRScan.phenotypes == 1] + 0.01,1)
beta[region.id[!up], DMRScan.phenotypes == 0] <- pmax(beta[region.id[!up], DMRScan.phenotypes == 0] - 0.01,0)
}

