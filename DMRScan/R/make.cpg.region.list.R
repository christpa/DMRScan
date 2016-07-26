make.cpg.region.list <- function(chr,pos,test_statistic,max.gap,min.cpg) {
## Based on ClusterMaker in "Bumphunter"

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

