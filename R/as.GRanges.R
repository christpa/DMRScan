




setGeneric("as.GRanges", function(x)
                    standardGeneric("as.GRanges"))


setMethod("as.GRanges", "Region",
      function(x){
			start <- min(x@position)
			stop  <- max(x@position)
			names <- paste("chr", x@chromosome,":",start, "-",stop, sep="")
			return(GRanges(
				seqnames = paste("chr",x@chromosome,sep=""),
				ranges   = IRanges(start, end = stop, names = names),
				nCpGs	 = length(x),
				pVal	 = min(pVal(x),1)
			))                 
                                    }
                                    )

#' Get p-values for a list of regions (RegionList)
#' @rdname as.GRanges
setMethod("as.GRanges", "RegionList",
          function(x){
                if(length(x) == 1)
                    return(as.GRanges(x[[1]]))
                else{
                  x     <- getRegions(x)
                  chr   <- sapply(x, function(x)x@chromosome)
                  start <- sapply(x,function(x)min(x@position))
				  stop  <- sapply(x,function(x)max(x@position))
				  names <- paste("chr", chr,":",start, "-",stop, sep="")
				  
                    return(GRanges(        
                    seqnames = chr,
                    ranges   = IRanges(start, end = stop, names = names),
                    nCpGs    = sapply(x,length),
                    pVal     = sapply(x,function(x)min(pVal(x),1))
                    ))
                }
                                    }
                                    )
