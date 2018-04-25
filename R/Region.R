# Set generic functions for Rt() and St() when 
## One for class "list" and one for class "Region"
 
#' Method Fixed window size scan for one window size 
#' @name oneWindowSizeScanner
#' @rdname Rt-methods
#' @aliases Rt
#' @param region Object of type GRanges 
#' @param windowThreshold Vector of window thresholds
#' @param windowSize Vector of window sizes to be tested on regions
#' @return A list of the windows that are significant
#' @examples
#' ## Not run
#'
setGeneric("oneWindowSizeScanner", function(region,windowThreshold,windowSize) 
     standardGeneric("oneWindowSizeScanner"))

#' Method Fixed window size scan for a sequence of window sizes
#' @name manyWindowSizeScanner
#' @rdname St-methods
#' @aliases St
#' @param region Object of type GRanges 
#' @param windowThreshold Vector of window thresholds
#' @param windowSize Vector of window sizes to be tested on regions
#' @return A list of the windows that are significant  
#' @examples
#' ## Not run
#'
setGeneric("manyWindowSizeScanner", function(region,windowThreshold,windowSize) 
     standardGeneric("manyWindowSizeScanner"))


setGeneric("chr",function(x)
	standardGeneric("chr"))

setGeneric("pos",function(x)
	standardGeneric("pos"))

setGeneric("tVal",function(x)
	standardGeneric("tVal"))

setGeneric("id",function(x)
	standardGeneric("id"))

#' @importFrom GenomeInfoDb seqnames
setMethod("chr", "GRangesList", function(x){
	return(rep(do.call(c,lapply(x, function(x)as.character(seqnames(x)@values))), 
				do.call(c,lapply(x, function(x)seqnames(x)@lengths))))
})

#' @importFrom GenomicRanges ranges 
setMethod("pos", "GRangesList", function(x){
	return(do.call(c,sapply(x, function(x)ranges(x)@start)))
})


setMethod("tVal", "GRanges", function(x){
	return(mcols(x)$tVal)
})

setMethod("tVal", "GRangesList", function(x){
	return(do.call(c,sapply(sapply(x, mcols),function(x)x$tVal)))
})


setMethod("id", "GRangesList", function(x){
	return(do.call(c,sapply(sapply(x, mcols),function(x)x$tVal)))
})

