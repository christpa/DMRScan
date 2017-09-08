#' Object of type Region
#' 
#' @name Region-class 
#' @description
#' Class \code{Region} is a collection of test statistics for a set of CpGs 
#' within a short genomic range 
#' @rdname Region-class
#' @exportClass Region
setClass(Class = "Region",
    representation = representation(
        tValues     = "numeric",
        id          = "character",
        position    = "integer",
        chromosome  = "character",
        pVal        = "numeric",
        nCpG        = "integer")
)


#' Shorthand for initializing region
#' @name Region
#' @rdname Region-class_init
#' @importFrom methods new
#' @param tValues A vector of test statistics
#' @param position A vector of position for each test statistc
#' @param chromosome An character describing the chromosome (1-22, X,Y)
#' @param pVal The P value of a region, set to numeric() if not given.
#' @param id The names of each probe in the region 
#' @importFrom methods validObject
#' @return An object of type Region
#' @return An object of type Region
#' @export
#' @examples
#' #Number of probes is n = 10
#' nCpG <- 10
#' region <- Region(tValues    = rnorm(nCpG),
#'                  position   = 1:nCpG,
#'                  chromosome = "3",
#'                  id         = paste("CpG",1:nCpG,sep="_"),
#'                  pVal       = runif(1)) 
Region <- function(tValues, position, chromosome, pVal, id){ 
   if(missing(pVal)) pVal <- numeric()
   if(missing(id)) id   <- character(length(tValues))

   if(class(chromosome) %in% c("integer", "numeric"))

    chromosome  <- as.character(chromosome)
    nCpG       <- length(position)
    object     <- new("Region", 
    tValues   = tValues,
    id        = id,
    position  = position,
    chromosome= chromosome,
    pVal      = pVal,
    nCpG      = nCpG)
    return(object)
}
#' Class RegionList
#' Class \code{RegionList} is a collection of Regions  
#' @name RegionList-class 
#' @rdname RegionList-class
#' @exportClass RegionList
setClass(Class = "RegionList",
    representation = representation(
        regions     = "array", 
        nRegions    = "integer"
    )
)

#' Shorthand for initializing RegionList
#' @name RegionList
#' @rdname RegionList-class_init
#' @param regions The regions to be included
#' @param nRegions The number of regions to be placed
#' @importFrom methods new
#' @return An object of type RegionList
#' @export
#' @examples 
#' # An empty list of 3 regions
#' RegionList(3L)
#' 
RegionList <- function(nRegions,regions){ 
    if(missing(regions)) regions <- array(list(),nRegions)
    if(missing(nRegions)) nRegions <- length(regions)
         else if(!is.integer(nRegions)) nRegions <- as.integer(nRegions)
    object <- new("RegionList",regions = regions, nRegions = nRegions)
    return(object)
}
## Set generic functions for Rt() and St() when 
## One for class "list" and one for class "Region"
 
#' Method Fixed window size scan for one window size 
#' @name oneWindowSizeScanner
#' @rdname Rt-methods
#' @aliases Rt
#' @param region Object of type Region or RegionList
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
#' @param region Object of type Region or RegionList
#' @param windowThreshold Vector of window thresholds
#' @param windowSize Vector of window sizes to be tested on regions
#' @return A list of the windows that are significant  
#' @examples
#' ## Not run
#'
setGeneric("manyWindowSizeScanner", function(region,windowThreshold,windowSize) 
     standardGeneric("manyWindowSizeScanner"))


#' Method get pvalue 
#' @name pVal
#' @rdname pVal
#' @exportMethod pVal
#' @param region An object of type Region or RegionList
#' @param n The number of digits to be presented. Default is 10
#' @return The reported p-value for a region
#' @examples
#' #Number of probes is n = 10
#' nCpG <- 10
#' region <- Region(tValues    = rnorm(nCpG),
#'                  position   = 1:nCpG,
#'                  chromosome = "3",
#'                  pVal       = runif(1))
#' ## Pvalues for Region
#' pVal(region)
setGeneric("pVal", function(region,n = 12) standardGeneric("pVal"))

#' Method pos 
#' @name pos
#' @rdname pos
#' @exportMethod pos
#' @param region An opbject of type Region or RegionList
#' @return An integer vector of positions for each probe site
#' @examples
#' #Number of probes is n = 10
#' nCpG <- 10
#' region <- Region(tValues    = rnorm(nCpG),
#'                  position   = 1:nCpG,
#'                  chromosome = "3")
#' ## Genomic coordinates for Region
#' pos(region)
setGeneric("pos", function(region) standardGeneric("pos"))

#' Method get T statistic for a region
#' @name tVal
#' @rdname tVal
#' @exportMethod tVal
#' @param region An opbject of type Region or RegionList
#' @param ... Index 
#' @return A numeric vector of t-values for a Region or RegionList
#' @examples
#' #Number of probes is n = 10
#' nCpG <- 10
#' region <- Region(tValues    = rnorm(nCpG),
#'                  position   = 1:nCpG,
#'                  chromosome = "3")
#' ## T values for Region
#' tVal(region)
#'
setGeneric("tVal", function(region,...) standardGeneric("tVal"))


#' Method setRegion   
#' @name setRegion   
#' @rdname setRegion 
#' @exportMethod setRegion 
#' @param x A region
#' @param i an index
#' @param ... To be pased to Region()
#' @return An updated version of RegionList x, with a new Region at index i
#' @examples
#' ## A region list with 3 regions
#' regList <- RegionList(3L)
#' #Number of probes in first is n = 10
#' nCpG <- 10
#' region <- Region(tValues    = rnorm(nCpG),
#'                  position   = 1:nCpG,
#'                  chromosome = "3")
#' ## Set first region in regList to region
#' regList <- setRegion(regList,i = 1, region)
#'
setGeneric("setRegion", function(x,i,...) standardGeneric("setRegion"))


#' Method getRegions        
#' @name getRegions
#' @param x An object of type RegionList
#' @rdname getRegions
#' @aliases getRegions
#' @exportMethod getRegions                
#' @examples 
#' someEmptyRegions <- RegionList(3L) 
#' # To get back three empty regions
#' getRegions(someEmptyRegions)
#' @return An object of type Region
#'
setGeneric("getRegions", function(x) standardGeneric("getRegions"))

#' Method nCpG               
#' @name nCpG               
#' @rdname nCpG
#' @exportMethod nCpG      
#' @examples 
#' someEmptyRegions <- RegionList(3L)
#' # The number of CpGs in this regions is 0
#' nCpG(someEmptyRegions)
#'
#' @param x An opbject of type Region or RegionList
#' @return The number of CpGs in an object 
setGeneric("nCpG", function(x) standardGeneric("nCpG"))         


#' Get p-values for a region
#' @rdname pVal
setMethod("pVal", "Region",
          function(region,n=12){
              return(round(region@pVal,digits = n))
        }
)

#' Get p-values for a list of regions (RegionList)
#' @rdname pVal
setMethod("pVal", "RegionList",
          function(region,n = 12){
              return(do.call(c,lapply(getRegions(region),pVal,n)))
        }
)

#' Get the chromosomal coordinates for a Region
#' @rdname pos
setMethod("pos", "Region",
          function(region){
              return(region@position)
        }
)

#' Get the chromosomal coordinates for a 
#' list of regions in a RegionList object
#' @rdname pos
setMethod("pos", "RegionList",
          function(region){
              return(do.call(c,lapply(getRegions(region),pos)))
        }
)

#' Get test statistic for an object of type Region
#' @rdname tVal
setMethod("tVal", "Region",
          function(region, index = NULL){
              if(is.null(index))
                  return(region@tValues)
              else
                  return(region@tValues[index])
        }
)

#' Get test statistic for all regins within a RegionList class
#' @rdname tVal
#' @param index Index to extract
setMethod("tVal", "RegionList",
          function(region, index = NULL){
              if(is.null(index))
                  return(do.call(c,lapply(getRegions(region),tVal)))
              else
                  return(do.call(c,lapply(getRegions(region)[index],tVal)))
        }
)

#' getRegions for Region List
#' @name getRegions
#' @aliases getRegions,RegionList-method
#' @usage NULL
#' @return A region from a RegionList
#' @docType methods
#' @rdname getRegions
setMethod("getRegions","RegionList",
          function(x){
              return(x@regions)
        }
)

#' Get Object Region
#' @name [
#' @aliases [,RegionList,ANY,ANY,ANY-method
#' @usage NULL 
#' @param x An object of type RegionList 
#' @param i Index, which region to extract
#' @param j (Not used)
#' @param ... (not used)
#' @param drop If drop is used
#' @return A region from a RegionList of class "list"
#' @docType methods
#' @rdname index
setMethod("[", signature(x = "RegionList", i = "ANY", j = "ANY"),
    function (x, i, j, ..., drop){
       x@regions[i]
    }
)

#' Get Object Region
#' @name [[
#' @aliases [[,RegionList-method
#' @param x An object of type RegionList 
#' @param i Index, which region to extract
#' @param j (Not used)
#' @param ... (not used)
#' @param drop If drop is used
#' @return A region from a RegionList with class "Region"
#' @docType methods
#' @rdname index_list
setMethod("[[", signature(x = "RegionList", i = "ANY", j="ANY"),
    function (x, i, j, ..., drop){
       x@regions[[i]]
    }
)

#' Update a RegionList object
#' @rdname setRegion
#' @param region An object of type Region to be inseted in RegionList
setMethod("setRegion","RegionList",
          function(x,i,region){
              x@regions[[i]] <- region
              return(x)
            }
)

#' Calculate the length of a region in terms of CpGs
#' @rdname length
#' @return The number of CpGs in a Region
setMethod("length", "Region",
            function(x){
                return(x@nCpG)
            }
)

#' Get the number of regions in a RegionList
#' @rdname length
#' @param x A RegionList object
#' @return The number of CpGs in a RegionList
setMethod("length", "RegionList",
          function(x){
              return(x@nRegions)
        }
)

#' Get the genomic position of a Region
#' @rdname range
#' @param x An object of type Region
#' @return A character giving the genomic position
setMethod("range", "Region",
          function(x){
              pos <- pos(x)
            return(paste("Chr",x@chromosome,":",max(pos),"-",min(pos),sep=""))
          }
)

#' Get the number of CpGs i a region
#' @rdname nCpG
setMethod("nCpG","Region",
          function(x){
              return(x@nCpG)
            }
)

#' Get the number of CpGs in a RegionList
#' @rdname nCpG
setMethod("nCpG","RegionList",
          function(x){
              return(do.call(sum,lapply(getRegions(x),length)))
            }
)

#' Print a region
#' @rdname print
#' @param x Object of type Region
#' @param ... Has no function
#' @return An print object of a Region class
setMethod("print", "Region",
          function(x,...){
            out <- paste("Region with ", x@nCpG, " nCpGs on chromosome", 
            x@chromosome, ":", min(x@position), "-" ,max (x@position), 
            " With P value ", x@pVal, sep = "")
          return(print(out))
          }
)

#' Print a number of regions in a RegionList
#' @rdname print
#' @return A printed object of all regions in a RegionList
setMethod("print", "RegionList",
          function(x){
           cat("|Genomic Coordinate \t\t| nCpGs | pVal |\n") 
           sapply(getRegions(x),show)
           return(invisible(0))
          }
)

#' Show a region
#' @rdname show
#' @param object The region to be desplied, of type Region
#' @importFrom methods show
#' @return Cat a region to screen
setMethod("show", "Region",
          function(object){
          out <- paste("|Chr",object@chromosome, ":",min(object@position),"-",
            max(object@position)," \t|", object@nCpG, "\t|", object@pVal,"|\n",sep="") 
          return(invisible(cat(out)))
          }
)

#' Cat the head of a list of regions in a RegionList object
#' @rdname head
#' @param x An object to be printed of type RegionList
#' @param n The number of regions to be printed when the RegionList is 
#'  longer than n
#' @return The top regins in a RegionList
setMethod("head","RegionList",
          function(x,n = 10L){
            if(length(x) > n)
                x <- RegionList(regions = x[1:n],nRegions = n)
          return(invisible(print(x)))
          }
)

#' Sort a set of regions on p-value in a RegionList object
#' @rdname sort
#' @param x An object of type RegionList
#' @param decreasing Inherited from base
#' @return An updated RegionList, sorted on empirical p-values
setMethod("sort","RegionList",
          function(x, decreasing = FALSE){
            order       <- order(pVal(x))
            x@regions   <- x@regions[order]    
            return(x)


          }
)

#' Get the names of all probes within a region
#' @rdname names
#' @param x An object of type Region
#' @return The names of individual CpGs in a Region
setMethod("names","Region",
          function(x){
              return(x@id)
        }
)

#' Get the names of all probes in a study
#' @rdname names
#' @return A character vector of all CpG ids in a RegionList
setMethod("names","RegionList",
          function(x){
            return(do.call(c,lapply(getRegions(x),names)))
          }
)

