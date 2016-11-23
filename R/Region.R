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

#' Constructor for class Region  
#' 
#' @name initialize-Region
#' @rdname initialize_Region 
#' @param .Object Empty
#' @param tValues A vector of test statistics
#' @param id A vector of id's for each test statistic
#' @param position A vector of position for each test statistc
#' @param chromosome An character describing the chromosome (1-22, X,Y)
#' @param pVal A numeric giving the estimated p-value for each significant DMR
#' @importFrom methods validObject
#' @return An object of type Region
setMethod("initialize", "Region", 
          function(.Object, tValues = "numeric", id = "character", 
                   position = "integer",chromosome = "character", 
                   pVal = "numeric"){
        
            if(!missing(pVal)) 
                .Object@pVal <- pVal
            
            if(class(chromosome) %in% c("integer", "numeric"))
                chromosome  <- as.character(chromosome)
            
            .Object@tValues     <- tValues
            .Object@position    <- position
            .Object@chromosome  <- chromosome[1]
            .Object@nCpG        <- length(tValues)
            
            if(!missing(id)){
                .Object@id      <- id
            }else{
                .Object@id      <- character(length(tValues))
            }

            validObject(.Object)
            return(.Object)
        }
)

#' Shorthand for initializing region
#' @name Region-init
#' @rdname Region-class_init
#' @importFrom methods new
#' @param ... Parameters to pased to new("Regions",...)
#' @return An object of type Region
#' @export
#' @examples
#' #Number of probes is n = 10
#' nCpG <- 10
#' region <- Region(tValues    = rnorm(nCpG),
#'                  position   = 1:nCpG,
#'                  chromosome = "3")
Region <- function(...) new("Region", ...)

#' Class RegionList                                                        
#'                                                                     
#' Class \code{RegionList} is a collection of Regions  
#'                                                                     
#' @name RegionList-class                                                  
#' @rdname RegionList-class                                               
#' @exportClass RegionList                                                 
setClass(Class = "RegionList",                                             
    representation = representation(                                   
        regions     = "array",                                       
        nRegions    = "integer"
    )
)

#' Constructor RegionList
#' @name initialize-RegionList 
#' @rdname initialize_RegionList
#' @param .Object Empty
#' @param nRegions The number of regions to be included in the RegionList object
#' @param regions An optional list of regions to be assigned to a RegionList
#' @importFrom methods validObject
#' @return An object of tpe RegionList
setMethod("initialize", "RegionList",    
          function(.Object, nRegions = "integer",regions = "array"){
            
            if(missing(regions)){
                .Object@regions <- array(list(),nRegions)
            }else{
                .Object@regions <- regions
            }
            if(missing(nRegions)){ 
                .Object@nRegions <- length(regions)
            }else{
                if(!is.integer(nRegions))nRegions <- as.integer(nRegions)
                .Object@nRegions <- nRegions
            }
            return(.Object)
        }
)

#' Shorthand for initializing RegionList
#' @name RegionList                                        
#' @rdname RegionList-class_init
#' @param ... Parameters to be past to new("RegionList",...)
#' @importFrom methods new
#' @return An object of type RegionList
#' @export
#' @examples 
#' # An empty list of 3 regions
#' RegionList(3L)
#' 
RegionList <- function(...){ 
                new("RegionList",...)         
}
 
## Set generic functions for Rt() and St() when 
## One for class "list" and one for class "Region"
 
#' Method Rt 
#' @name Rt
#' @rdname Rt-methods
#' @param region Object of type Region or RegionList
#' @param windowThreshold Vector of window thresholds
#' @param windowSize Vector of window sizes to be tested on regions
#' @return A list of which windows that are significant
#' @examples
#' ## Not run
#'
setGeneric("Rt", function(region,windowThreshold,windowSize) 
                    standardGeneric("Rt"))


#' Method St 
#' @name St
#' @rdname St-methods
#' @param region Object of type Region or RegionList
#' @param windowThreshold Vector of window thresholds
#' @param windowSize Vector of window sizes to be tested on regions
#' @return A list of which windows that are significant  
#' @examples
#' ## Not run
#'
setGeneric("St", function(region,windowThreshold,windowSize) 
                                standardGeneric("St"))


#' Method getP 
#' @name getP
#' @rdname getP
#' @exportMethod getP
#' @param region An object of type Region or RegionList
#' @param n The number of digits to be presented. Default is 10
#' @return A numeric vector of p-values
#' @examples
#' #Number of probes is n = 10
#' nCpG <- 10
#' region <- Region(tValues    = rnorm(nCpG),
#'                  position   = 1:nCpG,
#'                  chromosome = "3",
#'                  pVal       = runif(1))
#' ## Pvalues for Region
#' getP(region)
setGeneric("getP", function(region,n=10) standardGeneric("getP"))

#' Method getPos 
#' @name getPos
#' @rdname getPos
#' @exportMethod getPos
#' @param region An opbject of type Region or RegionList
#' @return An integer vector of positions for each probe site
#' @examples
#' #Number of probes is n = 10
#' nCpG <- 10
#' region <- Region(tValues    = rnorm(nCpG),
#'                  position   = 1:nCpG,
#'                  chromosome = "3")
#' ## Genomic coordinates for Region
#' getPos(region)
setGeneric("getPos", function(region) standardGeneric("getPos"))

#' Method getT 
#' @name getT
#' @rdname getT
#' @exportMethod getT
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
#' getT(region)
#'
setGeneric("getT", function(region,...) standardGeneric("getT"))


#' Method setRegion   
#' @name setRegion   
#' @rdname setRegion 
#' @exportMethod setRegion 
#' @param x A region
#' @param i an index
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
#' @rdname getRegions
#' @param x A RegionList object
#' @exportMethod getRegions                
#' @examples 
#' someEmptyRegions <- RegionList(3L) 
#' # To get back three empty regions
#' getRegions(someEmptyRegions)
#' 
#' @return An object of type Region
#' @examples
#' ## Not run
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
#' @rdname getP
setMethod("getP", "Region",
          function(region,n=10){
              return(print(region@pVal,digits = n))
        }
)

#' Get p-values for a list of regions (RegionList)
#' @rdname getP
setMethod("getP", "RegionList",
          function(region,n = 10){
              return(do.call(c,lapply(getRegions(region),getP,n)))
        }
)

#' Get the chromosomal coordinates for a Region
#' @rdname getPos
setMethod("getPos", "Region",
          function(region){
              return(region@position)
        }
)

#' Get the chromosomal coordinates for a 
#' list of regions in a RegionList object
#' @rdname getPos
setMethod("getPos", "RegionList",
          function(region){
              return(do.call(c,lapply(getRegions(region),getPos)))
        }
)

#' Get test statistic for an object of type Region
#' @rdname getT
setMethod("getT", "Region",
          function(region, index = NULL){
              if(is.null(index))
                  return(region@tValues)
              else
                  return(region@tValues[index])
        }
)

#' Get test statistic for all regins within a RegionList class
#' @rdname getT
#' @param index Index to extract
setMethod("getT", "RegionList",
          function(region, index = NULL){
              if(is.null(index))
                  return(do.call(c,lapply(getRegions(region),getT)))
              else
                  return(do.call(c,lapply(getRegions(region)[index],getT)))
        }
)

#' Get a individual region within an object of class RegionList
#' @name [
#' @rdname index
#' @param x An object of type RegionList
#' @param i Index, which region to extract
#' @return A region from a RegionList of class "list"
setMethod("[", signature(x = "RegionList", i = "ANY"),
            function(x,i){
                x@regions[i]
            }
)

#' Get a region from a RegionList object
#' @name [[
#' @rdname index_list 
#' @param x An object of type RegionList      
#' @param i Index, which region to extract
#' @return A region from a RegionList with class "Region"
setMethod("[[", signature(x = "RegionList", i = "ANY"),
            function(x,i){
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

#' Get the regions within an RegionList object
#' @rdname getRegions
setMethod("getRegions","RegionList",
          function(x){
             return(x@regions)
            }
)

#' Set a Region of type missing
#' @rdname setRegion
#' @param ... Arguments to be pased to Region(...)
setMethod("setRegion","missing",
          function(x,i,...){
              x@regions[[i]] <- Region(...)
              return(x)
            }
)

#' Calculate the length of a region in terms of CpGs
#' @rdname length
#' @param x An object of type Region
#' @return The number of CpGs in a Region
setMethod("length", "Region",
            function(x){
                return(x@nCpG)
            }
)

#' Get the number of regions in a RegionList
#' @rdname length
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
              pos <- getPos(x)
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
            out <- paste("Region with ", x@nCpG, " nCpGs on chromosome", x@chromosome, 
                         ":", min(x@position), "-" ,max (x@position), " With P value ", x@pVal, sep = "")
       #     print(out) 

          return(print(out))
          }
)

#' Print a number of regions in a RegionList
#' @rdname print
#' @return A printed object of all regions in a RegionList
setMethod("print", "RegionList",
          function(x){
           cat("|Genomic Coordinate \t\t| nCpGs | pVal |\n") 
           lapply(getRegions(x),show)
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
      #    cat(out)
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
            order       <- order(getP(x))
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

