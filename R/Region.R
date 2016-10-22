#' Class Region
#' 
#' Class \code{Region} is a collection of test statistics for a set of CpGs within a short genomic range 
#'  
#' @name Region-class 
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
#' @name Region
#' @rdname Region-class
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
#' @name Region
#' @rdname Region-class
#' @export 
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
                                              
                                                                                
#' Constructor for class RegionList
#'                                                                              
#' @name RegionList                                             
#' @rdname RegionList-class                                                         
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
#' @rdname RegionList-class                                                         
#' @export                                                                      
RegionList <- function(...){ 
                new("RegionList",...)         
}
 
## Set generic functions for Rt() and St() when 
## One for class "list" and one for class "Region"
 
#' Method Rt 
#' @name Rt
#' @rdname Rt-methods
#' @exportMethod Rt
setGeneric("Rt", function(region,windowThreshold,windowSize) standardGeneric("Rt"))


#' Method St 
#' @name St
#' @rdname St-methods
#' @exportMethod St
setGeneric("St", function(region,windowThreshold,windowSize) standardGeneric("St"))


#' Method getP 
#' @name getP
#' @rdname getP-methods
#' @exportMethod getP
setGeneric("getP", function(region) standardGeneric("getP"))

#' Method getPos 
#' @name getPos
#' @rdname getPos-methods
#' @exportMethod getPos
setGeneric("getPos", function(region) standardGeneric("getPos"))



#' Method getT 
#' @name getT
#' @rdname getT-methods
#' @exportMethod getT
setGeneric("getT", function(region,...) standardGeneric("getT"))


#' Method setRegion                                                                  
#' @name setRegion                                                             
#' @rdname setRegion-methods                                                         
#' @exportMethod setRegion                
setGeneric("setRegion", function(x,i,...) standardGeneric("setRegion"))


#' Method getRegions                                                                  
#' @name getRegions                                                             
#' @rdname getRegions-methods                                                         
#' @exportMethod getRegions                
setGeneric("getRegions", function(x) standardGeneric("getRegions"))

#' Method nCpG                                                                  
#' @name nCpG                                                              
#' @rdname nCpG-methods                                                         
#' @exportMethod nCpG                                                      
setGeneric("nCpG", function(x) standardGeneric("nCpG"))         


setMethod("getP", "Region",
          function(region){
              return(region@pVal)
        }
)
setMethod("getP", "RegionList",
          function(region){
              return(do.call(c,lapply(getRegions(region),getP)))
        }
)


setMethod("getPos", "Region",
          function(region){
              return(region@position)
        }
)


setMethod("getPos", "RegionList",
          function(region){
              return(do.call(c,lapply(getRegions(region),getPos)))
        }
)

setMethod("getT", "Region",
          function(region, index = NULL){
              if(is.null(index))
                  return(region@tValues)
              else
                  return(region@tValues[index])
        }
)

setMethod("getT", "RegionList",
          function(region, index = NULL){
              if(is.null(index))
                  return(do.call(c,lapply(getRegions(region),getT)))
              else
                  return(do.call(c,lapply(getRegions(region)[index],getT)))
        }
)

setMethod("[", "RegionList",
            function(x,i){
                x@regions[i]
            }
)

setMethod("[[", "RegionList",
            function(x,i){
                x@regions[[i]]
            }
)

setMethod("setRegion","RegionList",
          function(x,i,region){
              x@regions[[i]] <- region
              return(x)
            }
)

setMethod("getRegions","RegionList",
          function(x){
             return(x@regions)
            }
)

setMethod("setRegion","missing",
          function(x,i,...){
              x@regions[[i]] <- Region(...)
              return(x)
            }
)

setMethod("length", "Region",
            function(x){
                return(x@nCpG)
            }
)



setMethod("length", "RegionList",
          function(x){
              return(x@nRegions)
        }
)

## Same as length
setMethod("nCpG","Region",
          function(x){
              return(x@nCpG)
            }
)


setMethod("nCpG","RegionList",
          function(x){
              return(do.call(sum,lapply(getRegions(x),length)))
            }
)
                                               
setMethod("print", "Region",
          function(x,...){
            out <- paste("Region with ", x@nCpG, " nCpGs on chromosome", x@chromosome, 
                         ":", min(x@position), "-" ,max (x@position), " With P value ", x@pVal, sep = "")
       #     print(out) 

          return(print(out))
          }
)

setMethod("print", "RegionList",
          function(x){
           cat("|Genomic Coordinate \t\t| nCpGs | pVal |\n") 
           lapply(getRegions(x),show)
           return(invisible(0))
          }
)


setMethod("show", "Region",
          function(object){
          out <- paste("|Chr",object@chromosome, ":",min(object@position),"-",
                       max(object@position)," \t|", object@nCpG, "\t|", object@pVal,"|\n",sep="") 
      #    cat(out)
          return(invisible(cat(out)))
          }
)

setMethod("head","RegionList",
          function(x,n = 10L){
            if(length(x) > n)
                x <- RegionList(regions = x[1:n],nRegions = n)
            
          return(invisible(print(x)))
          }
)


setMethod("sort","RegionList",
          function(x, decreasing = FALSE){
            order       <- order(getP(x))
            x@regions   <- x@regions[order]    
            return(x)


          }
)

setMethod("names","Region",
          function(x){
              return(x@id)
        }
)

setMethod("names","RegionList",
          function(x){
            return(do.call(c,lapply(getRegions(x),names)))
          }
)


