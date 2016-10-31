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
#' @name initialize-Region
#' @rdname Region-class
#' @param .Object Empty
#' @param tValues A vector of test statistics
#' @param id A vector of id's for each test statistic
#' @param position A vector of position for each test statistc
#' @param chromosome An character describing the chromosome (1-22, X,Y)
#' @param pVal A numeric giving the estimated p-value for each significant DMR
#' @aliases Region-class,initialize
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
#' @rdname Region-class
#' @param ... Parameters to pased to new("Regions",...)
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
#' @name initialize-RegionList 
#' @rdname RegionList-class
#' @param .Object Empty
#' @param nRegions The number of regions to be included in the RegionList object
#' @param regions An optional list of regions to be assigned to a RegionList
#' @aliases RegionList-class,initialize
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
#' @param ... Parameters to be past to new("RegionList",...)
RegionList <- function(...){ 
                new("RegionList",...)         
}
 
## Set generic functions for Rt() and St() when 
## One for class "list" and one for class "Region"
 
#' Method Rt 
#' @name Rt
#' @rdname Rt-methods
#' @exportMethod Rt
#' @param region Object of type Region or RegionList
#' @param windowThreshold Vector of window thresholds
#' @param windowSize Vector of window sizes to be tested on regions
setGeneric("Rt", function(region,windowThreshold,windowSize) standardGeneric("Rt"))


#' Method St 
#' @name St
#' @rdname St-methods
#' @exportMethod St
#' @param region Object of type Region or RegionList
#' @param windowThreshold Vector of window thresholds
#' @param windowSize Vector of window sizes to be tested on regions
setGeneric("St", function(region,windowThreshold,windowSize) standardGeneric("St"))


#' Method getP 
#' @name getP
#' @rdname getP
#' @exportMethod getP
#' @param region An object of type Region or RegionList
setGeneric("getP", function(region) standardGeneric("getP"))

#' Method getPos 
#' @name getPos
#' @rdname getPos
#' @exportMethod getPos
#' @param region An opbject of type Region or RegionList
setGeneric("getPos", function(region) standardGeneric("getPos"))



#' Method getT 
#' @name getT
#' @rdname getT
#' @exportMethod getT
#' @param region An opbject of type Region or RegionList
#' @param ... Index 
setGeneric("getT", function(region,...) standardGeneric("getT"))


#' Method setRegion                                                                  
#' @name setRegion                                                             
#' @rdname setRegion                                                         
#' @exportMethod setRegion                
#' @param x A region
#' @param i an index
setGeneric("setRegion", function(x,i,...) standardGeneric("setRegion"))


#' Method getRegions                                                                  
#' @name getRegions                                                             
#' @rdname getRegions
#' @param x A RegionList object
#' @exportMethod getRegions                
setGeneric("getRegions", function(x) standardGeneric("getRegions"))

#' Method nCpG                                                                  
#' @name nCpG                                                              
#' @rdname nCpG
#' @exportMethod nCpG                                                     
#' @param x An opbject of type Region or RegionList
#' @aliases length
setGeneric("nCpG", function(x) standardGeneric("nCpG"))         

#' Get p-values for a region
#' @rdname getP
setMethod("getP", "Region",
          function(region){
              return(region@pVal)
        }
)

#' Get p-values for a list of regions (RegionList)
#' @rdname getP
setMethod("getP", "RegionList",
          function(region){
              return(do.call(c,lapply(getRegions(region),getP)))
        }
)

#' Get the chromosomal coordinates for a Region
#' @rdname getPos
setMethod("getPos", "Region",
          function(region){
              return(region@position)
        }
)

#' Get the chromosomal coordinates for a list of regions in a RegionList object
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
#' @name [<-
#' @rdname index
#' @param x An object of type RegionList
#' @param i Index, which region to extract
#' @param j Not used
#' @usage \S4method{[}{RegionList,ANY}(x, i)
#' @aliases [
setMethod("[<-", signature(x = "RegionList", i = "ANY"),
            function(x,i){
                x@regions[i]
            }
)

#' Get a region from a RegionList object
#' @name [[
#' @aliases [[, RegionList-method
#' @rdname index_list 
#' @param x An object of type RegionList                                        
#' @param i Index, which region to extract
#' @param j Not Used
setMethod("[[", signature(x = "RegionList", i = "ANY", j = "ANY"),
            function(x,i,j){
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
setMethod("length", "Region",
            function(x){
                return(x@nCpG)
            }
)

#' Get the number of regions in a RegionList
#' @rdname length
setMethod("length", "RegionList",
          function(x){
              return(x@nRegions)
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
#' @param n The number of regions to be printed when the RegionList is longer than n
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
setMethod("sort","RegionList",
          function(x, decreasing = FALSE){
            order       <- order(getP(x))
            x@regions   <- x@regions[order]    
            return(x)


          }
)

#' Get the names of all probes within a region
#' @rdname names
setMethod("names","Region",
          function(x){
              return(x@id)
        }
)

#' Get the names of all probes in a study
#' @rdname names
#' @param x An object of type RegionList
setMethod("names","RegionList",
          function(x){
            return(do.call(c,lapply(getRegions(x),names)))
          }
)

