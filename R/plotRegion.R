#' Plot DMRs of type Region
#' @name plot.Region
#' @param x A Region object to be ploted. Can be subsetted from RegionList
#' @param ... Inherited from plot()
#' @aliases plot
#' @import ggplot2
#' @export
#' @return A plot object 
plot.Region <- function(x,...){                
     dat <- data.frame(                             
                        tVal = abs(getT(x)),        
                        pos  = getPos(x)         
          )                                 
    title <- paste("Region with", nCpG(x), "CpGs on ", range(x), "with P value", 
                   getP(x,n = 8))
           p <-ggplot(dat, aes_(x = pos, y = tVal)) +     
              geom_smooth(se = FALSE, span = 2) + ## roaling average 
              geom_point() +                                     
              labs(list(title = title, x = "Genomic Position",   
                                          y = "Test Statistic"))    
         return(p)                                              
                                                               
}
