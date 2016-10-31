#' Plot DMRs of type Region
#' @name plot.Region
#' @param x A Region object to be ploted. Can be subsetted from RegionList
#' @aliases plot
#' @import methods 
#' @param ... Inherited from plot()
#' @export
plot.Region <- function(x,...){                                                 
     dat <- data.frame(                                                  
                        tVal = abs(getT(x)),                              
                        pos  = getPos(x)                                  
          )                                                               
    title <- strsplit(print(x),"With")[[1]][1]                            
     
    if("ggplot2" %in% installed.packages()[,1]){
           p <- ggplot2::ggplot(dat, aes(x = pos, y = tVal)) +                            
              ggplot2::geom_smooth(se = FALSE, span = 2) + ## roaling average                    
              ggplot2::geom_point() +                                                   
              ggplot2::labs(list(title = title, x = "Genomic Position",                 
                                          y = "Test Statistic"))               
         return(p)                                                             
    }else{
        plot(x = dat$pos, y = dat$tVal, type = "p", 
            main = title, xlab = "Genomic Position", ylab = "Test Statistic")
        lines(smooth.spline(x = dat$pos, y = dat$tVal))
    }                                                               
}
