#' @export
plot.Region <- function(x,...){                                                 
     dat <- data.frame(                                                  
                        tVal = abs(getT(x)),                              
                        pos  = getPos(x)                                  
          )                                                               
    title <- strsplit(print(x),"With")[[1]][1]                            
     
    if("ggplot2" %in% installed.packages()[,1]){
         require(ggplot2)
            p <- ggplot(dat, aes(x = pos, y = tVal)) +                            
               geom_smooth(se = FALSE, span = 2) + ## roaling average!                    
               geom_point() +                                                   
               labs(list(title = title, x = "Genomic Position",                 
                                           y = "Test Statistic"))               
          return(p)                                                             
    }else{
        plot(x = dat$pos, y = dat$tVal, type = "p", 
            main = title, xlab = "Genomic Position", ylab = "Test Statistic")
        lines(smooth.spline(x = dat$pos, y = dat$tVal))
    }                                                               
}
