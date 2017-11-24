na.omitt.dmrscan <- function(x){
        if(is.na(x[1,1])){
            out <- NULL
        }else{
            d1 <- max(which(!is.na(x[col(x) == row(x)])))    
            if(d1 < 2){
                out <- NULL
            }else{
                out <- x[1:d1,1:d1]
            }
        }
            return(out)
}
