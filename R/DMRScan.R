#' DMRScan: An R-package for identification of Differentially Metylated Regions 
#' @author Christian Page, \email{page.ntnu@gmail.com}
#' @references Not Published yet (Under revision) 
#' @aliases DMRScan_package
#' @keywords DMR, DMRScan
#' @docType package
#' @name DMRScan_package
#' @param observations An object of type RegionList                             
#' @param windowSize A sequence of windowSizes for the slidingWindow,           
#' must be an integer                                                           
#' @param windowThreshold Optional argument with corresponding cut-off for      
#' each window. Will be estimated if not supplied.                              
#' @param ... Optional arguments to be pased to \code{\link{estimateWindowThreshold}},     
#' if no grid is specified.                                                     
#' @return An object of type \code{\link{RegionList}} with signficantly differentially        
# methylated regions                                                            
#' @examples                                                                    
#' ## nProbeoad methylation data from chromosome 22                             
#' data(DMRScan.methylationData)                                                
#' ## nProbeoad phenotype (end-point for methylation data)                      
#' data(DMRScan.phenotypes)                                                     
#'                                                                              
#' ## Test for an association between phenotype and Methylation                 
#' test.statistics <- apply(DMRScan.methylationData,1,function(x,y)             
#'   summary(glm(y ~ x, family = binomial(link = "logit")))$coefficients[2,3],  
#'                                                    y = DMRScan.phenotypes)   
#' ## Set chromosomal position to each test-statistic                           
#' positions <- data.frame(matrix(as.integer(unlist(strsplit(names(test.statistics), split="chr|[.]"))), ncol = 3, byrow = TRUE))[,-1]
#' ## Set clustering features                                                   
#' min.cpg <- 4  ## Minimum number of CpGs in a tested cluster                  
#' ## Maxium distance (in base-pairs) within a cluster                          
#' ## before it is broken up into two seperate cluster                          
#' max.gap <- 750                                                               
#'                                                                              
#                                                                               
#' ## Identify all clusters, and generate a list for each cluster               
#' regions <- makeCpGregions(observations = test.statistics,                    
#'                           chr = positions[,1], pos = positions[,2],          
#'                           maxGap = max.gap, minCpG = min.cpg)                
#' ## Number of CpGs in the slidingWindows, can be either a single number       
#' ## or a sequence of windowSizes                                              
#' windowSizes <- 3:7                                                           
#' nCpG        <- nCpG(regions) ## Number of CpGs to be tested                  
#'                                                                              
#' # Estimate the windowThreshold, based on the number of CpGs and windowSizes  
#' windowThresholds <- estimateWindowThreshold(nProbe = nCpG,                   
#'                windowSize = windowSizes, method = "sampling", mcmc = 10000)  
#' ## Run the slidingWindow                                                     
#' DMRScanResults   <- dmrscan(observations = regions,                          
#'                             windowSize = windowSizes,                        
#'                             windowThreshold = windowThresholds)              
#' ## Print the result                                                          
#' print(DMRScanResults)                                                        
NULL
