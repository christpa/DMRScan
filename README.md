# DMRScan
R package for detection of differentially methylated regions with adjustment for multiple testing

# Important Info
## Authors
CM Page, L Vos, BK Andreasen

## Package structure
The function `DMRScan()` requires three key inputs; 
  - a vector of test statistics for each CpG site, 
  - the corresponding positions (chromosome and bp position). 

Additional inputs are the 
  - maximum allowed distance within a region
  - the minimum number of CpGs within a region.

<<<<<<< HEAD
### Installation 
For a stable release from Bioconductor use
```R
source("https://www.bioconductor.org/biocLite.R")
biocLite("DMRScan")
```
For the developmental version from Github, use 
```R
install_github("christmp/DMRScan")
```
Requred R-packages for DMRScan are:
* Matrix 
* MASS 
* RcppRoll 
* GenomicRanges
* IRanges
* methods
* mvtnorm
* stats
* parallel

=======
>>>>>>> upstream/master
### Usage
Please see the vignette.

# Citation
To cite the DMRScan package in publications use:

<<<<<<< HEAD
*Assessing Genome-Wide Significance for the Detection of Differentially Methylated Regions*. Page CM, Vos L, Rounge TB, Harbo HF, and Andreassen BK _Submitted_., (2017)
=======
*Assessing Genome-Wide Significance for the Detection of Differentially Methylated Regions*. Page CM, Vos L, Rounge TB, Harbo HF, and Andreassen BK _Submitted to BMC Bioinformatics_., (2016)
>>>>>>> upstream/master

A BibTeX entry for LaTeX users is
```BibTeX
@Article{DMRScanRpackage,
title = {Assessing Genome-Wide Significance for the Detection of Differentially Methylated Regions},
author = {Christian M Page and Linda Vos and Trine B Rounge and Hanne F Harbo and Bettina K Andreassen},
<<<<<<< HEAD
year = {2017}}
=======
journal = {BMC Bioinformatics (Submitted)},
year = {2016}}
>>>>>>> upstream/master
```
End of README
