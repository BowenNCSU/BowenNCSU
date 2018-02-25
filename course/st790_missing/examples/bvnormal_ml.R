################################################################
##
##  ML for bivariate normal data with missingness using 
##  the gls() function 
##
################################################################

library(nlme)
library(Matrix)
library(magic)

## Read in the data set

thedat <- read.table("bvnormal.R.dat")

##  Add an id variable for use with gls()

thedat <- cbind(seq(1,nrow(thedat),1),thedat)
colnames(thedat) <- c("id","y1","y2")

m <- nrow(thedat)   #  number of bivariate observations

##  The data set is in "wide" format with one record per individual,
##  so reconfigure to be in "long" format with one record per
##  observation as required for the gls() function

thedat.long <- reshape(thedat,varying=list(c("y1","y2")),v.names=c("y"),times=1:2,timevar="time",direction="long") 
thedat.long <- thedat.long[order(thedat.long$id),]

##  Maximum likelihood based on observed data
 
##  Need to have "time" be a classification variable to get a separate
##  mean for each component. Also need a "time variable" so that gls
##  can identify which observations below to which component in
##  constructing the covariance matrix

id1 <- as.factor(thedat.long$id)  
time1 <- as.factor(thedat.long$time)
thedat.long$comp <- as.numeric(factor(thedat.long$time))
  
mlfit <- gls(y ~ -1 + time1,correlation=corSymm(form = ~ comp| id1),
             weights=varIdent(form = ~ 1 | time1),method="ML",na.action=na.omit,data=thedat.long)

mu <- coef(mlfit)
Sigma <- getVarCov(mlfit)

## > mu
##   time11   time12 
##   4.991664 7.926341
## > Sigma
## Marginal variance covariance matrix
##         [,1]    [,2]
## [1,] 1.01620 0.51483
## [2,] 0.51483 1.09800  



  
