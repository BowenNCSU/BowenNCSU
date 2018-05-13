#########################################################################
#
#   Homework 4, Problem 3, Leprosy clinical trial
#   
#########################################################################

library(gee)
library(geepack)
library(ggplot2)

## Read in data in wide format

thedat.wide <- read.table("leprosy.dat")
colnames(thedat.wide) <- c("id","trt","pre","post")

#  Total number of individuals

m <- nrow(thedat.wide)

## Means and variances at baseline and follow-up for each treatment

mean.placebo <- apply(thedat.wide[thedat.wide$trt==0,3:4],2,mean)
var.placebo <- apply(thedat.wide[thedat.wide$trt==0,3:4],2,var)

mean.anti1 <- apply(thedat.wide[thedat.wide$trt==1,3:4],2,mean)
var.anti1 <- apply(thedat.wide[thedat.wide$trt==1,3:4],2,var)

mean.anti2 <- apply(thedat.wide[thedat.wide$trt==2,3:4],2,mean)
var.anti2 <- apply(thedat.wide[thedat.wide$trt==2,3:4],2,var)

## > rbind(mean.placebo,var.placebo)
##                   pre     post
## mean.placebo 12.90000 12.30000
## var.placebo  15.65556 51.12222
## > rbind(mean.anti1,var.anti1)
##                 pre     post
## mean.anti1  9.30000  5.30000
## var.anti1  22.67778 21.56667
## > rbind(mean.anti2,var.anti2)
##                 pre     post
## mean.anti2 10.00000  6.10000
## var.anti2  27.55556 37.87778

##  Clearly, variance is not equal to mean; there is apparent
##  overdispersion (except for the pretest values for placebo; note
##  also that the pretest mean was highest for placebo, despite
##  randomization

m <- nrow(thedat.wide)

## Reconfigure in long format

thedat <- reshape(thedat.wide,varying=c("pre","post"), 
            v.names="count",idvar="id",times=c(0,1),timevar="time",direction="long")

thedat <- thedat[order(thedat$id),]

##  Fit the model -- because there are only 2 time points, any
##  correlation model will give the same result; we use compound
##  symmetry

##  R will treat placebo (the "lowest" factor value) as the reference
##  treatment as in our desired model by default

tx <- as.factor(thedat$trt)

cs.gee <- gee(count ~ time + time:tx, id=id, family=poisson(link="log"),
              corstr="exchangeable",data=thedat,scale.fix=FALSE)

## This is the same result obtained in SAS using PROC GENMOD and PROC GEE

## > summary(cs.gee)

## Model:
##  Link:                      Logarithm 
##  Variance to Mean Relation: Poisson 
##  Correlation Structure:     Exchangeable 

## Call:
## gee(formula = count ~ time + time:tx, id = id, data = thedat, 
##     family = poisson(link = "log"), corstr = "exchangeable", 
##     scale.fix = FALSE)

## Coefficients:
##                Estimate Naive S.E.    Naive z Robust S.E.    Robust z
## (Intercept)  2.37335416  0.1035324 22.9237797  0.08013786 29.61588978
## time        -0.01382285  0.1110676 -0.1244544  0.15733975 -0.08785353
## time:tx1    -0.54058998  0.1817710 -2.9740158  0.21857634 -2.47323187
## time:tx2    -0.47906560  0.1779348 -2.6923657  0.22786326 -2.10242579

## Estimated Scale Parameter:  3.451505
## Number of Iterations:  5

## Working Correlation
##           [,1]      [,2]
## [1,] 1.0000000 0.7965695
## [2,] 0.7965695 1.0000000

##  This result and that using geeglm() with compound symmetry below
##  are not the same as the one obtained using gee().  Instead they
##  are the same as that using PROC GLIMMIX, so solving the quadratic
##  estimating equation to estimate alpha and sigma^2.  It appears
##  that gee() with unstructured and geeglm() with CS are using the
##  quadratic estimating equation instead of the simple moment
##  estimator.

un.gee <- gee(count ~ time + time:tx, id=id, family=poisson(link="log"),
              corstr="unstructured",data=thedat,scale.fix=FALSE)

## > summary(un.gee)

## Model:
##  Link:                      Logarithm 
##  Variance to Mean Relation: Poisson 
##  Correlation Structure:     Unstructured 

## Call:
## gee(formula = count ~ time + time:tx, id = id, data = thedat, 
##     family = poisson(link = "log"), corstr = "unstructured", 
##     scale.fix = FALSE)

## Coefficients:
##                 Estimate Naive S.E.     Naive z
## (Intercept)  2.373354164  0.1034103 22.95084777
## time        -0.002875909  0.1239320 -0.02320554
## time:tx1    -0.562568679  0.2023204 -2.78058272
## time:tx2    -0.495284414  0.1977550 -2.50453486
##             Robust S.E.    Robust z
## (Intercept)  0.08013786 29.61588978
## time         0.15700486 -0.01831732
## time:tx1     0.22198194 -2.53429933
## time:tx2     0.23420052 -2.11478786

## Estimated Scale Parameter:  3.443369
## Number of Iterations:  5

## Working Correlation
##           [,1]      [,2]
## [1,] 1.0000000 0.7383648
## [2,] 0.7383648 1.0000000

## We do not request empirical robust standard errors here, as this
## covariance structure is the only one possible.

cs.geeglm <- geeglm(count ~ time + time:tx,id=id,family=poisson,waves=time,
                    corstr="exch",data=thedat,scale.fix=FALSE)

## > summary(cs.geeglm)

## Call:
## geeglm(formula = count ~ time + time:tx, family = poisson, data = thedat, 
##     id = id, waves = time, corstr = "exch", scale.fix = FALSE)

##  Coefficients:
##              Estimate   Std.err    Wald Pr(>|W|)    
## (Intercept)  2.373354  0.080138 877.101   <2e-16 ***
## time        -0.002877  0.157005   0.000   0.9854    
## time:tx1    -0.562571  0.221983   6.423   0.0113 *  
## time:tx2    -0.495284  0.234201   4.472   0.0344 *  
## ---
## Signif. codes:  
## 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Estimated Scale Parameters:
##             Estimate Std.err
## (Intercept)    3.214  0.4998

## Correlation: Structure = exchangeable  Link = identity 

## Estimated Correlation Parameters:
##       Estimate Std.err
## alpha   0.7384 0.08149
## Number of clusters:   30   Maximum cluster size: 2 

## Using AR(1) also gives the same result

ar1.gee <- gee(count ~ time + time:tx, id=id, family=poisson(link="log"),
+ corstr="AR-M",Mv=1,data=thedat,scale.fix=FALSE)

## > summary(ar1.gee)

##  GEE:  GENERALIZED LINEAR MODELS FOR DEPENDENT DATA
##  gee S-function, version 4.13 modified 98/01/27 (1998) 

## Call:
## gee(formula = count ~ time + time:tx, id = id, data = thedat, 
##     family = poisson(link = "log"), corstr = "AR-M", Mv = 1, 
##     scale.fix = FALSE)


## Coefficients:
##              Estimate Naive S.E.  Naive z Robust S.E.
## (Intercept)  2.373354     0.1034 22.95085     0.08014
## time        -0.002876     0.1239 -0.02321     0.15700
## time:tx1    -0.562569     0.2023 -2.78058     0.22198
## time:tx2    -0.495284     0.1978 -2.50453     0.23420
##             Robust z
## (Intercept) 29.61589
## time        -0.01832
## time:tx1    -2.53430
## time:tx2    -2.11479

## Estimated Scale Parameter:  3.443
## Number of Iterations:  5

## Working Correlation
##        [,1]   [,2]
## [1,] 1.0000 0.7384
## [2,] 0.7384 1.0000
