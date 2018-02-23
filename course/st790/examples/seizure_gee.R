#########################################################################
#
#   CHAPTER 8, EXAMPLE 5, Epileptic Seizure Study 
# 
#   Fitting population-averaged models using the gee() function and 
#   the geeglm() function in geepack.  You will need to install.packages("gee")
#   and install.packages("geepack").
#
#   The syntax for each is similar.  A general call looks like
#
#   fit.object <- gee(model formula, id, corstr, family, data)
#
#   model formula is a usual R-type model specification for a generalized
#   linear model -- see the documentation for glm()
#
#   corstr is a specification for the working correlation matrix
#
#   family is a specification for the scaled exponential family that
#   would be relevant under independence; the canonical link is the default
#
#   Both of these implementations use moment methods to estimate the  
#   parameters in the working correlation structure.  It is not straightforward
#   to trick another R function into fitting the model with the covariance
#   parameters estimated using the quadratic estimating equation with
#   the Gaussian working structure.  This is possible using the gnls()
#   function in the nlme package, but requires some work on the part
#   of the user, so we don't show this.  
#
#########################################################################

library(gee)
library(geepack)

#  Read in the data

thedat <- read.table("seize.dat")
colnames(thedat) <- c("subj","seize","visit","trt","base","age")

#  Create other covariates

logage <- log(thedat$age)
o <- 8*(thedat$visit==0)+2*(thedat$visit>0)
logo <- log(o)
v <- as.numeric(thedat$visit>0)
v4 <- as.numeric(thedat$visit==4)

#  Basic models -- the unstructured fit using gee() does not converge
#  even with maxiter set to be much larger than the default

un.gee <- gee(seize ~ v*trt + offset(logo),id=subj,family=poisson,
              corstr="unstructured",data=thedat,maxiter=100)

cs.gee <- gee(seize ~ v*trt + offset(logo),id=subj,family=poisson,
              corstr="exchangeable",data=thedat)
## > summary(cs.gee)
##  GEE:  GENERALIZED LINEAR MODELS FOR DEPENDENT DATA
##  gee S-function, version 4.13 modified 98/01/27 (1998) 

## Model:
##  Link:                      Logarithm 
##  Variance to Mean Relation: Poisson 
##  Correlation Structure:     Exchangeable 

## Call:
## gee(formula = seize ~ v * trt + offset(logo), id = subj, data = thedat, 
##     family = poisson, corstr = "exchangeable")

## Summary of Residuals:
##        Min         1Q     Median         3Q        Max 
##  -4.299107  -1.299107   2.020161  10.374640 147.048387 


## Coefficients:
##                Estimate Naive S.E.    Naive z Robust S.E.   Robust z
## (Intercept)  1.34760922  0.1511851  8.9136359   0.1573571  8.5640166
## v            0.11079814  0.1547038  0.7161956   0.1160997  0.9543358
## trt          0.02651461  0.2072721  0.1279217   0.2218539  0.1195138
## v:trt       -0.10368067  0.2199500 -0.4713830   0.2136100 -0.4853736

## Estimated Scale Parameter:  19.70269
## Number of Iterations:  1

## Working Correlation
##          [,1]     [,2]     [,3]     [,4]     [,5]
## [1,] 1.000000 0.771588 0.771588 0.771588 0.771588
## [2,] 0.771588 1.000000 0.771588 0.771588 0.771588
## [3,] 0.771588 0.771588 1.000000 0.771588 0.771588
## [4,] 0.771588 0.771588 0.771588 1.000000 0.771588
## [5,] 0.771588 0.771588 0.771588 0.771588 1.000000

ar1.gee <- gee(seize ~ v*trt + offset(logo),id=subj,family=poisson,
              corstr="AR-M",Mv=1,data=thedat)
## > summary(ar1.gee)
## Coefficients:
##                Estimate Naive S.E.     Naive z Robust S.E.    Robust z
## (Intercept)  1.31202191  0.1429432  9.17862144   0.1618546  8.10617467
## v            0.15132391  0.1685943  0.89756232   0.1116287  1.35560051
## trt          0.01885819  0.1963247  0.09605613   0.2120029  0.08895252
## v:trt       -0.12824307  0.2409461 -0.53224803   0.2601825 -0.49289653

## Estimated Scale Parameter:  20.16395
## Number of Iterations:  3

## Working Correlation
##           [,1]      [,2]      [,3]      [,4]      [,5]
## [1,] 1.0000000 0.8102820 0.6565569 0.5319962 0.4310670
## [2,] 0.8102820 1.0000000 0.8102820 0.6565569 0.5319962
## [3,] 0.6565569 0.8102820 1.0000000 0.8102820 0.6565569
## [4,] 0.5319962 0.6565569 0.8102820 1.0000000 0.8102820
## [5,] 0.4310670 0.5319962 0.6565569 0.8102820 1.0000000

#  This fit is unreliable; note that one of the correlations is estimated
#  as being > 1

un.geeglm <- geeglm(seize ~ v*trt,id=subj,family=poisson("log"),
             offset=logo, corstr="unstructured",data=thedat)
## > summary(un.geeglm)

## Call:
## geeglm(formula = seize ~ v * trt, family = poisson("log"), data = thedat, 
##     offset = logo, id = subj, corstr = "unstructured")

##  Coefficients:
##             Estimate  Std.err  Wald Pr(>|W|)
## (Intercept)  1.18058  0.86813 1.849    0.174
## v            0.16106  0.31725 0.258    0.612
## trt          0.05042  1.05621 0.002    0.962
## v:trt       -0.11250  0.50257 0.050    0.823

## Estimated Scale Parameters:
##             Estimate Std.err
## (Intercept)    22.09   11.74

## Correlation: Structure = unstructured  Link = identity 

## Estimated Correlation Parameters:
##           Estimate Std.err
## alpha.1:2   1.0261  0.5257
## alpha.1:3   0.7339  0.3602
## alpha.1:4   0.8244  0.3413
## alpha.1:5   0.6900  0.3374
## alpha.2:3   0.8120  0.2945
## alpha.2:4   0.9541  0.2972
## alpha.2:5   0.7793  0.2704
## alpha.3:4   0.7080  0.1415
## alpha.3:5   0.5350  0.1486
## alpha.4:5   0.6813  0.1279
## Number of clusters:   59   Maximum cluster size: 5 

cs.geeglm <- geeglm(seize ~ v*trt,id=subj,family=poisson("log"),
             offset=logo, corstr="exchangeable",data=thedat)
## > summary(cs.geeglm)
## Call:
## geeglm(formula = seize ~ v * trt, family = poisson("log"), data = thedat, 
##     offset = logo, id = subj, corstr = "exchangeable")

##  Coefficients:
##             Estimate Std.err  Wald Pr(>|W|)    
## (Intercept)   1.3476  0.1574 73.34   <2e-16 ***
## v             0.1108  0.1161  0.91     0.34    
## trt           0.0265  0.2219  0.01     0.90    
## v:trt        -0.1037  0.2136  0.24     0.63    
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Estimated Scale Parameters:
##             Estimate Std.err
## (Intercept)     19.4     8.7

## Correlation: Structure = exchangeable  Link = identity 

## Estimated Correlation Parameters:
##       Estimate Std.err
## alpha    0.777  0.0753
## Number of clusters:   59   Maximum cluster size: 5

ar1.geeglm <- geeglm(seize ~ v*trt,id=subj,family=poisson("log"),
             offset=logo, corstr="ar1",data=thedat)
## > summary(ar1.geeglm)

## Call:
## geeglm(formula = seize ~ v * trt, family = poisson("log"), data = thedat, 
##     offset = logo, id = subj, corstr = "ar1")

##  Coefficients:
##             Estimate Std.err  Wald Pr(>|W|)    
## (Intercept)   1.3089  0.1622 65.14  6.7e-16 ***
## v             0.1554  0.1140  1.86     0.17    
## trt           0.0153  0.2118  0.01     0.94    
## v:trt        -0.1306  0.2676  0.24     0.63    
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Estimated Scale Parameters:
##             Estimate Std.err
## (Intercept)       20    8.96

## Correlation: Structure = ar1  Link = identity 

## Estimated Correlation Parameters:
##       Estimate Std.err
## alpha    0.893  0.0388
## Number of clusters:   59   Maximum cluster size: 5 

#  Fancier model with working compound symmetry

cs.fancy.gee <- gee(seize ~ logage + v4 + v*trt + v4:trt + offset(logo),
               id=subj,family=poisson,corstr="exchangeable",data=thedat)
## > summary(cs.fancy.gee)
## Coefficients:
##             Estimate Naive S.E. Naive z Robust S.E. Robust z
## (Intercept)   4.4254      1.357   3.262       1.375    3.218
## logage       -0.9250      0.408  -2.269       0.416   -2.224
## v4           -0.1009      0.160  -0.631       0.120   -0.839
## v             0.1351      0.153   0.881       0.133    1.015
## trt          -0.0264      0.204  -0.130       0.214   -0.124
## v:trt        -0.0769      0.216  -0.356       0.225   -0.342
## v4:trt       -0.1210      0.230  -0.526       0.130   -0.933

## Estimated Scale Parameter:  18.7
## Number of Iterations:  3

## Working Correlation
##       [,1]  [,2]  [,3]  [,4]  [,5]
## [1,] 1.000 0.767 0.767 0.767 0.767
## [2,] 0.767 1.000 0.767 0.767 0.767
## [3,] 0.767 0.767 1.000 0.767 0.767
## [4,] 0.767 0.767 0.767 1.000 0.767
## [5,] 0.767 0.767 0.767 0.767 1.000

cs.fancy.geeglm <- geeglm(seize ~ logage+v4+v*trt+v4:trt,id=subj,family=poisson("log"),
             offset=logo, corstr="exchangeable",data=thedat)
## > summary(cs.fancy.geeglm)
## Call:
## geeglm(formula = seize ~ logage + v4 + v * trt + v4:trt, family = poisson("log"), 
##     data = thedat, offset = logo, id = subj, corstr = "exchangeable")

##  Coefficients:
##             Estimate Std.err  Wald Pr(>|W|)   
## (Intercept)   4.4503  1.3775 10.44   0.0012 **
## logage       -0.9325  0.4165  5.01   0.0252 * 
## v4           -0.1009  0.1203  0.70   0.4017   
## v             0.1351  0.1330  1.03   0.3099   
## trt          -0.0269  0.2139  0.02   0.9000   
## v:trt        -0.0769  0.2249  0.12   0.7323   
## v4:trt       -0.1210  0.1297  0.87   0.3507   
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Estimated Scale Parameters:
##             Estimate Std.err
## (Intercept)     18.3    7.15

## Correlation: Structure = exchangeable  Link = identity 

## Estimated Correlation Parameters:
##       Estimate Std.err
## alpha    0.777  0.0766
## Number of clusters:   59   Maximum cluster size: 5 


