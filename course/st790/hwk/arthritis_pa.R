#########################################################################
#
#   Homework 4, arthritis clinical trial
#   
#########################################################################

library(gee)
library(geepack)
library(ggplot2)

#  Read in data in long format
                      
thedat <- read.table("arthritis.dat")
thedat <- thedat[,c(1:4,6)]
colnames(thedat) <- c("id","trt","age","month","dscale")

#  Summarize the proportions at each week for each drug level

thedat.wide <- reshape(thedat,v.names="dscale",idvar="id",
                       timevar="month",direction="wide")

prop.placebo <- apply(thedat.wide[thedat.wide$trt==0,4:7],2,mean)
prop.auranofin <-  apply(thedat.wide[thedat.wide$trt==1,4:7],2,mean)

## > prop.placebo
##  dscale.0  dscale.2  dscale.4  dscale.6 
## 0.3517241 0.2965517 0.2758621 0.2551724 
## > prop.auranofin
##  dscale.0  dscale.2  dscale.4  dscale.6 
## 0.3241379 0.1241379 0.2000000 0.1448276 

#  The raw proportions suggest that the probability of severe arthritis
#  symptoms shows a slight decline for placebo and a more dramatic but
#  erratic decline for the active drug.  

#  Fits with gee() of unstructured and compound symmetry

un.gee <- gee(dscale ~ month + month:trt,id=id,family=binomial,
              corstr="unstructured",scale.fix=TRUE,scale.value=1,data=thedat)

## > summary(un.gee)

## Model:
##  Link:                      Logit 
##  Variance to Mean Relation: Binomial 
##  Correlation Structure:     Unstructured 

## Coefficients:
##               Estimate Naive S.E.   Naive z Robust S.E.
## (Intercept) -0.8211133 0.11344041 -7.238278  0.11407868
## month       -0.0545369 0.03239288 -1.683607  0.03311225
## month:trt   -0.1230894 0.04618776 -2.664980  0.04688966
##              Robust z
## (Intercept) -7.197780
## month       -1.647031
## month:trt   -2.625087

## Estimated Scale Parameter:  1
## Number of Iterations:  3

## Working Correlation
##           [,1]      [,2]      [,3]      [,4]
## [1,] 1.0000000 0.3003101 0.3589448 0.3240567
## [2,] 0.3003101 1.0000000 0.2507039 0.3308577
## [3,] 0.3589448 0.2507039 1.0000000 0.4782445
## [4,] 0.3240567 0.3308577 0.4782445 1.0000000

cs.gee <- gee(dscale ~ month + month:trt,id=id,family=binomial,
              corstr="exchangeable",data=thedat,scale.fix=TRUE,scale.value=1)

## >  summary(cs.gee)

## Model:
##  Link:                      Logit 
##  Variance to Mean Relation: Binomial 
##  Correlation Structure:     Exchangeable 

## Coefficients:
##                Estimate Naive S.E.   Naive z Robust S.E.
## (Intercept) -0.81499859 0.11517430 -7.076219  0.11407562
## month       -0.05560793 0.03106225 -1.790209  0.03360620
## month:trt   -0.11269741 0.04351464 -2.589873  0.04621873
##              Robust z
## (Intercept) -7.144371
## month       -1.654692
## month:trt   -2.438349

## Estimated Scale Parameter:  1
## Number of Iterations:  3

## Working Correlation
##           [,1]      [,2]      [,3]      [,4]
## [1,] 1.0000000 0.3414471 0.3414471 0.3414471
## [2,] 0.3414471 1.0000000 0.3414471 0.3414471
## [3,] 0.3414471 0.3414471 1.0000000 0.3414471
## [4,] 0.3414471 0.3414471 0.3414471 1.0000000

# Fit the same models using geepack

un.geeglm <- geeglm(dscale ~ month + month:trt,id=id,family=binomial,
                    corstr="unstructured",data=thedat,scale.fix=TRUE,std.err="san.se")

## > summary(un.geeglm)

## Call:
## geeglm(formula = dscale ~ month + month:trt, family = binomial, 
##     data = thedat, id = id, corstr = "unstructured", scale.fix = TRUE, 
##     std.err = "san.se")

##  Coefficients:
##             Estimate  Std.err   Wald Pr(>|W|)    
## (Intercept) -0.82121  0.11408 51.817 6.09e-13 ***
## month       -0.05460  0.03311  2.719  0.09915 .  
## month:trt   -0.12298  0.04689  6.878  0.00873 ** 
## ---
## Signif. codes:  
## 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Scale is fixed.

## Correlation: Structure = unstructured  Link = identity 

## Estimated Correlation Parameters:
##           Estimate Std.err
## alpha.1:2   0.3019 0.05921
## alpha.1:3   0.3608 0.06366
## alpha.1:4   0.3258 0.06495
## alpha.2:3   0.2520 0.06262
## alpha.2:4   0.3326 0.06271
## alpha.3:4   0.4808 0.06609
## Number of clusters:   290   Maximum cluster size: 4

cs.geeglm <- geeglm(dscale ~ month + month:trt,id=id,family=binomial,
              corstr="exchangeable",data=thedat,scale.fix=TRUE,std.err="san.se")

## > summary(cs.geeglm)

## Call:
## geeglm(formula = dscale ~ month + month:trt, family = binomial, 
##     data = thedat, id = id, corstr = "exchangeable", scale.fix = TRUE, 
##     std.err = "san.se")

##  Coefficients:
##             Estimate Std.err  Wald Pr(>|W|)    
## (Intercept)  -0.8150  0.1141 51.04  9.1e-13 ***
## month        -0.0555  0.0336  2.73    0.098 .  
## month:trt    -0.1129  0.0462  5.96    0.015 *  
## ---
## Signif. codes:  
## 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Scale is fixed.

## Correlation: Structure = exchangeable  Link = identity 

## Estimated Correlation Parameters:
##       Estimate Std.err
## alpha     0.34  0.0383
## Number of clusters:   290   Maximum cluster size: 4 

## It looks like a compound symmetric structure is a reasonable
## choice.  Refit the basic model in alternative parameterization

trtfac <- as.factor(thedat$trt)

cs.gee.alt <- gee(dscale ~ month:trtfac,id=id,family=binomial,
              corstr="exchangeable",data=thedat,scale.fix=TRUE,scale.value=1)

## > summary(cs.gee.alt)

## Model:
##  Link:                      Logit 
##  Variance to Mean Relation: Binomial 
##  Correlation Structure:     Exchangeable 

## Coefficients:
##               Estimate Naive S.E. Naive z Robust S.E.
## (Intercept)    -0.8150     0.1152   -7.08      0.1141
## month:trtfac0  -0.0556     0.0311   -1.79      0.0336
## month:trtfac1  -0.1683     0.0359   -4.69      0.0378
##               Robust z
## (Intercept)      -7.14
## month:trtfac0    -1.65
## month:trtfac1    -4.45

## Estimated Scale Parameter:  1
## Number of Iterations:  3

## Working Correlation
##       [,1]  [,2]  [,3]  [,4]
## [1,] 1.000 0.341 0.341 0.341
## [2,] 0.341 1.000 0.341 0.341
## [3,] 0.341 0.341 1.000 0.341
## [4,] 0.341 0.341 0.341 1.000

## Based on this fit, obtain estimates of the probability of severe
## symptoms at month 6 and standard errors by the delta method

beta <- cs.gee.alt$coefficients
robust <- cs.gee.alt$robust.variance

##  function to get estimate of form L beta and standard errors and
##  conf intervals

est.L <- function(L,beta,Sig,alpha){
    est <- L%*%beta
    se <- sqrt(L%*%Sig%*%t(L))
    crit.val <- abs(qnorm(alpha/2))
    CI <- c(est-crit.val*se,est+crit.val*se)
    result <- matrix(c(est,se,CI),1,4,byrow=TRUE)
    colnames(result) <- c("Estimate","Std Error","Lower CI", "Upper CI")
    return(result)
}

#  expit function

expit <- function(u){
    exp(u)/(1+exp(u))
        }

##  Estimate linear predictor at month 6 for each treatment

L.plac <- matrix(c(1,6,0),1,3,byrow=TRUE)
L.aura <- matrix(c(1,0,6),1,3,byrow=TRUE)

lp.plac <- est.L(L.plac,beta,robust,0.05)
lp.aura <- est.L(L.aura,beta,robust,0.05)

## > lp.plac
##      Estimate Std Error Lower CI Upper CI
## [1,]   -1.149    0.1895    -1.52  -0.7773
## > lp.aura
##      Estimate Std Error Lower CI Upper CI
## [1,]   -1.825     0.206   -2.229   -1.421

##  estimated probabilities -- we need to use the delta method to
##  get approximate standard errors

p.plac <- expit(lp.plac[1])
p.aura <- expit(lp.aura[1])
se.lp.plac <- lp.plac[2]
se.lp.aura <- lp.aura[2]
se.p.plac <- p.plac*(1-p.plac)*se.lp.plac
se.p.aura <- p.plac*(1-p.aura)*se.lp.aura

##  estimated probabilities and delta method standard errors

## > c(p.plac,se.p.plac)
## [1] 0.24074 0.03463
## > c(p.aura,se.p.aura)
## [1] 0.1389 0.0427

##  function to compute chi squared test for linear hypothesis

chi2.test <- function(L,beta,Sig){
    df <- nrow(L)
    T <- t(L%*%beta)%*%solve(L%*%Sig%*%t(L))%*%(L%*%beta)
    p.value <- pchisq(T,df,lower.tail=FALSE)
    return(list(T=T,p.value=p.value,T.pvalue=c(T,round(p.value,4))))
}

##  compare value of linear predictor at 6 months

L <- matrix(c(0,6,-6),1,3,byrow=TRUE)
at6months <- chi2.test(L,beta,robust)

## > at6months$T.pvalue
## [1] 5.9455 0.0148


cs.geeglm.alt <- geeglm(dscale ~  month:trtfac,id=id,family=binomial,
              corstr="exchangeable",data=thedat,scale.fix=TRUE,std.err="san.se")

## > summary(cs.geeglm.alt)

## Call:
## geeglm(formula = dscale ~ month:trtfac, family = binomial, data = thedat, 
##     id = id, corstr = "exchangeable", scale.fix = TRUE, std.err = "san.se")

##  Coefficients:
##               Estimate  Std.err   Wald Pr(>|W|)    
## (Intercept)   -0.81497  0.11408 51.039 9.06e-13 ***
## month:trtfac0 -0.05553  0.03360  2.732   0.0984 .  
## month:trtfac1 -0.16840  0.03783 19.819 8.51e-06 ***
## ---
## Signif. codes:  
## 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Scale is fixed.

## Correlation: Structure = exchangeable  Link = identity 

## Estimated Correlation Parameters:
##       Estimate Std.err
## alpha   0.3396 0.03833
## Number of clusters:   290   Maximum cluster size: 4 

#  Effect of age at baseline

cs.gee.baseline <- gee(dscale ~ age + month + month:trt,id=id,family=binomial,
                       corstr="exchangeable",data=thedat,scale.fix=TRUE,scale.value=1)

## > summary(cs.gee.baseline)

## Model:
##  Link:                      Logit 
##  Variance to Mean Relation: Binomial 
##  Correlation Structure:     Exchangeable 

## Coefficients:
##             Estimate Naive S.E. Naive z Robust S.E.
## (Intercept) -1.02117    0.46305  -2.205     0.47960
## age          0.00408    0.00885   0.461     0.00916
## month       -0.05582    0.03108  -1.796     0.03357
## month:trt   -0.11236    0.04353  -2.581     0.04618
##             Robust z
## (Intercept)   -2.129
## age            0.446
## month         -1.663
## month:trt     -2.433

## Estimated Scale Parameter:  1
## Number of Iterations:  3

## Working Correlation
##       [,1]  [,2]  [,3]  [,4]
## [1,] 1.000 0.341 0.341 0.341
## [2,] 0.341 1.000 0.341 0.341
## [3,] 0.341 0.341 1.000 0.341
## [4,] 0.341 0.341 0.341 1.000

#  Effect of age at baseline and in slopes

cs.gee.age <- gee(dscale ~ age + month + month:age + month:trt + month:trt:age,
                  id=id,family=binomial,
                  corstr="exchangeable",data=thedat,scale.fix=TRUE,scale.value=1)

## > summary(cs.gee.age)

## Model:
##  Link:                      Logit 
##  Variance to Mean Relation: Binomial 
##  Correlation Structure:     Exchangeable 

## Coefficients:
##                Estimate Naive S.E. Naive z Robust S.E.
## (Intercept)   -0.803272    0.53781  -1.494     0.57360
## age           -0.000208    0.01040  -0.020     0.01102
## month         -0.088834    0.14598  -0.609     0.16756
## age:month      0.000648    0.00280   0.232     0.00320
## month:trt     -0.263973    0.21431  -1.232     0.23270
## age:month:trt  0.002985    0.00408   0.731     0.00443
##               Robust z
## (Intercept)    -1.4004
## age            -0.0189
## month          -0.5302
## age:month       0.2027
## month:trt      -1.1344
## age:month:trt   0.6735

## Estimated Scale Parameter:  1
## Number of Iterations:  3

## Working Correlation
##       [,1]  [,2]  [,3]  [,4]
## [1,] 1.000 0.342 0.342 0.342
## [2,] 0.342 1.000 0.342 0.342
## [3,] 0.342 0.342 1.000 0.342
## [4,] 0.342 0.342 0.342 1.000




