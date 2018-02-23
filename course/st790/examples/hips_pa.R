#########################################################################
#
#   CHAPTER 5, Hip Replacement Study
# 
#   Fitting population-averaged models uing the gls() function in nlme
#
#   See the program for the dental study program for details.
#                                        
#   As discussed in the notes, the missing values are almost entirely at
#   week 2, suggesting that something that had nothing to do with
#   the evolution of haematocrit or the gender or age of the patients, so
#  that it may be reasonable to assume MCAR missingness.  However, to be
#  safe, we use maximum likelihood as discussed below.
#
#   As noted in that program, gls() does not compute standard errors
#   correctly, nor does it compute robust standard errors.  We include
#   the function robust.cov() so we can get correct SEs for the
#   preferred model.  As in Section 5.6, we use model-ased SEs to be safe.
#
#########################################################################

library(nlme)
library(Matrix)
library(magic)

#  u <- model object from fit, m <- total # individuals

robust.cov <- function(u,m){
  form <-  formula(u)
  mf <- model.frame(form,getData(u))
  Xmat <- model.matrix(form,mf)
  Vlist <-  as.list(1:m)
  for (i in 1:m){
    Vlist[[i]] <- getVarCov(u,individual=i,type="marginal")
  }
  V <- Reduce(adiag,Vlist)
  Vinv <- solve(V)
  Sig.model <- solve(t(Xmat)%*%Vinv%*%Xmat)
  resid <- diag(residuals(u,type="response"))
  ones.list <- lapply(Vlist,FUN=function(u){matrix(1,nrow(u),ncol(u))})
  Ones <- Reduce(adiag,ones.list)
  meat <- t(Xmat)%*%Vinv%*%resid%*%Ones%*%resid%*%Vinv%*%Xmat
  Sig.robust <- Sig.model%*%meat%*%Sig.model
  se.robust <- sqrt(diag(Sig.robust))
  se.model <- sqrt(diag(Sig.model))
  return(list(Sig.model=Sig.model,se.model=se.model,Sig.robust=Sig.robust,se.robust=se.robust))
}

#  Read in the data -- they are in the "long" form required by gls()

thedat <- read.table("hips.dat")
colnames(thedat) <- c("pid","sex","age","week","h")

m <- max(thedat$pid)

#  create factors for use in gls() model and covariance

id <- factor(thedat$pid)
gender <- factor(thedat$sex)

week2 <- thedat$week^2

#  Create a time variable that will be used in the gls() calls to
#  align the times appropriately to account for the missing observations.
#  The times must be consecutive integer values, apparently starting
#  with 1.  The following assigns 1,2,3,4 to weeks 0,1,2,3

thedat$time <- as.numeric(factor(thedat$week))

#  Fit the basic quadratic model using ML.  To align the times properly
#  in the construction of the correlation matrix, we use the "time" variable
#  rather than 1 in the correlation and weight specifications

#  Common unstructured

hips.un <- gls(h ~ -1 + gender + week:gender + week2:gender,
               correlation=corSymm(form = ~ time | id),
               weights = varIdent(form = ~ 1 | week),
               data=thedat,method="ML")
beta.un <- coef(hips.un)
sebeta.un.incorrect <- summary(hips.un)$tTable[,"Std.Error"]
sebeta.un <- robust.cov(hips.un,m)$se.model

# patient 10 was observed at all 4 weeks; patient 1 is missing
# week 2

Gamma.un.1 <- corMatrix(hips.un$modelStruct$corStruct)[[1]]   
Gamma.un.10 <- corMatrix(hips.un$modelStruct$corStruct)[[10]]   
V.un.1 <- getVarCov(hips.un, individual=1)
V.un.10 <- getVarCov(hips.un, individual=10)
vars.un <- diag(V.un.10)

## > rbind(beta.un,sebeta.un)
##              gender0   gender1 gender0:week gender1:week gender0:week2
## beta.un   39.2553414 42.703366   -11.179391   -15.828690     2.8691480
## sebeta.un  0.9036174  1.023986     1.729341     1.942935     0.5357329
##           gender1:week2
## beta.un       4.2269624
## sebeta.un     0.6048256
# incorrect standard errors
## > sebeta.un.incorrect
##       gender0       gender1  gender0:week  gender1:week gender0:week2 
##     0.9323108     1.0565013     1.7842542     2.0046309     0.5527445 
## gender1:week2 
##     0.6240312 
## > Gamma.un.1
##          [,1]         [,2]         [,3]
## [1,] 1.0000000  0.258678035  0.234156597
## [2,] 0.2586780  1.000000000 -0.008103682
## [3,] 0.2341566 -0.008103682  1.000000000
## > V.un.1
## Marginal variance covariance matrix
##         [,1]     [,2]     [,3]
## [1,] 16.2700  4.10510  3.95860
## [2,]  4.1051 15.47900 -0.13363
## [3,]  3.9586 -0.13363 17.56600
##   Standard Deviations: 4.0336 3.9344 4.1912 
## > Gamma.un.10
##           [,1]         [,2]        [,3]         [,4]
## [1,]  1.0000000  0.258678035 -0.41810573  0.234156597
## [2,]  0.2586780  1.000000000 -0.01547345 -0.008103682
## [3,] -0.4181057 -0.015473451  1.00000000  0.727363531
## [4,]  0.2341566 -0.008103682  0.72736353  1.000000000
## > V.un.10
## Marginal variance covariance matrix
##          [,1]     [,2]      [,3]     [,4]
## [1,]  16.2700  4.10510 -14.19000  3.95860
## [2,]   4.1051 15.47900  -0.51224 -0.13363
## [3,] -14.1900 -0.51224  70.79800 25.65100
## [4,]   3.9586 -0.13363  25.65100 17.56600
##   Standard Deviations: 4.0336 3.9344 8.4142 4.1912 

#  Common compound symmetry, same variance all times

hips.cs <- gls(h ~ -1 + gender + week:gender + week2:gender,
               correlation=corCompSymm(form = ~ time | id),
               data=thedat,method="ML")
beta.cs <- coef(hips.cs)
sebeta.cs <- robust.cov(hips.cs,m)$se.model
Gamma.cs.1 <- corMatrix(hips.cs$modelStruct$corStruct)[[1]]   
Gamma.cs.10 <- corMatrix(hips.cs$modelStruct$corStruct)[[10]]   
V.cs.1 <- getVarCov(hips.cs, individual=1)
V.cs.10 <- getVarCov(hips.cs, individual=10)
vars.cs <- diag(V.cs.10)

## > rbind(beta.cs,sebeta.cs)
##             gender0   gender1 gender0:week gender1:week gender0:week2
## beta.cs   38.291744 42.185601    -9.673945   -14.245522     2.6078898
## sebeta.cs  1.014631  1.132319     1.608973     1.868087     0.5025518
##           gender1:week2
## beta.cs       3.8338897
## sebeta.cs     0.5874198
## > Gamma.cs.10
##          [,1]      [,2]      [,3]      [,4]
## [1,] 1.0000000 0.1945206 0.1945206 0.1945206
## [2,] 0.1945206 1.0000000 0.1945206 0.1945206
## [3,] 0.1945206 0.1945206 1.0000000 0.1945206
## [4,] 0.1945206 0.1945206 0.1945206 1.0000000
## > V.cs.10
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]
## [1,] 16.985  3.304  3.304  3.304
## [2,]  3.304 16.985  3.304  3.304
## [3,]  3.304  3.304 16.985  3.304
## [4,]  3.304  3.304  3.304 16.985
##   Standard Deviations: 4.1213 4.1213 4.1213 4.1213 

#  Common AR(1), same variance at all times

hips.ar1 <- gls(h ~ -1 + gender + week:gender + week2:gender,
               correlation=corAR1(form = ~ time | id),
               data=thedat,method="ML")
beta.ar1 <- coef(hips.ar1)
sebeta.ar1 <- robust.cov(hips.ar1,m)$se.model
Gamma.ar1.1 <- corMatrix(hips.ar1$modelStruct$corStruct)[[1]]   
Gamma.ar1.10 <- corMatrix(hips.ar1$modelStruct$corStruct)[[10]]   
V.ar1.1 <- getVarCov(hips.ar1, individual=1)
V.ar1.10 <- getVarCov(hips.ar1, individual=10)
vars.ar1 <- diag(V.ar1.10)

## > rbind(beta.ar1,sebeta.ar1)
##              gender0   gender1 gender0:week gender1:week gender0:week2
## beta.ar1   38.402047 42.331923    -9.877170   -14.574246     2.6482744
## sebeta.ar1  1.024547  1.140782     1.587993     1.825175     0.4954168
##            gender1:week2
## beta.ar1       3.9083660
## sebeta.ar1     0.5752724
## > Gamma.ar1.10
##            [,1]       [,2]       [,3]       [,4]
## [1,] 1.00000000 0.27368526 0.07490362 0.02050002
## [2,] 0.27368526 1.00000000 0.27368526 0.07490362
## [3,] 0.07490362 0.27368526 1.00000000 0.27368526
## [4,] 0.02050002 0.07490362 0.27368526 1.00000000
## > V.ar1.10
## Marginal variance covariance matrix
##          [,1]    [,2]    [,3]     [,4]
## [1,] 17.00100  4.6529  1.2734  0.34852
## [2,]  4.65290 17.0010  4.6529  1.27340
## [3,]  1.27340  4.6529 17.0010  4.65290
## [4,]  0.34852  1.2734  4.6529 17.00100
##   Standard Deviations: 4.1232 4.1232 4.1232 4.1232 

# One-dependent, same variance all times

hips.1d <- gls(h ~ -1 + gender + week:gender + week2:gender,
               correlation=corARMA(form = ~ time | id, q=1),
               data=thedat,method="ML")
beta.1d <- coef(hips.1d)
sebeta.1d <- robust.cov(hips.1d,m)$se.model
Gamma.1d.1 <- corMatrix(hips.1d$modelStruct$corStruct)[[1]]   
Gamma.1d.10 <- corMatrix(hips.1d$modelStruct$corStruct)[[10]]   
V.1d.1 <- getVarCov(hips.1d, individual=1)
V.1d.10 <- getVarCov(hips.1d, individual=10)
vars.1d <- diag(V.1d.10)
## > rbind(beta.1d,sebeta.1d)
##             gender0   gender1 gender0:week gender1:week gender0:week2
## beta.1d   38.532200 42.453961   -10.040464   -14.788859     2.6739058
## sebeta.1d  1.032581  1.148984     1.605236     1.833158     0.5054891
##           gender1:week2
## beta.1d       3.9513880
## sebeta.1d     0.5843124
## > Gamma.1d.10
##           [,1]      [,2]      [,3]      [,4]
## [1,] 1.0000000 0.3089592 0.0000000 0.0000000
## [2,] 0.3089592 1.0000000 0.3089592 0.0000000
## [3,] 0.0000000 0.3089592 1.0000000 0.3089592
## [4,] 0.0000000 0.0000000 0.3089592 1.0000000
## > V.1d.10
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]
## [1,] 17.1660  5.3035  0.0000  0.0000
## [2,]  5.3035 17.1660  5.3035  0.0000
## [3,]  0.0000  5.3035 17.1660  5.3035
## [4,]  0.0000  0.0000  5.3035 17.1660
##   Standard Deviations: 4.1431 4.1431 4.1431 4.1431 

#  Compare models via AIC and BIC -- one-dependent model is preferred

## > anova(hips.un,hips.cs,hips.ar1,hips.1d)
##          Model df      AIC      BIC    logLik   Test  L.Ratio p-value
## hips.un      1 16 574.4083 615.9302 -271.2042                        
## hips.cs      2  8 573.4241 594.1851 -278.7121 1 vs 2 15.01581  0.0588
## hips.ar1     3  8 573.3716 594.1325 -278.6858                        
## hips.1d      4  8 573.0199 593.7809 -278.5100           

u <- hips.1d

#  The model-based covariance matrix can be found as follows; however,
#  it is not correct

Sig.0 <- vcov(u)

#  we use the model=based covariance matrix computed by the function above
#  in the analysis below

Sig.model <- robust.cov(u,m)$Sig.model
sebeta.model <- robust.cov(u,m)$se.model
t.model <- beta.1d/sebeta.model
df <- nrow(thedat)-length(beta.1d)    #  this is the residual df used
                                      #  by gls; SAS uses a different approximation
p.value <- round(2*pt(-abs(t.model),df),4)

## > cbind(beta.1d,sebeta.model,t.model,p.value)
##                  beta.1d sebeta.model   t.model p.value
## gender0        38.532200    1.0325811 37.316390       0
## gender1        42.453961    1.1489845 36.949116       0
## gender0:week  -10.040464    1.6052360 -6.254821       0
## gender1:week  -14.788859    1.8331583 -8.067421       0
## gender0:week2   2.673906    0.5054891  5.289740       0
## gender1:week2   3.951388    0.5843124  6.762458       0

#  Fit the model incorporating age in the intercept

hips.1d.age <- gls(h ~ -1 + gender + gender:age + week:gender + week2:gender,
               correlation=corARMA(form = ~ time | id, q=1),
               data=thedat,method="ML")
beta.1d.age <- coef(hips.1d.age)
u <- hips.1d.age
Sig.model <- robust.cov(u,m)$Sig.model
sebeta.model <- robust.cov(u,m)$se.model
t.model <- beta.1d.age/sebeta.model
df <- nrow(thedat)-length(beta.1d.age)    
p.value <- round(2*pt(-abs(t.model),df),4)

## > cbind(beta.1d.age,sebeta.model,t.model,p.value)
##                beta.1d.age sebeta.model    t.model p.value
## gender0        36.19448323   4.17844132  8.6621973  0.0000
## gender1        40.29058620   6.13229657  6.5702279  0.0000
## gender0:age     0.03456706   0.06013804  0.5747952  0.5668
## gender1:age     0.03302359   0.09237520  0.3574941  0.7215
## gender0:week   -9.94833202   1.61076530 -6.1761524  0.0000
## gender1:week  -14.78348894   1.83411629 -8.0602789  0.0000
## gender0:week2   2.65394152   0.50614394  5.2434521  0.0000
## gender1:week2   3.95105019   0.58424891  6.7626146  0.0000

#  Function to contstruct F test statistics

f.test <- function(L,beta,Sig,dendf){
    numdf <- nrow(L)
    F <- t(L%*%beta)%*%solve(L%*%Sig%*%t(L))%*%(L%*%beta)/numdf
    p.value <- pf(F,numdf,dendf,lower.tail=FALSE)
    return(list(F=F,p.value=p.value,F.pvalue=c(F,round(p.value,4))))
}

#  Function to get estimate of form L beta and standard error
#  L has one row in this case; alpha is conf level for conf interval

est.L <- function(L,beta,Sig,alpha,df){
    est <- L%*%beta
    se <- sqrt(L%*%Sig%*%t(L))
    crit.val <- abs(qt(alpha/2,df))
    CI <- c(est-crit.val*se,est+crit.val*se)
    result <- matrix(c(est,se,CI),1,4,byrow=TRUE)
    colnames(result) <- c("Estimate","Std Error","Lower CI", "Upper CI")
    return(result)
}

#  Test diff in mean at 3 weeks between males and females 65 years old

L <- matrix(c(1,-1,65,-65,3,-3,9,-9),nrow=1,byrow=TRUE)

## > f.test(L,beta.1d.age,Sig.model,df)$F.pvalue
[1] 0.5596394 0.4563000

#  Estimate the mean for a 65 year old female at week 2.5

L <- matrix(c(1,0,65,0,2.5,0,6.25,0),nrow=1,byrow=TRUE)

## > est.L(L,beta.1d.age,Sig.model,0.05,df)
##     Estimate Std Error Lower CI Upper CI
## [1,] 30.15765  0.827914  28.5131  31.8022

#  Test if quadratic terms are needed at all

L <- matrix(c(0,0,0,0,0,0,1,0,
             0,0,0,0,0,0,0,1),nrow=2,byrow=TRUE)

## > f.test(L,beta.1d.age,Sig.model,df)$F.pvalue
## [1] 36.61337  0.00000


